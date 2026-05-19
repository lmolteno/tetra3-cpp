#!/usr/bin/env python3
"""End-to-end polar-alignment accuracy test.

Drives a real pi_tracker subprocess through a synthetic PA workflow:

  1. Pick a mount pole misaligned from the true celestial pole by a known
     amount in a randomised direction.
  2. Synthesize a sequence of camera frames as the operator rotates the
     mount around its (misaligned) RA axis — every frame is a noisy
     star-field image centred on a different point of the small circle
     traced by the camera.
  3. Feed those frames into pi_tracker via the --image-watch mechanism.
  4. Step the daemon: Tracking -> PASampling (collect arc) -> PAFix.
  5. Read the reported alt/az offset and compare to ground truth.

Sweeps noise level and prints accuracy statistics. Use --trials and
--noise-levels to control sweep size.

Usage:
    .venv/bin/python pa_accuracy.py [--trials 5] [--samples 8]
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from imaging import (generate_star_field, load_hipparcos_cached,
                     load_binary_catalog, save_pgm_16)


# ---------------------------------------------------------------------------
# Geometry helpers

def unit(ra_rad, dec_rad):
    cd = np.cos(dec_rad)
    return np.array([cd * np.cos(ra_rad), cd * np.sin(ra_rad), np.sin(dec_rad)])


def radec(v):
    dec = np.arcsin(np.clip(v[2], -1.0, 1.0))
    ra = np.arctan2(v[1], v[0]) % (2 * np.pi)
    return float(ra), float(dec)


def rotate_axis(v, axis, theta):
    """Rotate vector v around unit axis by theta (Rodrigues)."""
    c, s = np.cos(theta), np.sin(theta)
    return v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1.0 - c)


def perp_unit(v):
    """A unit vector perpendicular to v (any direction)."""
    if abs(v[2]) < 0.99:
        e = np.cross(v, np.array([0.0, 0.0, 1.0]))
    else:
        e = np.cross(v, np.array([1.0, 0.0, 0.0]))
    return e / np.linalg.norm(e)


def pick_mount_pole(true_pole_ra, true_pole_dec, offset_arcmin, az_dir_rad):
    """Pick a mount pole offset_arcmin from the true pole, in az_dir direction."""
    pole = unit(true_pole_ra, true_pole_dec)
    e1 = perp_unit(pole)
    e2 = np.cross(pole, e1)
    offset_rad = np.radians(offset_arcmin / 60.0)
    perp = np.cos(az_dir_rad) * e1 + np.sin(az_dir_rad) * e2
    v = pole * np.cos(offset_rad) + perp * np.sin(offset_rad)
    return radec(v)


def pa_arc_pointings(mount_pole_ra, mount_pole_dec, target_offset_deg,
                     arc_deg, n_samples):
    """Camera pointings as the mount rotates by `arc_deg` around its pole.

    Camera is at target_offset_deg from the mount pole (so it traces a
    small circle of that angular radius).
    """
    pole = unit(mount_pole_ra, mount_pole_dec)
    e = perp_unit(pole)
    init = pole * np.cos(np.radians(target_offset_deg)) + \
           e * np.sin(np.radians(target_offset_deg))
    arc_rad = np.radians(arc_deg)
    pts = []
    for i in range(n_samples):
        frac = 0.0 if n_samples == 1 else i / (n_samples - 1)
        cam = rotate_axis(init, pole, frac * arc_rad)
        pts.append(radec(cam))
    return pts


# ---------------------------------------------------------------------------
# Image IO

save_pgm = save_pgm_16


# ---------------------------------------------------------------------------
# pi_tracker driver

class Tracker:
    """Manages a running pi_tracker subprocess."""

    def __init__(self, tracker, solver, db_stars, db_patterns,
                 image_path, width, height, fov, frames_dir, state_dir):
        env = os.environ.copy()
        env["LIBCAMERA_LOG_LEVELS"] = "*:ERROR"
        self.proc = subprocess.Popen(
            [str(tracker),
             "--image", str(image_path), "--image-watch",
             "--solver", str(solver),
             "--db-stars", str(db_stars), "--db-patterns", str(db_patterns),
             "--width", str(width), "--height", str(height),
             "--fov", str(fov),
             "--max-fps", "30",
             "--frames-dir", str(frames_dir),
             "--state-dir", str(state_dir),
             "--no-distortion-k"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True, env=env,
        )

    def send(self, cmd):
        self.proc.stdin.write(json.dumps(cmd) + "\n")
        self.proc.stdin.flush()

    def read_tick(self, timeout=2.0):
        """Read one NDJSON line from the tracker. Returns dict or None."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            line = self.proc.stdout.readline()
            if not line:
                return None
            line = line.strip()
            if not line:
                continue
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                continue
        return None

    def wait_for_solve(self, fresh_after_ms=150, timeout=4.0):
        """Wait for a tick with a fresh solve (result_age below threshold)."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            tick = self.read_tick(timeout=max(0.1, deadline - time.monotonic()))
            if tick is None:
                continue
            if tick.get("solved") and (tick.get("timing_ms", {})
                                       .get("result_age", 1e9) < fresh_after_ms):
                return tick
        return None

    def wait_for_mode(self, mode, timeout=2.0):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            tick = self.read_tick(timeout=max(0.1, deadline - time.monotonic()))
            if tick is None:
                continue
            if tick.get("mode") == mode:
                return tick
        return None

    def drain(self, count=20, timeout=0.5):
        """Read up to `count` ticks and discard."""
        deadline = time.monotonic() + timeout
        for _ in range(count):
            if time.monotonic() >= deadline:
                break
            self.read_tick(timeout=0.1)

    def close(self):
        try:
            self.proc.stdin.close()
        except Exception:
            pass
        try:
            self.proc.terminate()
            self.proc.wait(timeout=2)
        except Exception:
            try:
                self.proc.kill()
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Trial

def run_trial(t, stars, args, noise_level, mount_offset_arcmin, az_dir_rad):
    """One PA trial: clear samples, sweep arc, enter PAFix, read offset.

    Returns dict with reported_total, reported_alt, reported_az, num_solves,
    or None on failure.
    """
    true_pole_dec = np.pi / 2 if args.northern else -np.pi / 2
    mount_ra, mount_dec = pick_mount_pole(
        0.0, true_pole_dec, mount_offset_arcmin, az_dir_rad)

    # Camera traces a small circle at args.target_offset° from the mount pole,
    # sweeping args.arc° in args.samples steps.
    pointings = pa_arc_pointings(
        mount_ra, mount_dec, args.target_offset, args.arc, args.samples)

    # Write the first image so the tracker has something to load.
    first_ra_deg = np.degrees(pointings[0][0])
    first_dec_deg = np.degrees(pointings[0][1])
    img = generate_star_field(
        stars, first_ra_deg, first_dec_deg,
        fov_deg=args.fov, width=args.width, height=args.height,
        noise_level=noise_level, psf_sigma=args.psf_sigma,
        flux_scale=args.flux_scale)
    save_pgm(img, args.image_path)

    # Reset tracker into PASampling (clears samples on entry).
    t.send({"mode": "tracking"})
    time.sleep(0.05)
    t.drain(count=10, timeout=0.2)
    t.send({"mode": "pa_sampling"})

    num_solved = 0
    for i, (ra, dec) in enumerate(pointings):
        if i > 0:
            img = generate_star_field(
                stars, np.degrees(ra), np.degrees(dec),
                fov_deg=args.fov, width=args.width, height=args.height,
                noise_level=noise_level, psf_sigma=args.psf_sigma,
                flux_scale=args.flux_scale)
            save_pgm(img, args.image_path)
        tick = t.wait_for_solve(timeout=5.0)
        if tick is None:
            return None
        if tick.get("solved"):
            num_solved += 1

    # Enter PAFix and read the offset.
    t.send({"mode": "pa_fix"})
    # First PAFix tick may not have live pole yet; wait for one with alt/az.
    deadline = time.monotonic() + 3.0
    last_pa = None
    while time.monotonic() < deadline:
        tick = t.read_tick(timeout=0.3)
        if tick is None:
            continue
        if tick.get("mode") != "pa_fix":
            continue
        pa = tick.get("pa", {})
        if "total_arcmin" in pa:
            last_pa = pa
            break
    if last_pa is None:
        return None

    return {
        "num_solved": num_solved,
        "reported_total": last_pa["total_arcmin"],
        "reported_alt":   last_pa["alt_arcmin"],
        "reported_az":    last_pa["az_arcmin"],
        "pole_ra":        last_pa.get("pole_ra_rad"),
        "pole_dec":       last_pa.get("pole_dec_rad"),
        "arc_deg":        last_pa.get("arc_deg"),
        "true_mount_ra":  mount_ra,
        "true_mount_dec": mount_dec,
    }


# ---------------------------------------------------------------------------
# Reporting

def fmt_stats(label, values):
    if not values:
        return f"  {label}: (no data)"
    arr = np.array(values)
    return (f"  {label}: n={len(arr):2d}  "
            f"mean={arr.mean():6.2f}'  "
            f"median={np.median(arr):6.2f}'  "
            f"std={arr.std():5.2f}'  "
            f"min={arr.min():5.2f}'  max={arr.max():6.2f}'")


def main():
    project_root = Path(__file__).resolve().parent / ".." / ".."

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tracker", default=str(project_root / "build/pi/pi_tracker"))
    p.add_argument("--solver",  default=str(project_root / "build/pi/solve_cli"))
    p.add_argument("--db-stars",    default=str(project_root / "tetra3_db_stars.bin"))
    p.add_argument("--db-patterns", default=str(project_root / "tetra3_db_patterns.bin"))
    p.add_argument("--width",  type=int, default=2028)
    p.add_argument("--height", type=int, default=1520)
    p.add_argument("--fov",    type=float, default=72.0)
    p.add_argument("--psf-sigma", type=float, default=1.0,
                   help="PSF sigma in pixels; 1.0 matches Pi Cam Module 3")
    p.add_argument("--flux-scale", type=float, default=500000,
                   help="Flux scaling for synthetic stars")
    p.add_argument("--samples", type=int, default=8,
                   help="Number of PA samples per trial")
    p.add_argument("--arc", type=float, default=90.0,
                   help="Mount rotation arc (degrees)")
    p.add_argument("--target-offset", type=float, default=5.0,
                   help="Camera offset from mount pole (degrees)")
    p.add_argument("--trials", type=int, default=3,
                   help="Trials per noise level")
    p.add_argument("--noise-levels", type=str, default="4,8,16,32",
                   help="Comma-separated noise levels to sweep "
                        "(noise=0 makes the detector's local-BG std collapse "
                        "and it fails to find catalog stars)")
    p.add_argument("--star-source", choices=("hipparcos", "binary"),
                   default="hipparcos",
                   help="hipparcos: full mag<5 set (more realistic, like a "
                        "real camera image). binary: just the solver's "
                        "catalog stars (easier to solve, fewer false dets)")
    p.add_argument("--offset-arcmin", type=float, default=30.0,
                   help="Ground-truth mount misalignment magnitude")
    p.add_argument("--northern", action="store_true",
                   help="Use NCP instead of SCP")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    noise_levels = [float(x) for x in args.noise_levels.split(",")]

    if args.star_source == "binary":
        print(f"Loading binary catalog from {args.db_stars}...", file=sys.stderr)
        stars = load_binary_catalog(args.db_stars)
    else:
        print("Loading Hipparcos catalogue...", file=sys.stderr)
        stars = load_hipparcos_cached()
    print(f"Loaded {len(stars)} stars", file=sys.stderr)

    print(f"\nGround truth: mount pole {args.offset_arcmin}' from "
          f"{'NCP' if args.northern else 'SCP'}")
    print(f"Sweep: {args.trials} trials × {len(noise_levels)} noise levels"
          f" × {args.samples} samples each")
    print(f"Arc: {args.arc}°, target offset: {args.target_offset}°,"
          f" FOV: {args.fov}°")
    print()

    # One tracker process for the whole sweep. Reset via mode commands per trial.
    import tempfile
    with tempfile.TemporaryDirectory() as tmpd:
        args.image_path = Path(tmpd) / "sim.pgm"
        # Seed the image so the tracker can start
        seed_img = generate_star_field(
            stars, 0.0, -85.0 if not args.northern else 85.0,
            fov_deg=args.fov, width=args.width, height=args.height,
            noise_level=0)
        save_pgm(seed_img, args.image_path)

        t = Tracker(
            tracker=args.tracker, solver=args.solver,
            db_stars=args.db_stars, db_patterns=args.db_patterns,
            image_path=args.image_path,
            width=args.width, height=args.height, fov=args.fov,
            frames_dir=Path(tmpd) / "frames",
            state_dir=Path(tmpd) / "state",
        )

        # Wait for the daemon to come up and produce its first tick
        t.read_tick(timeout=5.0)

        for noise in noise_levels:
            errors = []         # |reported - true|, in arcmin
            reported = []       # raw reported total_arcmin
            for trial in range(args.trials):
                az = rng.uniform(0, 2 * np.pi)
                result = run_trial(
                    t, stars, args, noise, args.offset_arcmin, az)
                if result is None:
                    print(f"  noise={noise:5.1f}  trial {trial+1}: FAILED",
                          file=sys.stderr)
                    continue
                err = abs(result["reported_total"] - args.offset_arcmin)
                errors.append(err)
                reported.append(result["reported_total"])
                print(f"  noise={noise:5.1f}  trial {trial+1}/{args.trials}: "
                      f"reported={result['reported_total']:6.2f}'  "
                      f"|err|={err:5.2f}'  solves={result['num_solved']}/"
                      f"{args.samples}", file=sys.stderr)

            print(f"noise={noise:5.1f}")
            print(fmt_stats("reported", reported))
            print(fmt_stats("|error|",  errors))
            print()

        t.close()


if __name__ == "__main__":
    main()
