#!/usr/bin/env python3
"""Full-sky solver coverage test.

Generates synthetic star fields on an even grid across the whole sky,
runs them through solve_cli in batch mode (DB loaded once per worker),
and reports solve rate, accuracy, and a sky map of failures.

Usage:
    cd tools/star_catalog
    .venv/bin/python coverage_test.py [--spacing 5] [--fov 75] [--jobs 4]
    .venv/bin/python coverage_test.py --spacing 1 --rolls 0,45,90 --jobs 14
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from multiprocessing import Process, Queue
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from devtool import load_hipparcos_cached, generate_star_field


def generate_grid(spacing_deg):
    """Generate evenly-spaced sky coordinates.

    Uses constant declination spacing with RA spacing adjusted by cos(dec)
    to get roughly uniform area coverage.
    """
    points = []
    dec = -90.0 + spacing_deg / 2.0
    while dec < 90.0:
        cos_dec = math.cos(math.radians(dec))
        ra_spacing = spacing_deg / cos_dec if cos_dec > 0.05 else 360.0
        ra = 0.0
        while ra < 360.0:
            points.append((ra, dec))
            ra += ra_spacing
        dec += spacing_deg
    return points


def worker(task_queue, result_queue, solver_path, db_stars, db_patterns, fov, rolls):
    """Worker process: keeps one solve_cli alive, generates images and pipes paths."""
    stars = load_hipparcos_cached()

    proc = subprocess.Popen(
        [solver_path,
         "--db-stars", db_stars,
         "--db-patterns", db_patterns,
         "--fov", str(fov),
         "--batch"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        bufsize=1,
    )

    tmp_dir = tempfile.mkdtemp(prefix="coverage_")

    try:
        while True:
            item = task_queue.get()
            if item is None:
                break  # poison pill

            ra, dec = item

            for roll in rolls:
                img = generate_star_field(stars, ra, dec, fov_deg=fov,
                                          roll_deg=roll, noise_level=0)

                tmp_path = os.path.join(tmp_dir, "img.pgm")
                with open(tmp_path, "wb") as f:
                    header = f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode()
                    f.write(header)
                    f.write(img.astype(">u2").tobytes())

                proc.stdin.write(tmp_path + "\n")
                proc.stdin.flush()
                line = proc.stdout.readline()

                try:
                    result = json.loads(line.strip()) if line.strip() else {"solved": False}
                except json.JSONDecodeError:
                    result = {"solved": False}

                err = None
                if result.get("solved"):
                    err_ra = (result["ra_deg"] - ra + 180) % 360 - 180
                    err_dec = result["dec_deg"] - dec
                    err = math.sqrt(
                        (err_ra * math.cos(math.radians(dec)))**2 + err_dec**2)
                    if err > 10:
                        err = None

                result_queue.put({
                    "ra": ra,
                    "dec": dec,
                    "roll": roll,
                    "solved": err is not None,
                    "error_deg": err,
                    "rmse": result.get("rmse", 0) if err is not None else None,
                    "time_ms": result.get("solve_time_ms", 0),
                    "num_detected": result.get("num_detected", 0),
                })
    finally:
        proc.stdin.close()
        proc.terminate()
        proc.wait(timeout=5)
        try:
            os.unlink(os.path.join(tmp_dir, "img.pgm"))
        except OSError:
            pass
        os.rmdir(tmp_dir)


def print_sky_map(results, spacing_deg):
    """Print an ASCII sky map showing solve success/failure."""
    dec_bins = int(180 / spacing_deg)
    ra_bins = int(360 / spacing_deg)

    col_scale = max(1, ra_bins // 72)
    row_scale = max(1, dec_bins // 36)
    map_cols = ra_bins // col_scale
    map_rows = dec_bins // row_scale

    grid = [[None for _ in range(map_cols)] for _ in range(map_rows)]

    for r in results:
        col = int(r["ra"] / 360.0 * map_cols) % map_cols
        row = int((r["dec"] + 90) / 180.0 * map_rows)
        row = min(row, map_rows - 1)
        if grid[row][col] is None:
            grid[row][col] = r["solved"]
        elif not r["solved"]:
            grid[row][col] = False

    print("\nSky map (· = solved, X = failed, - = no data):")
    print(f"  RA: 0° {'':>30s} 180° {'':>30s} 360°")
    for row_idx in range(map_rows - 1, -1, -1):
        dec_label = -90 + (row_idx + 0.5) * 180.0 / map_rows
        line = ""
        for cell in grid[row_idx]:
            if cell is None:
                line += "-"
            elif cell:
                line += "·"
            else:
                line += "X"
        print(f"{dec_label:+5.0f}° |{line}|")


def print_dec_bands(results, band_size=15):
    """Print solve rate by declination band."""
    bands = {}
    for r in results:
        band = int(r["dec"] / band_size) * band_size
        if band not in bands:
            bands[band] = {"solved": 0, "total": 0, "errors": [], "times": []}
        bands[band]["total"] += 1
        if r["solved"]:
            bands[band]["solved"] += 1
            bands[band]["errors"].append(r["error_deg"])
            bands[band]["times"].append(r["time_ms"])

    print(f"\n{'Dec band':<12s} {'Solved':>10s} {'Rate':>6s} {'Avg Err':>9s} {'Avg Time':>9s}")
    print("-" * 50)
    for band in sorted(bands.keys()):
        b = bands[band]
        pct = 100 * b["solved"] / b["total"]
        avg_err = np.mean(b["errors"]) if b["errors"] else float("nan")
        avg_time = np.mean(b["times"]) if b["times"] else float("nan")
        bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
        print(f"{band:+4d}° to {band+band_size:+4d}° {b['solved']:>4d}/{b['total']:<4d}"
              f" {pct:>5.1f}% {avg_err:>8.4f}° {avg_time:>7.0f}ms {bar}")


def print_roll_bands(results, rolls):
    """Print solve rate by roll angle."""
    if len(rolls) <= 1:
        return
    bands = {}
    for r in results:
        roll = r["roll"]
        if roll not in bands:
            bands[roll] = {"solved": 0, "total": 0, "errors": []}
        bands[roll]["total"] += 1
        if r["solved"]:
            bands[roll]["solved"] += 1
            bands[roll]["errors"].append(r["error_deg"])

    print(f"\n{'Roll':<8s} {'Solved':>10s} {'Rate':>6s} {'Avg Err':>9s}")
    print("-" * 40)
    for roll in sorted(bands.keys()):
        b = bands[roll]
        pct = 100 * b["solved"] / b["total"]
        avg_err = np.mean(b["errors"]) if b["errors"] else float("nan")
        print(f"{roll:+6.0f}°  {b['solved']:>4d}/{b['total']:<4d} {pct:>5.1f}% {avg_err:>8.4f}°")


def main():
    script_dir = Path(__file__).parent
    project_root = script_dir / ".." / ".."

    parser = argparse.ArgumentParser(description="Full-sky solver coverage test")
    parser.add_argument(
        "--solver",
        default=str(project_root / "build" / "pi" / "solve_cli"),
    )
    parser.add_argument(
        "--db-stars",
        default=str(project_root / "tetra3_db_stars.bin"),
    )
    parser.add_argument(
        "--db-patterns",
        default=str(project_root / "tetra3_db_patterns.bin"),
    )
    parser.add_argument("--fov", type=float, default=75.0)
    parser.add_argument("--spacing", type=float, default=5.0,
                        help="Grid spacing in degrees (default: 5)")
    parser.add_argument("--rolls", type=str, default="0",
                        help="Comma-separated roll angles to test (default: 0)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="Parallel worker processes")
    args = parser.parse_args()

    rolls = [float(r) for r in args.rolls.split(",")]
    points = generate_grid(args.spacing)
    total_solves = len(points) * len(rolls)

    print(f"Coverage test: {len(points)} fields × {len(rolls)} rolls = {total_solves} solves")
    print(f"Solver: {args.solver}")
    print(f"FOV: {args.fov}°, Rolls: {rolls}, Workers: {args.jobs}")
    print(f"Using batch mode (DB loaded once per worker)")
    print()

    task_queue = Queue()
    result_queue = Queue()

    # Start workers
    workers = []
    for _ in range(args.jobs):
        p = Process(target=worker, args=(
            task_queue, result_queue,
            args.solver, args.db_stars, args.db_patterns,
            args.fov, rolls,
        ))
        p.start()
        workers.append(p)

    # Feed tasks
    for ra, dec in points:
        task_queue.put((ra, dec))
    # Poison pills
    for _ in workers:
        task_queue.put(None)

    # Collect results
    results = []
    t0 = time.monotonic()
    while len(results) < total_solves:
        r = result_queue.get()
        results.append(r)
        done = len(results)
        if done % 200 == 0 or done == total_solves:
            elapsed = time.monotonic() - t0
            rate = done / elapsed
            eta = (total_solves - done) / rate if rate > 0 else 0
            solved_so_far = sum(1 for r in results if r["solved"])
            print(f"\r  {done}/{total_solves} ({100*done/total_solves:.0f}%) "
                  f"solved={solved_so_far} "
                  f"rate={rate:.1f}/s ETA={eta:.0f}s", end="", flush=True)

    for p in workers:
        p.join()

    print()
    elapsed = time.monotonic() - t0

    # Summary
    solved = sum(1 for r in results if r["solved"])
    errors = [r["error_deg"] for r in results if r["solved"]]
    times = [r["time_ms"] for r in results]

    print(f"\n{'='*60}")
    print(f"Results: {solved}/{len(results)} solved ({100*solved/len(results):.1f}%)")
    print(f"Time: {elapsed:.1f}s ({len(results)/elapsed:.1f} fields/s)")
    if errors:
        print(f"Position error: mean={np.mean(errors):.4f}° "
              f"median={np.median(errors):.4f}° "
              f"max={np.max(errors):.4f}° "
              f"p95={np.percentile(errors, 95):.4f}°")
        print(f"Solve time: mean={np.mean(times):.0f}ms "
              f"median={np.median(times):.0f}ms "
              f"max={np.max(times):.0f}ms")
    print(f"{'='*60}")

    print_dec_bands(results)
    print_roll_bands(results, rolls)
    print_sky_map(results, args.spacing)

    # Save raw results
    out_path = script_dir / "coverage_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nRaw results saved to {out_path}")


if __name__ == "__main__":
    main()
