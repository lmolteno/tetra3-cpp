#!/usr/bin/env python3
"""Solver accuracy test.

Generates synthetic star fields at various coordinates, runs them through
solve_cli, and reports solve rate and positional accuracy.

Usage:
    .venv/bin/python accuracy_test.py [--solver PATH] [--trials N]
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from devtool import load_hipparcos_cached, generate_star_field

# Test coordinates: (ra_deg, dec_deg, label)
TEST_FIELDS = [
    (186.0, -63.0, "Southern Cross"),
    (84.05, -1.2, "Orion"),
    (280.0, -26.0, "Sagittarius"),
    (0.0, -89.0, "South Pole"),
    (300.0, 45.0, "Cygnus"),
    (80.0, 46.0, "Auriga"),
    (201.0, -43.0, "Centaurus"),
    (247.0, -26.0, "Scorpius"),
    (23.0, 15.0, "Aries"),
    (170.0, 12.0, "Leo"),
]


def run_test(solver_path, db_stars, db_patterns, fov, trials, fields):
    stars = load_hipparcos_cached()
    print(f"Loaded {len(stars)} Hipparcos stars")
    print(f"Solver: {solver_path}")
    print(f"FOV: {fov}°, Trials per field: {trials}")
    print()
    print(f"{'Field':<20s} {'Solved':>6s} {'Avg Err':>9s} {'Max Err':>9s} {'Avg RMSE':>9s} {'Avg Time':>9s}")
    print("-" * 72)

    total_solved = 0
    total_trials = 0
    all_errors = []

    for ra, dec, label in fields:
        errors = []
        rmses = []
        times = []

        for _ in range(trials):
            img = generate_star_field(stars, ra, dec, fov_deg=fov)

            with tempfile.NamedTemporaryFile(suffix=".pgm", delete=False) as f:
                tmp_path = f.name
                header = f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode()
                f.write(header)
                f.write(img.astype(">u2").tobytes())

            try:
                proc = subprocess.run(
                    [solver_path,
                     "--db-stars", db_stars,
                     "--db-patterns", db_patterns,
                     "--fov", str(fov),
                     tmp_path],
                    capture_output=True, text=True, timeout=30,
                )
                lines = [l for l in proc.stdout.strip().splitlines() if l.strip()]
                result = json.loads(lines[-1]) if lines else {"solved": False}
            except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError) as e:
                result = {"solved": False}
            finally:
                os.unlink(tmp_path)

            if result.get("solved"):
                err_ra = (result["ra_deg"] - ra + 180) % 360 - 180
                err_dec = result["dec_deg"] - dec
                # Correct RA difference for cos(dec) to get true angular error
                err = np.sqrt((err_ra * np.cos(np.radians(dec)))**2 + err_dec**2)
                # Reject obvious false matches
                if err < 10:
                    errors.append(err)
                    rmses.append(result.get("rmse", 0))
                    times.append(result.get("solve_time_ms", 0))
                else:
                    errors.append(None)
            else:
                errors.append(None)

        solved = sum(1 for e in errors if e is not None)
        total_solved += solved
        total_trials += trials
        good = [e for e in errors if e is not None]
        all_errors.extend(good)

        avg_err = np.mean(good) if good else float("nan")
        max_err = np.max(good) if good else float("nan")
        avg_rmse = np.mean(rmses) if rmses else float("nan")
        avg_time = np.mean(times) if times else float("nan")

        status = f"{solved}/{trials}"
        print(f"{label:<20s} {status:>6s} {avg_err:>8.4f}° {max_err:>8.4f}° {avg_rmse:>7.1f}\" {avg_time:>7.1f}ms")

    print("-" * 72)
    pct = 100 * total_solved / total_trials if total_trials else 0
    avg_all = np.mean(all_errors) if all_errors else float("nan")
    print(f"{'TOTAL':<20s} {total_solved}/{total_trials} ({pct:.0f}%)  avg err={avg_all:.4f}°")


def main():
    script_dir = Path(__file__).parent
    project_root = script_dir / ".." / ".."

    parser = argparse.ArgumentParser(description="Solver accuracy test")
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
    parser.add_argument("--trials", type=int, default=5)
    args = parser.parse_args()

    run_test(args.solver, args.db_stars, args.db_patterns, args.fov, args.trials, TEST_FIELDS)


if __name__ == "__main__":
    main()
