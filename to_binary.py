#!/usr/bin/env python3
"""Convert a tetra3 .npz star catalog to the two .bin files the C++ solver
loads (tetra3_db_stars.bin, tetra3_db_patterns.bin).

tetra3 ships a `data/default_database.npz` in its source tree. Get the
upstream project and point this script at the file:

    git clone https://github.com/esa/tetra3 ../tetra3
    tools/star_catalog/.venv/bin/python to_binary.py ../tetra3/tetra3/data/default_database.npz

If you have tetra3 installed as a Python package, the bundled .npz is
located automatically; pass no argument.

Run from the repo root so the output files land where the solver expects.
"""

import argparse
import importlib.util
import struct
import sys
from pathlib import Path

import numpy as np


# Where to look for the bundled .npz when no path is given.
SEARCH_HINTS = [
    Path.home() / "git/tetra3/tetra3/data/default_database.npz",
    Path.home() / "git/OpenAstroExplorer/tracker/tetra3/tetra3/data/default_database.npz",
]


def find_tetra3_npz():
    """Locate default_database.npz: try the installed package, then a few
    common clone locations."""
    try:
        spec = importlib.util.find_spec("tetra3")
        if spec and spec.submodule_search_locations:
            for root in spec.submodule_search_locations:
                p = Path(root) / "data" / "default_database.npz"
                if p.exists():
                    return p
    except Exception:
        pass
    for hint in SEARCH_HINTS:
        if hint.exists():
            return hint
    return None


def convert(npz_path: Path, stars_out: Path, patterns_out: Path):
    data = np.load(npz_path)
    stars = data["star_table"].astype(np.float32)
    patterns = data["pattern_catalog"]   # uint16 (npatterns, 4)
    props = data["props_packed"]

    print(f"source:   {npz_path}", file=sys.stderr)
    print(f"  props:  {props}", file=sys.stderr)
    print(f"  stars:  {stars.shape} {stars.dtype}", file=sys.stderr)
    print(f"  pats:   {patterns.shape} {patterns.dtype}", file=sys.stderr)

    if patterns.dtype != np.uint16:
        raise SystemExit(
            f"pattern catalog must be uint16, got {patterns.dtype}")

    with open(stars_out, "wb") as f:
        f.write(struct.pack("I", stars.shape[0]))
        f.write(stars.tobytes())
    with open(patterns_out, "wb") as f:
        f.write(struct.pack("I", patterns.shape[0]))
        f.write(patterns.tobytes())

    print(f"wrote {stars_out} ({stars_out.stat().st_size:,} bytes)",
          file=sys.stderr)
    print(f"wrote {patterns_out} ({patterns_out.stat().st_size:,} bytes)",
          file=sys.stderr)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source", nargs="?", default=None,
                   help="Path to tetra3 .npz file. If omitted, search the "
                        "installed `tetra3` package + common clone paths.")
    p.add_argument("--stars-out",    default="tetra3_db_stars.bin")
    p.add_argument("--patterns-out", default="tetra3_db_patterns.bin")
    args = p.parse_args()

    if args.source is None:
        path = find_tetra3_npz()
        if path is None:
            print("ERROR: pass a .npz path, or clone tetra3 to one of:",
                  file=sys.stderr)
            for hint in SEARCH_HINTS:
                print(f"  {hint}", file=sys.stderr)
            sys.exit(1)
    else:
        path = Path(args.source)
        if not path.exists():
            print(f"ERROR: {path} not found.", file=sys.stderr)
            sys.exit(1)

    convert(path, Path(args.stars_out), Path(args.patterns_out))


if __name__ == "__main__":
    main()
