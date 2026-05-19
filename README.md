# tetra3-cpp

A C++ port of [ESA's tetra3](https://github.com/esa/tetra3) plate solver, plus
a small standalone star-tracker that runs on a Raspberry Pi Zero (and an
earlier ESP32 prototype).

The motivation: do polar alignment in the field without a laptop. Plate-solve
a few frames as the mount rotates around its mechanical axis, fit a circle on
the celestial sphere, and read the alt/az offset directly off an OLED.

## What's here

| Path                | What it is                                                            |
|---------------------|-----------------------------------------------------------------------|
| `lib/solver/`       | Core plate-solver. Loads tetra3's binary catalog, matches centroids.  |
| `pi/`               | Pi Zero daemon (`pi_tracker`) + `solve_cli` + `oled_test`.            |
| `esp/`              | ESP32-CAM prototype. Earlier exploration — see notes below.           |
| `tools/star_catalog/` | Python tooling: catalog gen, viewers, accuracy/coverage tests.      |
| `tetra3_db_*.bin`   | Pre-built pattern + star database (mag-7, 4-star patterns, ~10° FOV). |

## Hardware

### Pi Zero 2 W + Camera Module 3 NoIR

The intended target. The Module 3 sees enough stars at 100 ms exposures to
solve reliably with a 72° FOV. Solve RMSE is around 12 arcseconds in field
tests — vs. ~100″ for astrometry.net's `solve-field` on the same images at
default settings. See [`pi/BUILD.md`](pi/BUILD.md) for cross-compilation.

### ESP32-CAM (earlier prototype)

The southern-sky catalog fits in PSRAM and solves run in a couple of seconds.
The OV-series camera modules turned out to be the limit — limited control,
soft lenses. Kept in `esp/` for reference; not actively developed.

## Build

The intended build path is the Docker cross-compile in
[`pi/BUILD.md`](pi/BUILD.md) — it produces statically-linked `pi_tracker` and
`solve_cli` binaries for ARMv7 hardfloat (Pi 2 and newer running the 32-bit
OS).

For development on x86 you can build the solver library and tests natively:

```sh
cmake -S . -B build
cmake --build build --target tetra3_solver test_polar_align oled_test -j
```

The full `pi_tracker` / `solve_cli` native build needs `stb_image.h`
(provided by Alpine's `stb` package inside the Docker image; on Debian/Ubuntu
fetch it manually from <https://github.com/nothings/stb>).

## Run

Replay a saved frame through the daemon:

```sh
./build/pi/pi_tracker \
  --image some_frame.jpg \
  --db-stars tetra3_db_stars.bin \
  --db-patterns tetra3_db_patterns.bin \
  --solver ./build/pi/solve_cli \
  --fov 75
```

Output is NDJSON on stdout: solve result, framebuffer (base64), timing, mode.
The terminal viewer renders it:

```sh
./build/pi/pi_tracker ... | tools/star_catalog/view.py
```

## Catalogs

The solver loads two binary files at runtime:

- `tetra3_db_stars.bin` — flat array of `(ra, dec, x, y, z, mag)` floats.
- `tetra3_db_patterns.bin` — hash table of 4-star pattern indices into
  the star table, ~100 MB at mag-7. Too large to keep in git; build it
  locally with:

```sh
git clone https://github.com/esa/tetra3 ../tetra3   # one-time
tools/star_catalog/.venv/bin/python to_binary.py \
    ../tetra3/tetra3/data/default_database.npz
```

`to_binary.py` with no argument also auto-locates the bundled
`default_database.npz` from common clone paths or an installed
`tetra3` package. The result is the same catalog tetra3 itself uses —
8818 stars to magnitude 7, ~12 M 4-star patterns for FOVs in 10°–30°.

`generate_catalog_headers.py` builds smaller catalogs flashed directly
into the ESP32 firmware (`esp/`).

## Tests

```sh
cmake --build build --target test_polar_align
./build/pi/test_polar_align
```

Covers the polar-alignment pipeline (`pi/src/app/polar_align.cpp`) with
synthetic samples around a known fake pole.

## License

Apache 2.0. tetra3 itself is Apache-2.0; this port inherits.
