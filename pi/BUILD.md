# Cross-building for the Pi (armhf)

`pi/Dockerfile.armhf` cross-builds `pi_tracker` and `solve_cli` for ARMv7
hardfloat (Pi 2 Model B and every newer Pi running the 32-bit OS) inside an
Alpine container. The host runs Docker with `qemu-arm-static` so no toolchain
or sysroot needs to be installed locally.

## Build

From the repo root:

```sh
docker build -f pi/Dockerfile.armhf -t pi_tracker_armhf .
```

Cold build is ~5 min on this machine. A BuildKit cache mount on `/src/build`
makes incremental rebuilds (after editing source) finish in ~30 s.

## Extract the binaries

```sh
id=$(docker create pi_tracker_armhf)
mkdir -p build-armhf
docker cp "$id:/out/pi_tracker" build-armhf/
docker cp "$id:/out/solve_cli"  build-armhf/
docker rm "$id"
```

Lands two ELF binaries in `build-armhf/`:

- `pi_tracker` — the daemon. Statically links libstdc++, libgcc, zlib;
  dynamically links libcamera (which `dlopen`s its IPA modules at runtime
  and so can't be static-linked cleanly).
- `solve_cli` — the plate solver. Fully static, single-file deploy.

## Deploy to the Pi

```sh
scp build-armhf/pi_tracker build-armhf/solve_cli pi@<host>:./
scp tetra3_db_stars.bin tetra3_db_patterns.bin pi@<host>:./
```

On the Pi, install the runtime deps once (Alpine ARMv7 with libcamera 0.7.x,
matching what the Dockerfile builds against):

```sh
apk add libcamera libcamera-raspberrypi
```

Then run:

```sh
./pi_tracker --solver ./solve_cli \
  --db-stars ./tetra3_db_stars.bin --db-patterns ./tetra3_db_patterns.bin \
  --gain 16 --exposure-us 1000000 --fov 75
```

## Live viewer

Pipe the daemon's NDJSON output through the terminal viewer to see the
SH1106-sized framebuffer + status:

```sh
ssh -tt <host> "./pi_tracker ..." | tools/star_catalog/view.py
```

`-tt` is needed so Ctrl-C forwards to the Pi process.

## libcamera SOVERSION pin

The Dockerfile tracks `arm32v7/alpine:edge` so the libcamera version it links
against matches what current Alpine on the Pi ships (0.7.x at time of
writing). If the Pi's Alpine release falls out of sync with edge, the binary
will fail to load — pin to a stable Alpine release once the SOVERSIONs align.
