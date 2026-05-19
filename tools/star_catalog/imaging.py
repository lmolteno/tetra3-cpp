"""Synthetic star-field image generation.

Pulled out of devtool.py so other tools (pa_accuracy, pa_simulator) can
import the imaging helpers without dragging in the IMU server / aiohttp /
textual dependencies.

Two star sources are supported:

  load_hipparcos_cached()  — fetches Hipparcos up to a given magnitude from
                             Vizier (cached to .cache/). Big superset.
  load_binary_catalog()    — reads tetra3_db_stars.bin (the catalog the C++
                             solver actually loads). Subset of Hipparcos.

For tests that drive solve_cli, use load_binary_catalog: pattern matching
requires the synthetic stars to be in the *binary* star table since the
pattern catalog references them by index there.
"""

import struct
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Hipparcos star cache
# ---------------------------------------------------------------------------
CACHE_DIR = Path(__file__).parent / ".cache"
CACHE_FILE = CACHE_DIR / "hipparcos_5.0.npz"


def load_hipparcos_cached(mag_limit=5.0):
    """Load Hipparcos catalogue, caching to disk."""
    if CACHE_FILE.exists():
        data = np.load(CACHE_FILE)
        stars = {}
        for i in range(len(data["hip"])):
            stars[int(data["hip"][i])] = (
                float(data["ra"][i]),
                float(data["dec"][i]),
                float(data["mag"][i]),
            )
        return stars

    # Fetch from Vizier
    from astroquery.vizier import Vizier

    print(f"Fetching Hipparcos (mag < {mag_limit}) from Vizier (first time only)...")
    v = Vizier(columns=["HIP", "RAICRS", "DEICRS", "Vmag"], row_limit=-1)
    result = v.query_constraints(catalog="I/239/hip_main", Vmag=f"<{mag_limit}")
    if not result:
        print("ERROR: Vizier query failed", file=sys.stderr)
        sys.exit(1)

    hip_ids, ras, decs, mags = [], [], [], []
    for row in result[0]:
        try:
            h = int(row["HIP"])
            ra = float(row["RAICRS"])
            dec = float(row["DEICRS"])
            mag = float(row["Vmag"])
            if not (np.isnan(ra) or np.isnan(dec) or np.isnan(mag)):
                hip_ids.append(h)
                ras.append(np.radians(ra))
                decs.append(np.radians(dec))
                mags.append(mag)
        except (ValueError, KeyError):
            continue

    CACHE_DIR.mkdir(exist_ok=True)
    np.savez(
        CACHE_FILE,
        hip=np.array(hip_ids, dtype=np.int32),
        ra=np.array(ras, dtype=np.float64),
        dec=np.array(decs, dtype=np.float64),
        mag=np.array(mags, dtype=np.float64),
    )
    print(f"Cached {len(hip_ids)} stars to {CACHE_FILE}")

    stars = {}
    for i in range(len(hip_ids)):
        stars[hip_ids[i]] = (ras[i], decs[i], mags[i])
    return stars


def load_binary_catalog(path):
    """Read tetra3_db_stars.bin into a {index: (ra_rad, dec_rad, mag)} dict.

    The binary catalog is the *operational* star set the C++ solver uses.
    Pattern matching only works against these stars; synthesizing from a
    larger superset (e.g. raw Hipparcos) puts most "stars" in the image
    that the solver can't recognise — pattern hashes match nothing.
    """
    with open(path, "rb") as f:
        n = struct.unpack("I", f.read(4))[0]
        table = np.frombuffer(f.read(n * 24), dtype=np.float32).reshape(n, 6)
    # Columns: ra_rad, dec_rad, x, y, z, mag
    return {i: (float(table[i, 0]), float(table[i, 1]), float(table[i, 5]))
            for i in range(n)}


# ---------------------------------------------------------------------------
# Synthetic star field image generation
# ---------------------------------------------------------------------------

def _prepare_star_arrays(stars_dict):
    """Convert star dict to contiguous arrays for vectorized operations."""
    n = len(stars_dict)
    ra = np.empty(n, dtype=np.float64)
    dec = np.empty(n, dtype=np.float64)
    mag = np.empty(n, dtype=np.float64)
    for i, (r, d, m) in enumerate(stars_dict.values()):
        ra[i] = r
        dec[i] = d
        mag[i] = m
    return ra, dec, mag


_star_arrays_cache = {}


def _get_star_arrays(stars_dict):
    key = id(stars_dict)
    if key not in _star_arrays_cache:
        _star_arrays_cache[key] = _prepare_star_arrays(stars_dict)
    return _star_arrays_cache[key]


def generate_star_field(
    stars_dict,
    center_ra_deg,
    center_dec_deg,
    fov_deg=75.0,
    roll_deg=0.0,
    width=2028,
    height=1520,
    noise_level=8,
    bg_level=50,
    psf_sigma=2.0,
    flux_scale=200000.0,
):
    """Generate a 16-bit synthetic star field image.

    Uses gnomonic projection, Gaussian PSFs, and Gaussian read noise.
    roll_deg rotates the camera around the optical axis (degrees).
    Returns a uint16 numpy array of shape (height, width).
    """
    center_ra = np.radians(center_ra_deg)
    center_dec = np.radians(center_dec_deg)
    fov_rad = np.radians(fov_deg)
    roll_rad = np.radians(roll_deg)

    scale = width / (2.0 * np.tan(fov_rad / 2.0))
    cx = width / 2.0
    cy = height / 2.0
    cos_roll = np.cos(roll_rad)
    sin_roll = np.sin(roll_rad)
    sin_d0 = np.sin(center_dec)
    cos_d0 = np.cos(center_dec)

    all_ra, all_dec, all_mag = _get_star_arrays(stars_dict)

    sin_d = np.sin(all_dec)
    cos_d = np.cos(all_dec)
    dra = all_ra - center_ra
    cos_dra = np.cos(dra)
    sin_dra = np.sin(dra)

    cos_c = sin_d0 * sin_d + cos_d0 * cos_d * cos_dra
    visible = cos_c > 0.01
    idx = np.where(visible)[0]

    cos_c_v = cos_c[idx]
    x_proj = (cos_d[idx] * sin_dra[idx]) / cos_c_v
    y_proj = (cos_d0 * sin_d[idx] - sin_d0 * cos_d[idx] * cos_dra[idx]) / cos_c_v

    if roll_rad != 0.0:
        rx = x_proj * cos_roll - y_proj * sin_roll
        ry = x_proj * sin_roll + y_proj * cos_roll
        x_proj = rx
        y_proj = ry

    px = cx - x_proj * scale
    py = cy - y_proj * scale

    on_sensor = (px > -20) & (px < width + 20) & (py > -20) & (py < height + 20)
    sel = np.where(on_sensor)[0]
    px = px[sel]
    py = py[sel]
    mags = all_mag[idx[sel]]

    flux = np.power(10, -0.4 * mags) * flux_scale

    # Stamp each star with a Gaussian PSF
    sigma = float(psf_sigma)
    r = int(4 * sigma)
    norm = 1.0 / (2 * np.pi * sigma**2)
    two_sigma2 = 2.0 * sigma * sigma

    img = np.full((height, width), bg_level, dtype=np.float32)

    for i in range(len(px)):
        sx, sy, sf = px[i], py[i], flux[i]
        x0 = max(0, int(sx) - r)
        x1 = min(width, int(sx) + r + 1)
        y0 = max(0, int(sy) - r)
        y1 = min(height, int(sy) + r + 1)
        if x0 >= x1 or y0 >= y1:
            continue

        xs = np.arange(x0, x1, dtype=np.float32) - sx
        ys = np.arange(y0, y1, dtype=np.float32) - sy
        gx = np.exp(-xs * xs / two_sigma2)
        gy = np.exp(-ys * ys / two_sigma2)
        img[y0:y1, x0:x1] += (sf * norm) * np.outer(gy, gx)

    if noise_level > 0:
        img += np.random.normal(0, noise_level, img.shape).astype(np.float32)

    np.clip(img, 0, 65535, out=img)
    return img.astype(np.uint16)


def save_pgm_16(img, path):
    """Atomic write of a 16-bit PGM file."""
    import os
    tmp = str(path) + ".tmp"
    with open(tmp, "wb") as f:
        f.write(f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode())
        f.write(img.astype(">u2").tobytes())
    os.replace(tmp, path)
