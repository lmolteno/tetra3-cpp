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


# ---------------------------------------------------------------------------
# Pi Camera Module 3 NoIR (IMX708) calibration constants
#
# Empirical / approximate fit so synthetic images "look like" what the daemon
# sees with `--gain 16 --exposure-us 100000`. Numbers tuned so:
#
#   - V≈0 stars saturate (e.g. Sirius, Vega).
#   - V≈5 stars have peak ~1500 ADU above a ~1000 ADU light-polluted sky.
#   - V≈7 stars peak ~200 ADU above bg, near the detection floor.
#   - PSF is sharp (sigma ≈ 1 px) — Pi Cam at infinity focus.
#
# Cool stars (M/K class) read brighter on a NoIR sensor than V-mag suggests
# because the IR cut filter is missing — that effect is not modelled here
# (catalog has no B−V data), so this is an idealised V-band response.
PI_CAM_BASE = dict(
    bg_at_gain1_per_sec     = 400.0,   # ADU/s/pixel of light-polluted sky
    read_noise_at_gain1     = 2.0,     # ADU stddev at unity gain
    electrons_per_adu_at_g1 = 8.0,     # rough sensor gain at unity (e/ADU)
    flux_v0_at_gain1_per_s  = 240000.0,  # total ADU/s from a V=0 star
    psf_sigma_px            = 1.0,     # at infinity focus
)
# At gain=16, exposure=0.1s, those numbers give:
#   - sky bg            ≈ 640 ADU
#   - read noise stddev ≈ 32 ADU
#   - V=0 total flux    ≈ 384000 ADU → peak ≈ 61000 ADU (just saturating)
#   - V=5 total flux    ≈ 3840 ADU   → peak ≈ 611 ADU above bg
#   - V=7 total flux    ≈ 605 ADU    → peak ≈ 96 ADU above bg (~3σ vs noise)


def generate_star_field(
    stars_dict,
    center_ra_deg,
    center_dec_deg,
    fov_deg=75.0,
    roll_deg=0.0,
    width=2028,
    height=1520,
    *,
    # New Pi-Cam-like model (default). Pass exposure_ms / gain instead of
    # bg_level / noise_level / flux_scale.
    exposure_ms=100.0,
    gain=16.0,
    saturation_adu=65535,
    # Legacy knobs. If `noise_level` is given, fall back to the old simple
    # Gaussian-noise model so existing tests keep working unchanged.
    noise_level=None,
    bg_level=50,
    psf_sigma=None,
    flux_scale=None,
):
    """Generate a 16-bit synthetic star field image.

    Default behaviour models the Pi Cam Module 3 NoIR at 100 ms / gain 16:
    flux scales as 10^(-0.4·mag), per-pixel shot noise + read noise
    contribute realistic variance, and bright stars saturate. The PNG
    preview the simulator emits will then show a wide range of brightnesses.

    Pass `noise_level` to use the legacy Gaussian-noise model (used by
    the existing accuracy/coverage tests).
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

    # Two calibration paths.
    use_legacy = noise_level is not None
    if use_legacy:
        # Old behaviour — used by accuracy_test.py and coverage_test.py.
        sigma = float(psf_sigma if psf_sigma is not None else 2.0)
        bg = float(bg_level)
        fs = float(flux_scale if flux_scale is not None else 200000.0)
        flux = np.power(10, -0.4 * mags) * fs
        noise_sd = float(noise_level)
        use_shot = False
    else:
        # Pi Cam Module 3 NoIR model.
        cal = PI_CAM_BASE
        exp_s = float(exposure_ms) / 1000.0
        bg = cal["bg_at_gain1_per_sec"] * exp_s * gain
        # Total ADU per V=0 star, scaled by exposure and gain.
        fs = cal["flux_v0_at_gain1_per_s"] * exp_s * gain
        flux = np.power(10, -0.4 * mags) * fs
        sigma = float(psf_sigma if psf_sigma is not None else cal["psf_sigma_px"])
        # Read noise scales linearly with gain in this idealised model.
        read_noise_sd = cal["read_noise_at_gain1"] * gain
        # Electrons-per-ADU scales inversely with gain.
        electrons_per_adu = cal["electrons_per_adu_at_g1"] / gain
        noise_sd = None
        use_shot = True

    r = int(4 * sigma)
    norm = 1.0 / (2 * np.pi * sigma ** 2)
    two_sigma2 = 2.0 * sigma * sigma

    img = np.full((height, width), bg, dtype=np.float32)

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

    if use_legacy:
        if noise_sd > 0:
            img += np.random.normal(0, noise_sd, img.shape).astype(np.float32)
    else:
        # Shot noise: Poisson on electrons, scaled back to ADU. Approximate
        # with a Gaussian whose sd is sqrt(signal × electrons_per_ADU) / gain
        # — same variance as the Poisson at the electron count, expressed
        # back in ADU. Add read noise in quadrature.
        signal_clip = np.maximum(img, 0.0)
        var_adu = signal_clip / electrons_per_adu + read_noise_sd ** 2
        sd_adu = np.sqrt(var_adu, dtype=np.float32)
        img += np.random.normal(0.0, 1.0, img.shape).astype(np.float32) * sd_adu

    np.clip(img, 0, float(saturation_adu), out=img)
    return img.astype(np.uint16)


def save_pgm_16(img, path):
    """Atomic write of a 16-bit PGM file."""
    import os
    tmp = str(path) + ".tmp"
    with open(tmp, "wb") as f:
        f.write(f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode())
        f.write(img.astype(">u2").tobytes())
    os.replace(tmp, path)
