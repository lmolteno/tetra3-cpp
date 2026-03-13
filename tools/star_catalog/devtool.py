#!/usr/bin/env python3
"""Star tracker development tool.

Generates synthetic star field images from Hipparcos, runs them through
the C++ solver pipeline, and shows the OLED display output in the terminal.

Usage:
    .venv/bin/python devtool.py [--solver-path ../../build/pi/pi_tracker]
"""

import argparse
import asyncio
import base64
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Static, Input, Label, Button
from textual.containers import Horizontal, Vertical
from textual.reactive import reactive
from textual import on

from celestial import imu_to_celestial, quaternion_to_alt_az_roll
from imu_server import IMUServer

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


# ---------------------------------------------------------------------------
# Synthetic star field image generation
# ---------------------------------------------------------------------------

def _prepare_star_arrays(stars_dict):
    """Convert star dict to contiguous arrays for vectorized operations.

    Call once, reuse across many generate_star_field calls.
    Returns (ra, dec, mag) arrays.
    """
    n = len(stars_dict)
    ra = np.empty(n, dtype=np.float64)
    dec = np.empty(n, dtype=np.float64)
    mag = np.empty(n, dtype=np.float64)
    for i, (r, d, m) in enumerate(stars_dict.values()):
        ra[i] = r
        dec[i] = d
        mag[i] = m
    return ra, dec, mag


# Module-level cache so arrays are built once
_star_arrays_cache = {}


def _get_star_arrays(stars_dict):
    """Get or create cached star arrays."""
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
):
    """Generate a 16-bit synthetic star field image.

    Uses gnomonic projection, Gaussian PSFs, and light read noise.
    roll_deg rotates the camera around the optical axis (degrees).
    Returns a numpy uint16 array.
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

    # Vectorized star projection
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

    flux = np.power(10, -0.4 * mags) * 200000

    # Stamp each star with a Gaussian PSF
    sigma = 2.0
    r = int(4 * sigma)  # 8 pixels
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

        # Use arange instead of mgrid — avoids creating 2D index arrays
        xs = np.arange(x0, x1, dtype=np.float32) - sx
        ys = np.arange(y0, y1, dtype=np.float32) - sy
        gx = np.exp(-xs * xs / two_sigma2)
        gy = np.exp(-ys * ys / two_sigma2)
        img[y0:y1, x0:x1] += (sf * norm) * np.outer(gy, gx)

    if noise_level > 0:
        img += np.random.normal(0, noise_level, img.shape).astype(np.float32)

    np.clip(img, 0, 65535, out=img)
    return img.astype(np.uint16)


def save_image_png(img, path):
    """Save a uint16 numpy array as 16-bit PNG."""
    try:
        from PIL import Image

        pil_img = Image.fromarray(img, mode="I;16")
        pil_img.save(path)
    except Exception:
        # Fallback: use raw binary that OpenCV can read
        import struct

        # Write as raw pgm (P5)
        with open(path, "wb") as f:
            header = f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode()
            f.write(header)
            f.write(img.astype(">u2").tobytes())


# ---------------------------------------------------------------------------
# OLED framebuffer rendering (braille characters)
# ---------------------------------------------------------------------------

def fb_to_halfblocks(fb_bytes):
    """Convert a 128x64 1-bit framebuffer to a half-block character string.

    The framebuffer uses column-page format: byte at offset (x + page*128)
    contains 8 vertical pixels for column x, rows page*8..page*8+7, LSB=top.

    Each character cell maps to 1x2 pixels using half-block characters:
      ▀ = top on,  bottom off
      ▄ = top off, bottom on
      █ = both on
      ' ' = both off
    Result: 128 columns x 32 rows.
    """
    if len(fb_bytes) < 1024:
        return f"Invalid framebuffer ({len(fb_bytes)} bytes)"
    fb_bytes = fb_bytes[:1024]

    def get_pixel(x, y):
        if x < 0 or x >= 128 or y < 0 or y >= 64:
            return False
        byte_idx = x + (y // 8) * 128
        return bool(fb_bytes[byte_idx] & (1 << (y % 8)))

    lines = []
    for row in range(0, 64, 2):
        line = ""
        for col in range(128):
            top = get_pixel(col, row)
            bot = get_pixel(col, row + 1)
            if top and bot:
                line += "\u2588"   # █
            elif top:
                line += "\u2580"   # ▀
            elif bot:
                line += "\u2584"   # ▄
            else:
                line += " "
        lines.append(line)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Textual App
# ---------------------------------------------------------------------------

class OledPreview(Static):
    """Renders the 128x64 OLED framebuffer using half-block characters."""

    DEFAULT_CSS = """
    OledPreview {
        width: 134;
        height: 36;
        background: black;
        color: white;
        border: round $primary;
        padding: 1 2;
    }
    """

    fb_data = reactive(b"\x00" * 1024)

    def watch_fb_data(self, value):
        self.update(fb_to_halfblocks(value))


class SkyPreview(Static):
    """Shows a small text summary of the generated sky field."""

    DEFAULT_CSS = """
    SkyPreview {
        width: 40;
        height: 18;
        border: round $accent;
        padding: 1;
    }
    """

    info = reactive("")

    def watch_info(self, value):
        self.update(value or "Enter RA, Dec and press Go")


_VIEWER_SCRIPT = """\
import sys, os, time, numpy as np, matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

path = sys.argv[1]
img = np.load(path)
fig, ax = plt.subplots(figsize=(10, 7))
vmax = max(np.percentile(img, 99.9), 1)
im = ax.imshow(img, cmap="gray", vmin=0, vmax=vmax)
ax.set_title("Synthetic star field")
fig.colorbar(im, ax=ax, label="ADU")
last_mtime = os.path.getmtime(path)

def poll_reload(event=None):
    global last_mtime
    try:
        mt = os.path.getmtime(path)
        if mt != last_mtime:
            last_mtime = mt
            new_img = np.load(path)
            vmax = max(np.percentile(new_img, 99.9), 1)
            im.set_data(new_img)
            im.set_clim(0, vmax)
            fig.canvas.draw_idle()
    except Exception:
        pass
    fig.canvas.manager.window.after(500, poll_reload)

fig.canvas.manager.window.after(500, poll_reload)
plt.show()
"""


class DevToolApp(App):
    """Star tracker development tool."""

    TITLE = "tetra3 devtool"
    CSS = """
    Screen {
        layout: vertical;
    }
    #top-row {
        height: auto;
        padding: 1;
    }
    #input-row {
        height: auto;
        padding: 0 1;
    }
    #input-row Input {
        width: 20;
        margin: 0 1;
    }
    #input-row Button {
        margin: 0 1;
    }
    #preview-row {
        height: auto;
        padding: 1;
    }
    #status {
        height: 3;
        padding: 0 1;
        color: $text-muted;
    }
    #mode-row {
        height: auto;
        padding: 0 1;
    }
    #mode-row Button {
        margin: 0 1;
    }
    """

    BINDINGS = [
        ("ctrl+q", "quit", "Quit"),
        ("ctrl+g", "go", "Solve"),
    ]

    def __init__(self, tracker_path, solver_path, db_stars, db_patterns, fov,
                 imu_port=8080, **kwargs):
        super().__init__(**kwargs)
        self.tracker_path = tracker_path
        self.solver_path = solver_path
        self.db_stars = db_stars
        self.db_patterns = db_patterns
        self.fov = fov
        self.stars = None
        self.tracker_proc = None
        self.last_image = None  # last generated star field (numpy uint16)
        self._viewer_proc = None
        # IMU state
        self._imu_mode = False
        self._imu_server = IMUServer(
            on_orientation_update=self._on_imu_orientation,
            on_command=self._on_phone_command,
            port=imu_port,
        )
        self._imu_ra_rad = 0.0
        self._imu_dec_rad = 0.0
        self._imu_roll_rad = 0.0
        self._imu_last_send = 0.0  # monotonic time of last tracker send
        self._imu_solver_task = None
        self._imu_last_fb_b64 = ""
        self._imu_last_solve = None  # last solver result dict

    def compose(self) -> ComposeResult:
        yield Header()
        yield Label("Synthetic star field → solve → OLED preview", id="top-row")
        with Horizontal(id="input-row"):
            yield Input(placeholder="RA (deg)", id="ra-input", value="84.05")
            yield Input(placeholder="Dec (deg)", id="dec-input", value="-1.2")
            yield Input(placeholder="FOV (deg)", id="fov-input", value="75")
            yield Input(placeholder="Roll (deg)", id="roll-input", value="0")
            yield Button("Go", variant="primary", id="go-btn")
            yield Button("Show Image", id="show-img-btn")
            yield Button("No Solution", id="nosol-btn")
        with Horizontal(id="mode-row"):
            yield Button("Tracking", id="mode-tracking")
            yield Button("PA Sample", id="mode-pa-sampling")
            yield Button("PA Fix", id="mode-pa-fix")
            yield Button("IMU", id="mode-imu")
            yield Input(placeholder="Lat", id="lat-input", value="")
            yield Input(placeholder="Lon", id="lon-input", value="")
        with Horizontal(id="preview-row"):
            yield OledPreview(id="oled")
            yield SkyPreview(id="sky")
        yield Static("Ready", id="status")
        yield Footer()

    def on_mount(self):
        self.set_status("Loading Hipparcos catalogue...")
        self.call_later(self._load_stars)
        # Start IMU server (always running, even if IMU mode not active yet)
        asyncio.create_task(self._start_imu_server())

    def _load_stars(self):
        self.stars = load_hipparcos_cached()
        self.set_status(f"Loaded {len(self.stars)} stars. Enter coordinates and press Go.")
        self._start_tracker()

    def _start_tracker(self):
        """Start pi_tracker in --stdin-json mode as a subprocess."""
        if self.tracker_proc and self.tracker_proc.poll() is None:
            self.tracker_proc.terminate()

        cmd = [self.tracker_path, "--stdin-json", "--display", "none"]
        try:
            self.tracker_proc = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )
            self.set_status("Tracker process started")
            # Send initial "no solution" so the OLED shows its boot state
            resp = self._send_to_tracker({"solved": False})
            self._update_oled(resp)
        except FileNotFoundError:
            self.set_status(f"ERROR: {self.tracker_path} not found. Build pi_tracker first.")
            self.tracker_proc = None

    def _send_to_tracker(self, msg):
        """Send a JSON line to the tracker and read the response."""
        if not self.tracker_proc or self.tracker_proc.poll() is not None:
            self._start_tracker()
            if not self.tracker_proc:
                return None

        try:
            self.tracker_proc.stdin.write(json.dumps(msg) + "\n")
            self.tracker_proc.stdin.flush()
            line = self.tracker_proc.stdout.readline()
            if not line:
                return None
            return json.loads(line.strip())
        except (BrokenPipeError, json.JSONDecodeError, OSError) as e:
            self.set_status(f"Tracker error: {e}")
            return None

    def _update_oled(self, response):
        """Update the OLED preview from a tracker response."""
        if not response:
            return
        fb_b64 = response.get("fb", "")
        if fb_b64:
            try:
                fb_bytes = base64.b64decode(fb_b64)
                self.query_one("#oled", OledPreview).fb_data = fb_bytes
            except Exception:
                pass

    def set_status(self, text):
        self.query_one("#status", Static).update(text)

    @on(Button.Pressed, "#go-btn")
    def on_go_btn(self):
        self._do_solve()

    def action_go(self):
        self._do_solve()

    @on(Button.Pressed, "#show-img-btn")
    def do_show_image(self):
        if self.last_image is None:
            self.set_status("No image yet — press Go first")
            return
        self.set_status("Opening image viewer...")
        self.call_later(self._show_image)

    def _show_image(self):
        # Save to temp file. A persistent viewer process watches for changes.
        tmp = os.path.join(tempfile.gettempdir(), "devtool_starfield.npy")
        np.save(tmp, self.last_image)

        if self._viewer_proc is None or self._viewer_proc.poll() is not None:
            self._viewer_proc = subprocess.Popen(
                [sys.executable, "-c", _VIEWER_SCRIPT, tmp],
                start_new_session=True,
            )
        self.set_status("Image viewer updated")

    @on(Button.Pressed, "#nosol-btn")
    def do_nosol(self):
        resp = self._send_to_tracker({"solved": False})
        self._update_oled(resp)
        self.set_status("Sent: no solution")
        self.query_one("#sky", SkyPreview).info = "No solution sent"

    @on(Button.Pressed, "#mode-tracking")
    def do_mode_tracking(self):
        resp = self._send_to_tracker({"mode": "tracking"})
        self._update_oled(resp)
        self.set_status("Mode: Tracking")

    @on(Button.Pressed, "#mode-pa-sampling")
    def do_mode_pa_sampling(self):
        resp = self._send_to_tracker({"mode": "pa_sampling"})
        self._update_oled(resp)
        self.set_status("Mode: PA Sampling")

    @on(Button.Pressed, "#mode-pa-fix")
    def do_mode_pa_fix(self):
        resp = self._send_to_tracker({"mode": "pa_fix"})
        self._update_oled(resp)
        self.set_status("Mode: PA Fix")

    def _do_solve(self):
        if not self.stars:
            self.set_status("Stars not loaded yet")
            return

        try:
            ra_deg = float(self.query_one("#ra-input", Input).value)
            dec_deg = float(self.query_one("#dec-input", Input).value)
            fov_deg = float(self.query_one("#fov-input", Input).value)
            roll_deg = float(self.query_one("#roll-input", Input).value)
        except ValueError:
            self.set_status("Invalid coordinates")
            return

        sky = self.query_one("#sky", SkyPreview)
        sky.info = f"Generating...\nRA={ra_deg:.2f} Dec={dec_deg:.2f}\nFOV={fov_deg:.1f}deg Roll={roll_deg:.1f}deg"
        self.set_status("Generating star field...")
        self.call_later(lambda: self._run_solve(ra_deg, dec_deg, fov_deg, roll_deg))

    def _run_solve(self, ra_deg, dec_deg, fov_deg, roll_deg=0.0):
        sky = self.query_one("#sky", SkyPreview)

        # Generate synthetic image
        img = generate_star_field(
            self.stars, ra_deg, dec_deg, fov_deg=fov_deg,
            roll_deg=roll_deg, width=2028, height=1520,
        )
        self.last_image = img

        # Update the live viewer if it's running
        if self._viewer_proc and self._viewer_proc.poll() is None:
            tmp_viewer = os.path.join(tempfile.gettempdir(), "devtool_starfield.npy")
            np.save(tmp_viewer, img)

        # Save to temp file
        with tempfile.NamedTemporaryFile(suffix=".pgm", delete=False) as f:
            tmp_path = f.name
            header = f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode()
            f.write(header)
            f.write(img.astype(">u2").tobytes())

        self.set_status("Running solver...")

        try:
            # Run solve_cli
            cmd = [
                self.solver_path,
                "--db-stars", self.db_stars,
                "--db-patterns", self.db_patterns,
                "--fov", str(self.fov),
                tmp_path,
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            # The solver library prints debug info to stdout before the JSON line.
            # Grab the last non-empty line which contains the JSON result.
            lines = [l for l in proc.stdout.strip().splitlines() if l.strip()]
            json_line = lines[-1] if lines else ""
            result = json.loads(json_line) if json_line else {"solved": False}
        except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError) as e:
            self.set_status(f"Solver error: {e}")
            sky.info = f"Solver error:\n{e}"
            os.unlink(tmp_path)
            return
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

        # Show result info
        if result.get("solved"):
            roll_deg_result = result.get('roll_rad', 0) * 180.0 / np.pi
            sky.info = (
                f"SOLVED\n"
                f"RA:  {result['ra_deg']:.4f} deg\n"
                f"Dec: {result['dec_deg']:.4f} deg\n"
                f"Roll: {roll_deg_result:.2f} deg\n"
                f"FOV: {result['fov_deg']:.2f} deg\n"
                f"RMSE: {result.get('rmse', 0):.2f}\"\n"
                f"Stars: {result.get('num_detected', '?')}\n"
                f"Time: {result.get('solve_time_ms', 0):.0f}ms"
            )
            self.set_status(
                f"Solved: RA={result['ra_deg']:.2f} Dec={result['dec_deg']:.2f}"
            )
        else:
            sky.info = (
                f"NO SOLUTION\n"
                f"Input: RA={ra_deg:.2f} Dec={dec_deg:.2f}\n"
                f"Stars detected: {result.get('num_detected', 0)}\n"
                f"Time: {result.get('solve_time_ms', 0):.0f}ms"
            )
            self.set_status("No solution found")

        # Send to tracker for OLED rendering
        resp = self._send_to_tracker(result)
        self._update_oled(resp)

    # ------------------------------------------------------------------
    # IMU mode
    # ------------------------------------------------------------------

    async def _start_imu_server(self):
        try:
            await self._imu_server.start()
            self.set_status(f"IMU server: {self._imu_server.url}")
        except Exception as e:
            self.set_status(f"IMU server failed: {e}")

    @on(Button.Pressed, "#mode-imu")
    def do_mode_imu(self):
        self._imu_mode = not self._imu_mode
        btn = self.query_one("#mode-imu", Button)
        if self._imu_mode:
            btn.variant = "success"
            self.set_status(f"IMU mode ON — open {self._imu_server.url} on phone")
            # Start periodic solver
            if self._imu_solver_task is None or self._imu_solver_task.done():
                self._imu_solver_task = asyncio.create_task(self._imu_solver_loop())
        else:
            btn.variant = "default"
            self.set_status("IMU mode OFF")
            if self._imu_solver_task and not self._imu_solver_task.done():
                self._imu_solver_task.cancel()

    async def _on_imu_orientation(self, q, lat, lon):
        """Called by IMUServer when phone sends orientation data."""
        if not self._imu_mode:
            return

        # Default GPS: use phone value, then devtool input, then give up
        if lat is None or lon is None:
            lat, lon = self._get_fallback_location()
            if lat is None:
                return

        try:
            alt_deg, az_deg, roll_deg_local = quaternion_to_alt_az_roll(q)
            ra_rad, dec_rad, roll_rad = imu_to_celestial(
                q, lat, lon, datetime.now(timezone.utc)
            )
        except Exception:
            return

        self._imu_ra_rad = ra_rad
        self._imu_dec_rad = dec_rad
        self._imu_roll_rad = roll_rad
        self._imu_alt_deg = alt_deg
        self._imu_az_deg = az_deg

        # Throttle tracker updates to ~10 Hz
        now = time.monotonic()
        if now - self._imu_last_send < 0.1:
            return
        self._imu_last_send = now

        fov_deg = self._get_fov_deg()
        ra_deg = math.degrees(ra_rad)
        dec_deg = math.degrees(dec_rad)
        roll_deg = math.degrees(roll_rad)

        # Send to pi_tracker (expects fov_rad)
        msg = {
            "solved": True,
            "ra_rad": ra_rad,
            "dec_rad": dec_rad,
            "roll_rad": roll_rad,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "fov_rad": math.radians(fov_deg),
        }
        resp = self._send_to_tracker(msg)
        self._update_oled(resp)

        # Extract framebuffer for phone broadcast
        fb_b64 = resp.get("fb", "") if resp else ""
        if fb_b64:
            self._imu_last_fb_b64 = fb_b64

        # Update RA/Dec inputs to reflect current pointing
        self._update_inputs_from_imu(ra_deg, dec_deg, roll_deg)

        # Broadcast to phone
        await self._imu_server.broadcast_state(
            fb_b64=self._imu_last_fb_b64,
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            roll_deg=roll_deg,
            fov_deg=fov_deg,
            alt_deg=alt_deg,
            az_deg=az_deg,
            solve_status=self._imu_last_solve,
        )

    def _update_inputs_from_imu(self, ra_deg, dec_deg, roll_deg):
        """Update the input fields with current IMU-derived coordinates."""
        try:
            self.query_one("#ra-input", Input).value = f"{ra_deg:.2f}"
            self.query_one("#dec-input", Input).value = f"{dec_deg:.2f}"
            self.query_one("#roll-input", Input).value = f"{roll_deg:.1f}"
        except Exception:
            pass

    async def _on_phone_command(self, cmd, value):
        """Called by IMUServer when phone sends a command."""
        if cmd == "solve":
            self._trigger_imu_solve()
        elif cmd == "mode":
            mode_map = {
                "tracking": "mode-tracking",
                "pa_sampling": "mode-pa-sampling",
                "pa_fix": "mode-pa-fix",
            }
            btn_id = mode_map.get(value)
            if btn_id:
                self.query_one(f"#{btn_id}", Button).press()
        elif cmd == "fov":
            try:
                fov = float(value)
                self.query_one("#fov-input", Input).value = str(fov)
            except (ValueError, TypeError):
                pass

    def _trigger_imu_solve(self):
        """Trigger an immediate solver run with current IMU coordinates."""
        if not self.stars:
            return
        ra_deg = math.degrees(self._imu_ra_rad)
        dec_deg = math.degrees(self._imu_dec_rad)
        roll_deg = math.degrees(self._imu_roll_rad)
        fov_deg = self._get_fov_deg()
        self.call_later(lambda: self._run_imu_solve(ra_deg, dec_deg, fov_deg, roll_deg))

    async def _imu_solver_loop(self):
        """Periodic solver loop: run every ~3 seconds while IMU mode is active."""
        while self._imu_mode:
            await asyncio.sleep(3.0)
            if not self._imu_mode or not self.stars:
                continue
            ra_deg = math.degrees(self._imu_ra_rad)
            dec_deg = math.degrees(self._imu_dec_rad)
            roll_deg = math.degrees(self._imu_roll_rad)
            fov_deg = self._get_fov_deg()
            # Run solver in thread to avoid blocking the event loop
            loop = asyncio.get_event_loop()
            await loop.run_in_executor(
                None, self._run_imu_solve, ra_deg, dec_deg, fov_deg, roll_deg
            )

    def _run_imu_solve(self, ra_deg, dec_deg, fov_deg, roll_deg):
        """Run the solver for IMU validation and update sky preview."""
        sky = self.query_one("#sky", SkyPreview)
        sky.info = (
            f"IMU Solve...\n"
            f"RA={ra_deg:.2f} Dec={dec_deg:.2f}\n"
            f"Roll={roll_deg:.1f} FOV={fov_deg:.1f}"
        )

        img = generate_star_field(
            self.stars, ra_deg, dec_deg, fov_deg=fov_deg,
            roll_deg=roll_deg, width=2028, height=1520,
        )
        self.last_image = img

        with tempfile.NamedTemporaryFile(suffix=".pgm", delete=False) as f:
            tmp_path = f.name
            header = f"P5\n{img.shape[1]} {img.shape[0]}\n65535\n".encode()
            f.write(header)
            f.write(img.astype(">u2").tobytes())

        try:
            cmd = [
                self.solver_path,
                "--db-stars", self.db_stars,
                "--db-patterns", self.db_patterns,
                "--fov", str(self.fov),
                tmp_path,
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            lines = [l for l in proc.stdout.strip().splitlines() if l.strip()]
            json_line = lines[-1] if lines else ""
            result = json.loads(json_line) if json_line else {"solved": False}
        except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError) as e:
            sky.info = f"IMU solve error:\n{e}"
            return
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

        self._imu_last_solve = result

        if result.get("solved"):
            solver_ra = result["ra_deg"]
            solver_dec = result["dec_deg"]
            delta_ra = solver_ra - ra_deg
            delta_dec = solver_dec - dec_deg
            sky.info = (
                f"IMU vs Solver\n"
                f"IMU:    RA={ra_deg:.2f} Dec={dec_deg:.2f}\n"
                f"Solver: RA={solver_ra:.2f} Dec={solver_dec:.2f}\n"
                f"Delta:  dRA={delta_ra:+.2f} dDec={delta_dec:+.2f}\n"
                f"RMSE: {result.get('rmse', 0):.2f}\"\n"
                f"Time: {result.get('solve_time_ms', 0):.0f}ms"
            )
            # Send corrected result to tracker
            resp = self._send_to_tracker(result)
            self._update_oled(resp)
        else:
            sky.info = (
                f"IMU: RA={ra_deg:.2f} Dec={dec_deg:.2f}\n"
                f"Solver: NO SOLUTION\n"
                f"Stars: {result.get('num_detected', 0)}\n"
                f"Time: {result.get('solve_time_ms', 0):.0f}ms"
            )

    def _get_fov_deg(self):
        """Get current FOV from input field."""
        try:
            return float(self.query_one("#fov-input", Input).value)
        except (ValueError, Exception):
            return 75.0

    def _get_fallback_location(self):
        """Get lat/lon from the devtool input fields, or None."""
        try:
            lat = float(self.query_one("#lat-input", Input).value)
            lon = float(self.query_one("#lon-input", Input).value)
            return lat, lon
        except (ValueError, Exception):
            return None, None

    def on_unmount(self):
        if self.tracker_proc and self.tracker_proc.poll() is None:
            self.tracker_proc.terminate()
        if self._viewer_proc and self._viewer_proc.poll() is None:
            self._viewer_proc.terminate()
        if self._imu_solver_task and not self._imu_solver_task.done():
            self._imu_solver_task.cancel()
        asyncio.create_task(self._imu_server.stop())


def main():
    script_dir = Path(__file__).parent
    project_root = script_dir / ".." / ".."

    parser = argparse.ArgumentParser(description="Star tracker dev tool")
    parser.add_argument(
        "--tracker",
        default=str(project_root / "build" / "pi" / "pi_tracker"),
        help="Path to pi_tracker binary",
    )
    parser.add_argument(
        "--solver",
        default=str(project_root / "build" / "pi" / "solve_cli"),
        help="Path to solve_cli binary",
    )
    parser.add_argument(
        "--db-stars",
        default=str(project_root / "tetra3_db_stars.bin"),
        help="Star database path",
    )
    parser.add_argument(
        "--db-patterns",
        default=str(project_root / "tetra3_db_patterns.bin"),
        help="Pattern database path",
    )
    parser.add_argument("--fov", type=float, default=75.0, help="FOV estimate for solver (deg)")
    parser.add_argument("--imu-port", type=int, default=8080, help="IMU server port")
    args = parser.parse_args()

    app = DevToolApp(
        tracker_path=args.tracker,
        solver_path=args.solver,
        db_stars=args.db_stars,
        db_patterns=args.db_patterns,
        fov=args.fov,
        imu_port=args.imu_port,
    )
    app.run()


if __name__ == "__main__":
    main()
