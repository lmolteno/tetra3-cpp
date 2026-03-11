#!/usr/bin/env python3
"""Star tracker development tool.

Generates synthetic star field images from Hipparcos, runs them through
the C++ solver pipeline, and shows the OLED display output in the terminal.

Usage:
    .venv/bin/python devtool.py [--solver-path ../../build/pi/pi_tracker]
"""

import argparse
import base64
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Static, Input, Label, Button
from textual.containers import Horizontal, Vertical
from textual.reactive import reactive
from textual import on

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

def generate_star_field(
    stars_dict,
    center_ra_deg,
    center_dec_deg,
    fov_deg=75.0,
    width=2028,
    height=1520,
    noise_level=8,
    bg_level=50,
):
    """Generate a 16-bit synthetic star field image.

    Uses gnomonic projection, Gaussian PSFs, and light read noise.
    Returns a numpy uint16 array.
    """
    center_ra = np.radians(center_ra_deg)
    center_dec = np.radians(center_dec_deg)
    fov_rad = np.radians(fov_deg)

    # Pixels per radian
    scale = width / (2.0 * np.tan(fov_rad / 2.0))

    cx = width / 2.0
    cy = height / 2.0

    # Start with background
    img = np.full((height, width), bg_level, dtype=np.float64)

    sin_d0 = np.sin(center_dec)
    cos_d0 = np.cos(center_dec)

    for hip, (ra, dec, mag) in stars_dict.items():
        # Quick angular distance check
        cos_c = sin_d0 * np.sin(dec) + cos_d0 * np.cos(dec) * np.cos(ra - center_ra)
        if cos_c <= 0.01:
            continue

        # Gnomonic projection
        dra = ra - center_ra
        sin_d = np.sin(dec)
        cos_d = np.cos(dec)
        cos_dra = np.cos(dra)
        sin_dra = np.sin(dra)

        x_proj = (cos_d * sin_dra) / cos_c
        y_proj = (cos_d0 * sin_d - sin_d0 * cos_d * cos_dra) / cos_c

        px = cx - x_proj * scale
        py = cy - y_proj * scale

        if px < -20 or px >= width + 20 or py < -20 or py >= height + 20:
            continue

        # Star brightness: flux ~ 10^(-0.4 * mag)
        # Scale so a mag 0 star peaks at ~40000 ADU, mag 5 at ~400 ADU
        flux = 10 ** (-0.4 * mag) * 200000

        # Gaussian PSF (sigma ~2.0 pixels for realistic seeing)
        sigma = 2.0
        r = int(4 * sigma)
        x0 = max(0, int(px) - r)
        x1 = min(width, int(px) + r + 1)
        y0 = max(0, int(py) - r)
        y1 = min(height, int(py) + r + 1)

        if x0 >= x1 or y0 >= y1:
            continue

        yy, xx = np.mgrid[y0:y1, x0:x1]
        gauss = np.exp(-((xx - px) ** 2 + (yy - py) ** 2) / (2 * sigma**2))
        gauss *= flux / (2 * np.pi * sigma**2)

        img[y0:y1, x0:x1] += gauss

    # Add light read noise only (no Poisson — keeps background clean)
    img += np.random.normal(0, noise_level, img.shape)

    # Clip and convert to uint16
    img = np.clip(img, 0, 65535).astype(np.uint16)

    return img


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

    def __init__(self, tracker_path, solver_path, db_stars, db_patterns, fov, **kwargs):
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

    def compose(self) -> ComposeResult:
        yield Header()
        yield Label("Synthetic star field → solve → OLED preview", id="top-row")
        with Horizontal(id="input-row"):
            yield Input(placeholder="RA (deg)", id="ra-input", value="84.05")
            yield Input(placeholder="Dec (deg)", id="dec-input", value="-1.2")
            yield Input(placeholder="FOV (deg)", id="fov-input", value="75")
            yield Button("Go", variant="primary", id="go-btn")
            yield Button("Show Image", id="show-img-btn")
            yield Button("No Solution", id="nosol-btn")
        with Horizontal(id="mode-row"):
            yield Button("Tracking", id="mode-tracking")
            yield Button("PA Sample", id="mode-pa-sampling")
            yield Button("PA Fix", id="mode-pa-fix")
        with Horizontal(id="preview-row"):
            yield OledPreview(id="oled")
            yield SkyPreview(id="sky")
        yield Static("Ready", id="status")
        yield Footer()

    def on_mount(self):
        self.set_status("Loading Hipparcos catalogue...")
        self.call_later(self._load_stars)

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
        except ValueError:
            self.set_status("Invalid coordinates")
            return

        sky = self.query_one("#sky", SkyPreview)
        sky.info = f"Generating...\nRA={ra_deg:.2f} Dec={dec_deg:.2f}\nFOV={fov_deg:.1f}deg"
        self.set_status("Generating star field...")
        self.call_later(lambda: self._run_solve(ra_deg, dec_deg, fov_deg))

    def _run_solve(self, ra_deg, dec_deg, fov_deg):
        sky = self.query_one("#sky", SkyPreview)

        # Generate synthetic image
        img = generate_star_field(
            self.stars, ra_deg, dec_deg, fov_deg=fov_deg,
            width=2028, height=1520,
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
            sky.info = (
                f"SOLVED\n"
                f"RA:  {result['ra_deg']:.4f} deg\n"
                f"Dec: {result['dec_deg']:.4f} deg\n"
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

    def on_unmount(self):
        if self.tracker_proc and self.tracker_proc.poll() is None:
            self.tracker_proc.terminate()
        if self._viewer_proc and self._viewer_proc.poll() is None:
            self._viewer_proc.terminate()


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
    args = parser.parse_args()

    app = DevToolApp(
        tracker_path=args.tracker,
        solver_path=args.solver,
        db_stars=args.db_stars,
        db_patterns=args.db_patterns,
        fov=args.fov,
    )
    app.run()


if __name__ == "__main__":
    main()
