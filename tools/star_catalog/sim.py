#!/usr/bin/env python3
"""Interactive star-tracker simulator — web edition.

Fake hardware for testing pi_tracker without a Pi or camera. Open a browser,
point the synthetic camera around the sky, watch the OLED framebuffer
update in real time, press buttons that mirror the two GPIO buttons on the
real device.

Useful for:
  - Walking through every UI state (tracking, PA sampling, PA fix).
  - End-to-end PA workflow testing — flip on "mount" mode and the camera
    follows a virtual misaligned RA axis. Turn the alt/az knobs and watch
    the live PAFixTracker chase the offset down.
  - Sanity-checking solve accuracy under noise.

Run:
    cd <repo>
    tools/star_catalog/.venv/bin/python tools/star_catalog/sim.py

Then open http://localhost:8765 in a browser.
"""

import argparse
import asyncio
import base64
import io
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from aiohttp import web, WSMsgType

sys.path.insert(0, str(Path(__file__).parent))
from imaging import (generate_star_field, load_binary_catalog,
                     load_hipparcos_cached, save_pgm_16)


# ---------------------------------------------------------------------------
# Sky geometry helpers

def unit(ra_rad, dec_rad):
    cd = np.cos(dec_rad)
    return np.array([cd * np.cos(ra_rad), cd * np.sin(ra_rad), np.sin(dec_rad)])


def radec(v):
    dec = float(np.arcsin(np.clip(v[2], -1.0, 1.0)))
    ra = float(np.arctan2(v[1], v[0]) % (2 * np.pi))
    return ra, dec


def rotate_axis(v, axis, theta):
    c, s = np.cos(theta), np.sin(theta)
    return v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1.0 - c)


def perp_unit(v):
    if abs(v[2]) < 0.99:
        e = np.cross(v, np.array([0.0, 0.0, 1.0]))
    else:
        e = np.cross(v, np.array([1.0, 0.0, 0.0]))
    return e / np.linalg.norm(e)


# ---------------------------------------------------------------------------
# Image -> PNG (downsampled, contrast-stretched preview for the browser).

def img_to_png(img_u16, max_dim=512):
    from PIL import Image

    h, w = img_u16.shape
    scale = max(w / max_dim, h / max_dim, 1.0)
    if scale > 1.0:
        new_w = int(w / scale)
        new_h = int(h / scale)
        # Block-average downsample for speed; PIL.LANCZOS would be slightly
        # nicer but a few percent of overall round-trip latency is unnoticeable.
        sy = h // new_h
        sx = w // new_w
        img_u16 = img_u16[: sy * new_h, : sx * new_w].reshape(
            new_h, sy, new_w, sx).mean(axis=(1, 3))

    arr = np.asarray(img_u16, dtype=np.float32)
    # Contrast-stretch by percentile so stars pop without being saturated.
    lo = float(np.percentile(arr, 10))
    hi = float(np.percentile(arr, 99.9))
    if hi <= lo:
        hi = lo + 1.0
    out = np.clip((arr - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)
    im = Image.fromarray(out, mode="L")
    buf = io.BytesIO()
    im.save(buf, format="PNG")
    return buf.getvalue()


# ---------------------------------------------------------------------------
# pi_tracker subprocess

class Tracker:
    def __init__(self, *, tracker, solver, db_stars, db_patterns,
                 image_path, frames_dir, state_dir, width, height, fov):
        env = os.environ.copy()
        env["LIBCAMERA_LOG_LEVELS"] = "*:ERROR"
        self.proc = subprocess.Popen(
            [str(tracker),
             "--image", str(image_path), "--image-watch",
             "--solver", str(solver),
             "--db-stars", str(db_stars), "--db-patterns", str(db_patterns),
             "--width", str(width), "--height", str(height),
             "--fov", str(fov),
             "--max-fps", "10",
             "--frames-dir", str(frames_dir),
             "--state-dir", str(state_dir),
             "--no-distortion-k"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True, env=env,
        )

    def send(self, cmd: dict):
        try:
            self.proc.stdin.write(json.dumps(cmd) + "\n")
            self.proc.stdin.flush()
        except (BrokenPipeError, OSError):
            pass

    def close(self):
        for fn in (self.proc.stdin.close, self.proc.terminate):
            try: fn()
            except Exception: pass
        try: self.proc.wait(timeout=2)
        except Exception:
            try: self.proc.kill()
            except Exception: pass


# ---------------------------------------------------------------------------
# Simulator state

class Sim:
    def __init__(self, args, stars):
        self.args = args
        self.stars = stars

        # Camera state.
        self.cam_ra_deg   = 0.0
        self.cam_dec_deg  = 0.0
        self.cam_roll_deg = 0.0
        self.cam_fov_deg  = args.fov
        self.noise        = args.noise

        # Mount state. When mount_enabled, the camera direction is derived
        # from the (misaligned) mount + rotation around its RA axis.
        self.mount_enabled               = False
        self.mount_alt_offset_arcmin     = 30.0
        self.mount_az_offset_arcmin      = 0.0
        self.mount_ra_axis_deg           = 0.0
        self.target_offset_deg           = 8.0

        # Temp dir for the watched frame + tracker working area.
        self._tmp = Path(tempfile.mkdtemp(prefix="tetra3_sim_"))
        self.image_path = self._tmp / "sim.pgm"
        (self._tmp / "frames").mkdir()
        (self._tmp / "state").mkdir()

        # Image cache (PNG bytes) and version counter for cache busting.
        self.image_png = b""
        self.image_version = 0

        # Latest daemon state tick.
        self.last_state = {}

        # WS clients to push state updates to.
        self.subscribers = set()

        # Lock to serialise image regeneration.
        self._regen_lock = asyncio.Lock()

        self.tracker = Tracker(
            tracker=args.tracker, solver=args.solver,
            db_stars=args.db_stars, db_patterns=args.db_patterns,
            image_path=self.image_path,
            frames_dir=self._tmp / "frames",
            state_dir=self._tmp / "state",
            width=args.width, height=args.height, fov=args.fov,
        )

    # ---- camera direction ----

    def cam_radec_now(self):
        """Return the synthesizer's camera (RA, Dec) in degrees."""
        if not self.mount_enabled:
            return self.cam_ra_deg, self.cam_dec_deg
        # Derive from mount.
        true_pole_dec = np.pi / 2 if self.args.latitude >= 0 else -np.pi / 2
        offset_total = np.hypot(self.mount_alt_offset_arcmin,
                                self.mount_az_offset_arcmin)
        offset_rad = np.radians(offset_total / 60.0)
        if offset_total < 1e-6:
            mount_pole_v = unit(0.0, true_pole_dec)
        else:
            az_dir = np.arctan2(self.mount_az_offset_arcmin,
                                self.mount_alt_offset_arcmin)
            pole = unit(0.0, true_pole_dec)
            e1 = perp_unit(pole)
            e2 = np.cross(pole, e1)
            perp = np.cos(az_dir) * e1 + np.sin(az_dir) * e2
            mount_pole_v = pole * np.cos(offset_rad) + perp * np.sin(offset_rad)
        init = (mount_pole_v * np.cos(np.radians(self.target_offset_deg))
              + perp_unit(mount_pole_v) * np.sin(np.radians(self.target_offset_deg)))
        cam_v = rotate_axis(init, mount_pole_v,
                            np.radians(self.mount_ra_axis_deg))
        ra_rad, dec_rad = radec(cam_v)
        return float(np.degrees(ra_rad)), float(np.degrees(dec_rad))

    # ---- image regeneration ----

    async def regenerate_image(self):
        async with self._regen_lock:
            ra_deg, dec_deg = self.cam_radec_now()
            # Heavy lifting off-thread.
            loop = asyncio.get_event_loop()
            img = await loop.run_in_executor(
                None,
                lambda: generate_star_field(
                    self.stars, ra_deg, dec_deg,
                    fov_deg=self.cam_fov_deg, roll_deg=self.cam_roll_deg,
                    width=self.args.width, height=self.args.height,
                    noise_level=self.noise, psf_sigma=1.0, flux_scale=500000,
                ))
            await loop.run_in_executor(
                None, lambda: save_pgm_16(img, self.image_path))
            self.image_png = await loop.run_in_executor(None, img_to_png, img)
            self.image_version += 1

    # ---- WS broadcast ----

    async def broadcast(self, msg: dict):
        dead = []
        for ws in list(self.subscribers):
            try:
                await ws.send_json(msg)
            except Exception:
                dead.append(ws)
        for ws in dead:
            self.subscribers.discard(ws)

    def public_state(self):
        s = self.last_state
        ra_deg, dec_deg = self.cam_radec_now()
        pa = s.get("pa") or {}
        return {
            "cam": {
                "ra_deg":  ra_deg, "dec_deg": dec_deg,
                "fov_deg": self.cam_fov_deg, "roll_deg": self.cam_roll_deg,
                "noise":   self.noise,
            },
            "mount": {
                "enabled":          self.mount_enabled,
                "alt_offset_arcmin": self.mount_alt_offset_arcmin,
                "az_offset_arcmin":  self.mount_az_offset_arcmin,
                "ra_axis_deg":       self.mount_ra_axis_deg,
                "target_offset_deg": self.target_offset_deg,
            },
            "daemon": {
                "mode":    s.get("mode", "?"),
                "solved":  bool(s.get("solved")),
                "ra_deg":  float(np.degrees(s["ra_rad"])) if s.get("solved") else None,
                "dec_deg": float(np.degrees(s["dec_rad"])) if s.get("solved") else None,
                "roll_deg":float(np.degrees(s["roll_rad"])) if s.get("solved") else None,
                "rmse":    s.get("rmse"),
                "matches": s.get("num_matches"),
            },
            "pa": {
                "num_samples":   pa.get("num_samples", 0),
                "pole_valid":    pa.get("pole_valid", False),
                "arc_deg":       pa.get("arc_deg"),
                "alt_arcmin":    pa.get("alt_arcmin"),
                "az_arcmin":     pa.get("az_arcmin"),
                "total_arcmin":  pa.get("total_arcmin"),
            },
            "fb_b64": s.get("fb", ""),
            "image_version": self.image_version,
        }

    # ---- tracker reader ----

    async def tracker_reader(self):
        loop = asyncio.get_event_loop()
        proc = self.tracker.proc
        while True:
            line = await loop.run_in_executor(None, proc.stdout.readline)
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            try:
                self.last_state = json.loads(line)
            except json.JSONDecodeError:
                continue
            await self.broadcast({"type": "tick", "state": self.public_state()})

    # ---- commands from the browser ----

    async def handle_cmd(self, cmd: dict):
        kind = cmd.get("cmd")
        if kind == "pan":
            if self.mount_enabled:
                return  # ignored
            self.cam_ra_deg  = float(cmd.get("ra",  self.cam_ra_deg))  % 360
            self.cam_dec_deg = float(np.clip(cmd.get("dec", self.cam_dec_deg),
                                             -89.5, 89.5))
            await self.regenerate_image()
        elif kind == "pan_delta":
            if self.mount_enabled:
                return
            self.cam_ra_deg  = (self.cam_ra_deg  + float(cmd.get("dra",  0))) % 360
            self.cam_dec_deg = float(np.clip(
                self.cam_dec_deg + float(cmd.get("ddec", 0)), -89.5, 89.5))
            await self.regenerate_image()
        elif kind == "roll":
            self.cam_roll_deg = (float(cmd.get("deg", self.cam_roll_deg))
                                 + 360.0) % 360.0
            await self.regenerate_image()
        elif kind == "fov":
            self.cam_fov_deg = float(np.clip(cmd.get("deg", self.cam_fov_deg),
                                              5.0, 120.0))
            await self.regenerate_image()
        elif kind == "noise":
            self.noise = int(cmd.get("value", self.noise))
            await self.regenerate_image()
        elif kind == "mount_toggle":
            self.mount_enabled = bool(cmd.get("on", not self.mount_enabled))
            await self.regenerate_image()
        elif kind == "mount_adjust":
            self.mount_alt_offset_arcmin += float(cmd.get("dalt", 0))
            self.mount_az_offset_arcmin  += float(cmd.get("daz",  0))
            await self.regenerate_image()
        elif kind == "mount_rotate":
            self.mount_ra_axis_deg = (self.mount_ra_axis_deg
                                     + float(cmd.get("ddeg", 0))) % 360.0
            await self.regenerate_image()
        elif kind == "mount_set":
            if "alt_offset_arcmin" in cmd:
                self.mount_alt_offset_arcmin = float(cmd["alt_offset_arcmin"])
            if "az_offset_arcmin" in cmd:
                self.mount_az_offset_arcmin = float(cmd["az_offset_arcmin"])
            if "ra_axis_deg" in cmd:
                self.mount_ra_axis_deg = float(cmd["ra_axis_deg"]) % 360.0
            if "target_offset_deg" in cmd:
                self.target_offset_deg = float(cmd["target_offset_deg"])
            await self.regenerate_image()
        elif kind == "button":
            which = cmd.get("which")
            s = self.last_state
            mode = s.get("mode", "tracking")
            pa_samples = (s.get("pa") or {}).get("num_samples", 0)
            if which == "forward":
                if mode == "tracking":
                    self.tracker.send({"mode": "pa_sampling"})
                elif mode == "pa_sampling":
                    if pa_samples >= 3:
                        self.tracker.send({"mode": "pa_fix"})
                elif mode == "pa_fix":
                    self.tracker.send({"mode": "tracking"})
            elif which == "back":
                if mode != "tracking":
                    self.tracker.send({"mode": "tracking"})
        elif kind == "mode":
            # Direct mode setting (debug / scripted use).
            self.tracker.send({"mode": cmd.get("mode", "tracking")})

        # Push fresh state to all subscribers after any command.
        await self.broadcast({"type": "tick", "state": self.public_state()})


# ---------------------------------------------------------------------------
# HTML / JS

INDEX_HTML = r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>tetra3 sim</title>
<style>
  body { font: 13px/1.4 sans-serif; background:#181820; color:#eee;
         margin:0; padding:0 12px 12px; }
  h1 { font-size: 14px; margin: 8px 0; color:#8ad; font-weight:600; }
  .row { display:flex; gap:12px; align-items:flex-start; flex-wrap:wrap; }
  .panel { background:#22222a; border-radius:6px; padding:8px 12px; }
  #sky { image-rendering: pixelated; max-width: 512px; border:1px solid #444;
         background:#000; }
  #oled-wrap { width: 384px; }
  #oled { width: 384px; height: 192px; image-rendering: pixelated;
          background:#000; border:1px solid #444; }
  .state-grid { display:grid; grid-template-columns:auto auto auto auto;
                gap:2px 14px; font-family:monospace; }
  .state-grid b { color:#aaf; font-weight:400; }
  button { background:#345; color:#eee; border:1px solid #678;
           border-radius:4px; padding:4px 10px; font:13px/1.2 sans-serif;
           cursor:pointer; }
  button:hover { background:#456; }
  button.primary { background:#3a6; border-color:#5b8; color:#fff; }
  button.primary:hover { background:#4b7; }
  button.back { background:#a44; border-color:#c66; color:#fff; }
  button.back:hover { background:#b55; }
  .btnrow { display:flex; gap:6px; flex-wrap:wrap; }
  label { color:#9ad; }
  input[type=number] { width:5em; background:#181820; color:#eee;
                       border:1px solid #555; padding:2px 4px; border-radius:3px; }
  input[type=range] { vertical-align:middle; }
  .keys { color:#888; font-size:11px; font-family:monospace; }
  #status { color:#9c9; min-height:18px; font-family:monospace; }
  .pa-good { color: #6f6; }
  .pa-bad  { color: #f88; }
</style>
</head>
<body>

<h1>tetra3 simulator <span id=status></span></h1>

<div class=row>
  <div class=panel>
    <div>Synthetic camera <span class=keys>(2028×1520, scaled)</span></div>
    <img id=sky width=512>
  </div>
  <div class=panel id=oled-wrap>
    <div>OLED framebuffer <span class=keys>(128×64, 3×)</span></div>
    <canvas id=oled width=384 height=192></canvas>
    <div class=btnrow style="margin-top:6px">
      <button class=primary onclick="cmd({cmd:'button', which:'forward'})">
        ⬛ Forward (GPIO 1)
      </button>
      <button class=back onclick="cmd({cmd:'button', which:'back'})">
        ⬛ Back (GPIO 2)
      </button>
    </div>
    <div class=keys style="margin-top:4px">
      shortcut: SPACE / BACKSPACE
    </div>
  </div>

  <div class=panel id=statepanel>
    <div>State</div>
    <div class=state-grid id=stategrid></div>
  </div>
</div>

<div class=row style="margin-top:10px">
  <div class=panel>
    <div>Free pointing <span class=keys>(disabled in mount mode)</span></div>
    <table>
      <tr><td>RA</td>
          <td><input id=ra type=number step=0.1></td>
          <td><button onclick="bump('ra', -5)">−5°</button>
              <button onclick="bump('ra', -1)">−1°</button>
              <button onclick="bump('ra', +1)">+1°</button>
              <button onclick="bump('ra', +5)">+5°</button></td></tr>
      <tr><td>Dec</td>
          <td><input id=dec type=number step=0.1></td>
          <td><button onclick="bump('dec', -5)">−5°</button>
              <button onclick="bump('dec', -1)">−1°</button>
              <button onclick="bump('dec', +1)">+1°</button>
              <button onclick="bump('dec', +5)">+5°</button></td></tr>
      <tr><td>Roll</td>
          <td><input id=roll type=number step=0.1></td>
          <td><button onclick="setRoll(roll.value - 5)">−5°</button>
              <button onclick="setRoll(roll.value - 1)">−1°</button>
              <button onclick="setRoll(roll.value + 1)">+1°</button>
              <button onclick="setRoll(roll.value + 5)">+5°</button></td></tr>
      <tr><td>FOV</td>
          <td><input id=fov type=number step=1></td>
          <td><button onclick="setFov(fov.value - 5)">−5°</button>
              <button onclick="setFov(fov.value + 5)">+5°</button></td></tr>
      <tr><td>Noise</td>
          <td colspan=2>
            <select id=noise onchange="cmd({cmd:'noise', value:+this.value})">
              <option>0</option><option>4</option>
              <option selected>8</option><option>16</option><option>32</option>
            </select>
            <span class=keys>(σ; 0 fails — detector needs some BG noise)</span>
          </td></tr>
    </table>
    <div class=keys>Arrow keys pan, Shift = 5°, Ctrl = 0.1°</div>
  </div>

  <div class=panel>
    <div>Mount <input type=checkbox id=mounten
      onchange="cmd({cmd:'mount_toggle', on:this.checked})"></div>
    <table>
      <tr><td>Alt offset</td>
          <td><span id=alt_off>30.0</span>′</td>
          <td><button onclick="cmd({cmd:'mount_adjust', dalt:-5})">−5′</button>
              <button onclick="cmd({cmd:'mount_adjust', dalt:-1})">−1′</button>
              <button onclick="cmd({cmd:'mount_adjust', dalt:-0.2})">−0.2′</button>
              <button onclick="cmd({cmd:'mount_adjust', dalt:+0.2})">+0.2′</button>
              <button onclick="cmd({cmd:'mount_adjust', dalt:+1})">+1′</button>
              <button onclick="cmd({cmd:'mount_adjust', dalt:+5})">+5′</button></td></tr>
      <tr><td>Az offset</td>
          <td><span id=az_off>0.0</span>′</td>
          <td><button onclick="cmd({cmd:'mount_adjust', daz:-5})">−5′</button>
              <button onclick="cmd({cmd:'mount_adjust', daz:-1})">−1′</button>
              <button onclick="cmd({cmd:'mount_adjust', daz:-0.2})">−0.2′</button>
              <button onclick="cmd({cmd:'mount_adjust', daz:+0.2})">+0.2′</button>
              <button onclick="cmd({cmd:'mount_adjust', daz:+1})">+1′</button>
              <button onclick="cmd({cmd:'mount_adjust', daz:+5})">+5′</button></td></tr>
      <tr><td>RA axis</td>
          <td><span id=ra_axis>0.0</span>°</td>
          <td><button onclick="cmd({cmd:'mount_rotate', ddeg:-15})">⟲ 15°</button>
              <button onclick="cmd({cmd:'mount_rotate', ddeg:-5})">⟲ 5°</button>
              <button onclick="cmd({cmd:'mount_rotate', ddeg:+5})">5° ⟳</button>
              <button onclick="cmd({cmd:'mount_rotate', ddeg:+15})">15° ⟳</button></td></tr>
      <tr><td colspan=3 class=keys>
        Mount mode locks RA/Dec; the camera follows a virtual mount pivoting
        around its (mis)aligned RA axis. Drive the offset down with the
        alt/az "knobs" and watch the PA fix chart converge.
      </td></tr>
    </table>
  </div>
</div>

<script>
const sky = document.getElementById('sky');
const oled = document.getElementById('oled');
const ctx = oled.getContext('2d');
const status_el = document.getElementById('status');
const grid = document.getElementById('stategrid');
const ra = document.getElementById('ra');
const dec = document.getElementById('dec');
const roll = document.getElementById('roll');
const fov = document.getElementById('fov');
const noise = document.getElementById('noise');
const mounten = document.getElementById('mounten');
const alt_off = document.getElementById('alt_off');
const az_off = document.getElementById('az_off');
const ra_axis = document.getElementById('ra_axis');

const proto = location.protocol === 'https:' ? 'wss:' : 'ws:';
const ws = new WebSocket(proto + '//' + location.host + '/ws');

let lastVersion = -1;

ws.onopen  = () => { status_el.textContent = 'connected'; };
ws.onclose = () => { status_el.textContent = 'disconnected — reload to retry'; };

ws.onmessage = (ev) => {
  const m = JSON.parse(ev.data);
  if (m.type === 'tick') {
    render(m.state);
  }
};

function cmd(o) { ws.send(JSON.stringify(o)); }
function bump(field, delta) {
  const cur = parseFloat(document.getElementById(field).value) || 0;
  document.getElementById(field).value = (cur + delta).toFixed(1);
  cmd({cmd:'pan_delta', dra: field==='ra' ? delta : 0,
                        ddec: field==='dec' ? delta : 0});
}
function setRoll(v) { cmd({cmd:'roll', deg: parseFloat(v)}); }
function setFov(v)  { cmd({cmd:'fov',  deg: parseFloat(v)}); }

// Allow direct input editing
ra.addEventListener('change',
  () => cmd({cmd:'pan', ra:+ra.value, dec:+dec.value}));
dec.addEventListener('change',
  () => cmd({cmd:'pan', ra:+ra.value, dec:+dec.value}));
roll.addEventListener('change', () => setRoll(roll.value));
fov.addEventListener('change',  () => setFov(fov.value));

// Keyboard shortcuts: arrows pan, space/backspace = buttons
document.addEventListener('keydown', (e) => {
  if (e.target.tagName === 'INPUT') return;
  let dra = 0, ddec = 0;
  const step = e.shiftKey ? 5 : (e.ctrlKey ? 0.1 : 1);
  if      (e.key === 'ArrowLeft')  dra = -step;
  else if (e.key === 'ArrowRight') dra = +step;
  else if (e.key === 'ArrowDown')  ddec = -step;
  else if (e.key === 'ArrowUp')    ddec = +step;
  else if (e.key === ' ')   { cmd({cmd:'button', which:'forward'}); e.preventDefault(); return; }
  else if (e.key === 'Backspace') { cmd({cmd:'button', which:'back'}); e.preventDefault(); return; }
  else if (e.key === '[')   { setRoll(parseFloat(roll.value)-1); return; }
  else if (e.key === ']')   { setRoll(parseFloat(roll.value)+1); return; }
  else if (e.key === '+' || e.key === '=') { setFov(parseFloat(fov.value)-2); return; }
  else if (e.key === '-')   { setFov(parseFloat(fov.value)+2); return; }
  else return;
  cmd({cmd:'pan_delta', dra, ddec});
  e.preventDefault();
});

// OLED framebuffer rendering (SH1106 page-major).
function drawFb(b64) {
  const raw = atob(b64);
  if (raw.length < 1024) return;
  const id = ctx.createImageData(128, 64);
  for (let y = 0; y < 64; y++) {
    const page = y >> 3;
    const bit = 1 << (y & 7);
    for (let x = 0; x < 128; x++) {
      const on = raw.charCodeAt(x + page * 128) & bit;
      const off = (y * 128 + x) * 4;
      const v = on ? 255 : 0;
      id.data[off] = v;
      id.data[off+1] = v;
      id.data[off+2] = v;
      id.data[off+3] = 255;
    }
  }
  // Draw 1:1 onto a 128x64 backing, then scale up via the canvas's size.
  const tmp = document.createElement('canvas');
  tmp.width = 128; tmp.height = 64;
  tmp.getContext('2d').putImageData(id, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(tmp, 0, 0, 384, 192);
}

function fmt(n, d=2) { return (n == null) ? '—' : n.toFixed(d); }

function render(s) {
  // Update sky image only on version change.
  if (s.image_version !== lastVersion) {
    lastVersion = s.image_version;
    sky.src = '/sky.png?v=' + lastVersion;
  }
  if (s.fb_b64) drawFb(s.fb_b64);

  // Sync input fields (without firing change events).
  if (document.activeElement !== ra)   ra.value   = fmt(s.cam.ra_deg, 3);
  if (document.activeElement !== dec)  dec.value  = fmt(s.cam.dec_deg, 3);
  if (document.activeElement !== roll) roll.value = fmt(s.cam.roll_deg, 1);
  if (document.activeElement !== fov)  fov.value  = fmt(s.cam.fov_deg, 1);
  noise.value = String(s.cam.noise);
  mounten.checked = s.mount.enabled;
  alt_off.textContent = fmt(s.mount.alt_offset_arcmin, 1);
  az_off.textContent  = fmt(s.mount.az_offset_arcmin, 1);
  ra_axis.textContent = fmt(s.mount.ra_axis_deg, 1);

  const rows = [];
  rows.push(['Mode',   s.daemon.mode]);
  if (s.daemon.solved) {
    rows.push(['Solve',
      `RA ${fmt(s.daemon.ra_deg,3)}°  Dec ${fmt(s.daemon.dec_deg,3)}°`]);
    rows.push(['',
      `roll ${fmt(s.daemon.roll_deg,1)}°  RMSE ${fmt(s.daemon.rmse,1)}″`
      + `  matches ${s.daemon.matches}`]);
  } else {
    rows.push(['Solve', '—']);
  }
  rows.push(['PA samples', String(s.pa.num_samples) +
    (s.pa.arc_deg != null ? `  (arc ${fmt(s.pa.arc_deg,1)}°)` : '')]);
  if (s.pa.total_arcmin != null) {
    rows.push(['PA offset',
      `<span class="${s.pa.total_arcmin < 5 ? 'pa-good' : 'pa-bad'}">`
      + `alt ${fmt(s.pa.alt_arcmin,1)}′  az ${fmt(s.pa.az_arcmin,1)}′  `
      + `total ${fmt(s.pa.total_arcmin,2)}′</span>`]);
  }
  grid.innerHTML = rows.map(([k, v]) =>
    `<b>${k}</b><span>${v}</span><b></b><span></span>`).join('');
}
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# HTTP / WS handlers

async def index(request):
    return web.Response(text=INDEX_HTML, content_type="text/html")


async def sky_png(request):
    sim = request.app["sim"]
    return web.Response(body=sim.image_png, content_type="image/png",
                         headers={"Cache-Control": "no-store"})


async def ws_handler(request):
    sim = request.app["sim"]
    ws = web.WebSocketResponse()
    await ws.prepare(request)
    sim.subscribers.add(ws)
    # Initial snapshot.
    await ws.send_json({"type": "tick", "state": sim.public_state()})
    try:
        async for msg in ws:
            if msg.type == WSMsgType.TEXT:
                try:
                    await sim.handle_cmd(json.loads(msg.data))
                except Exception as e:
                    print(f"cmd error: {e}", file=sys.stderr)
            elif msg.type == WSMsgType.ERROR:
                break
    finally:
        sim.subscribers.discard(ws)
    return ws


async def init_app(args):
    if args.star_source == "binary":
        stars = load_binary_catalog(args.db_stars)
    else:
        stars = load_hipparcos_cached()
    sim = Sim(args, stars)
    # Initial image.
    await sim.regenerate_image()
    # Kick off the tracker reader.
    asyncio.create_task(sim.tracker_reader())

    app = web.Application()
    app["sim"] = sim
    app.router.add_get("/", index)
    app.router.add_get("/sky.png", sky_png)
    app.router.add_get("/ws", ws_handler)
    return app


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
    p.add_argument("--noise",  type=int,   default=8)
    p.add_argument("--latitude", type=float, default=-33.86)
    p.add_argument("--star-source", choices=("hipparcos", "binary"),
                   default="hipparcos")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8765)
    args = p.parse_args()

    async def runner():
        app = await init_app(args)
        runner = web.AppRunner(app)
        await runner.setup()
        site = web.TCPSite(runner, args.host, args.port)
        await site.start()
        print(f"\nSimulator running at http://{args.host}:{args.port}\n",
              file=sys.stderr)
        try:
            while True:
                await asyncio.sleep(3600)
        except asyncio.CancelledError:
            pass
        finally:
            await runner.cleanup()
            app["sim"].tracker.close()

    asyncio.run(runner())


if __name__ == "__main__":
    main()
