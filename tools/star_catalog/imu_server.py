"""IMU WebSocket/HTTP server.

Serves the phone UI and receives orientation quaternions via WebSocket.
Designed to run inside Textual's asyncio event loop.
"""

import asyncio
import json
import logging
import socket

from aiohttp import web, WSMsgType

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Phone HTML — embedded as string so no extra files needed
# ---------------------------------------------------------------------------

PHONE_HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=no">
<title>tetra3 IMU</title>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: -apple-system, system-ui, sans-serif;
    background: #1a1a2e; color: #e0e0e0;
    padding: 12px; max-width: 540px; margin: 0 auto;
    -webkit-user-select: none; user-select: none;
  }
  h1 { font-size: 18px; text-align: center; margin-bottom: 8px; color: #7fbbff; }

  .status-bar {
    display: flex; gap: 8px; justify-content: center;
    margin-bottom: 10px; font-size: 13px;
  }
  .status-dot {
    width: 10px; height: 10px; border-radius: 50%;
    display: inline-block; margin-right: 4px; vertical-align: middle;
  }
  .dot-ok { background: #4caf50; }
  .dot-err { background: #f44336; }
  .dot-wait { background: #ff9800; }

  #oled-canvas {
    display: block; margin: 0 auto 10px;
    border: 1px solid #444; border-radius: 4px;
    image-rendering: pixelated;
    width: 100%; max-width: 512px; aspect-ratio: 2/1;
  }

  .readout {
    display: grid; grid-template-columns: 1fr 1fr 1fr;
    gap: 4px; margin-bottom: 10px; text-align: center; font-size: 13px;
  }
  .readout div { background: #16213e; border-radius: 4px; padding: 6px 4px; }
  .readout .label { font-size: 10px; color: #888; }
  .readout .value { font-size: 16px; font-weight: 600; font-variant-numeric: tabular-nums; }

  .controls { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 10px; justify-content: center; }
  .controls button {
    padding: 8px 14px; border: 1px solid #444; border-radius: 6px;
    background: #16213e; color: #e0e0e0; font-size: 13px; cursor: pointer;
    transition: background 0.15s;
  }
  .controls button:active, .controls button.active { background: #0f3460; border-color: #7fbbff; }

  .fov-row {
    display: flex; align-items: center; gap: 8px;
    margin-bottom: 10px; padding: 0 8px;
  }
  .fov-row label { font-size: 13px; white-space: nowrap; }
  .fov-row input[type=range] { flex: 1; }
  .fov-row .fov-val { font-size: 14px; font-weight: 600; width: 50px; text-align: right; }

  .solve-info {
    font-size: 12px; color: #888; text-align: center;
    margin-bottom: 8px; min-height: 32px;
  }

  .perm-btn {
    display: block; width: 100%; padding: 14px;
    background: #0f3460; color: white; border: none; border-radius: 8px;
    font-size: 16px; cursor: pointer; margin-bottom: 10px;
  }
  .perm-btn:active { background: #1a5276; }
  .hidden { display: none; }
</style>
</head>
<body>

<h1>tetra3 IMU</h1>

<button id="perm-btn" class="perm-btn">Tap to enable sensors</button>

<div id="main-ui" class="hidden">
  <div class="status-bar">
    <span><span id="ws-dot" class="status-dot dot-wait"></span>Server</span>
    <span><span id="imu-dot" class="status-dot dot-wait"></span>IMU</span>
    <span><span id="gps-dot" class="status-dot dot-wait"></span>GPS</span>
  </div>

  <canvas id="oled-canvas" width="128" height="64"></canvas>

  <div class="readout">
    <div><div class="label">RA</div><div class="value" id="ra-val">—</div></div>
    <div><div class="label">Dec</div><div class="value" id="dec-val">—</div></div>
    <div><div class="label">Roll</div><div class="value" id="roll-val">—</div></div>
  </div>
  <div class="readout">
    <div><div class="label">Alt</div><div class="value" id="alt-val">—</div></div>
    <div><div class="label">Az</div><div class="value" id="az-val">—</div></div>
    <div><div class="label">α/β/γ</div><div class="value" id="euler-val" style="font-size:11px">—</div></div>
  </div>
  <div class="readout" id="gps-readout">
    <div><div class="label">GPS</div><div class="value" id="gps-val" style="font-size:13px">waiting...</div></div>
    <div style="grid-column: span 2; display:flex; align-items:center; gap:4px; justify-content:center">
      <input type="number" id="manual-lat" placeholder="Lat" step="0.01" style="width:80px; background:#0d1b2a; border:1px solid #444; border-radius:4px; color:#e0e0e0; padding:4px; font-size:13px">
      <input type="number" id="manual-lon" placeholder="Lon" step="0.01" style="width:80px; background:#0d1b2a; border:1px solid #444; border-radius:4px; color:#e0e0e0; padding:4px; font-size:13px">
      <button id="btn-set-gps" style="padding:4px 8px; background:#0f3460; border:1px solid #444; border-radius:4px; color:#e0e0e0; font-size:12px">Set</button>
    </div>
  </div>

  <div class="fov-row">
    <label>FOV</label>
    <input type="range" id="fov-slider" min="10" max="150" value="75" step="1">
    <span class="fov-val" id="fov-val">75°</span>
  </div>

  <div class="controls">
    <button id="btn-solve">Solve Now</button>
    <button id="btn-tracking" class="active">Tracking</button>
    <button id="btn-pa-sample">PA Sample</button>
    <button id="btn-pa-fix">PA Fix</button>
  </div>

  <div class="solve-info" id="solve-info">Waiting for solver...</div>
</div>

<script>
const DEG2RAD = Math.PI / 180;
const RAD2DEG = 180 / Math.PI;

let ws = null;
let lat = null, lon = null;
let imuActive = false;
let sendInterval = null;
let lastQuat = [1, 0, 0, 0];

// --- Canvas ---
const canvas = document.getElementById('oled-canvas');
const ctx = canvas.getContext('2d');
ctx.imageSmoothingEnabled = false;

function renderFB(b64) {
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  if (bytes.length < 1024) return;
  const imgData = ctx.createImageData(128, 64);
  for (let y = 0; y < 64; y++) {
    for (let x = 0; x < 128; x++) {
      const byteIdx = x + Math.floor(y / 8) * 128;
      const bit = (bytes[byteIdx] >> (y % 8)) & 1;
      const px = (y * 128 + x) * 4;
      imgData.data[px] = imgData.data[px+1] = imgData.data[px+2] = bit ? 255 : 0;
      imgData.data[px+3] = 255;
    }
  }
  ctx.putImageData(imgData, 0, 0);
}

// --- Quaternion from Euler (ZXY intrinsic) ---
function eulerToQuaternion(alpha, beta, gamma) {
  const a = alpha * DEG2RAD / 2;
  const b = beta * DEG2RAD / 2;
  const g = gamma * DEG2RAD / 2;
  const ca = Math.cos(a), sa = Math.sin(a);
  const cb = Math.cos(b), sb = Math.sin(b);
  const cg = Math.cos(g), sg = Math.sin(g);
  return [
    ca*cb*cg - sa*sb*sg,  // w
    ca*sb*cg - sa*cb*sg,  // x
    ca*cb*sg + sa*sb*cg,  // y
    ca*sb*sg + sa*cb*cg,  // z
  ];
}

// --- Status dots ---
function setDot(id, state) {
  const el = document.getElementById(id);
  el.className = 'status-dot ' + ({ok:'dot-ok', err:'dot-err', wait:'dot-wait'}[state] || 'dot-wait');
}

// --- WebSocket ---
function connectWS() {
  const proto = location.protocol === 'https:' ? 'wss:' : 'ws:';
  ws = new WebSocket(proto + '//' + location.host + '/ws');
  ws.onopen = () => { setDot('ws-dot', 'ok'); startSending(); };
  ws.onclose = () => { setDot('ws-dot', 'err'); stopSending(); setTimeout(connectWS, 2000); };
  ws.onerror = () => { ws.close(); };
  ws.onmessage = (ev) => {
    try {
      const msg = JSON.parse(ev.data);
      if (msg.fb) renderFB(msg.fb);
      if (msg.ra !== undefined) document.getElementById('ra-val').textContent = msg.ra.toFixed(2) + '°';
      if (msg.dec !== undefined) document.getElementById('dec-val').textContent = msg.dec.toFixed(2) + '°';
      if (msg.roll !== undefined) document.getElementById('roll-val').textContent = msg.roll.toFixed(1) + '°';
      if (msg.alt !== undefined) document.getElementById('alt-val').textContent = msg.alt.toFixed(1) + '°';
      if (msg.az !== undefined) document.getElementById('az-val').textContent = msg.az.toFixed(1) + '°';
      if (msg.fov !== undefined) {
        document.getElementById('fov-slider').value = msg.fov;
        document.getElementById('fov-val').textContent = Math.round(msg.fov) + '°';
      }
      if (msg.solve) {
        const s = msg.solve;
        if (s.solved) {
          document.getElementById('solve-info').textContent =
            'Solved: RA=' + s.ra_deg.toFixed(2) + ' Dec=' + s.dec_deg.toFixed(2) +
            ' RMSE=' + (s.rmse||0).toFixed(1) + '" (' + (s.solve_time_ms||0).toFixed(0) + 'ms)';
        } else {
          document.getElementById('solve-info').textContent = 'No solution (' + (s.solve_time_ms||0).toFixed(0) + 'ms)';
        }
      }
    } catch(e) {}
  };
}

function startSending() {
  if (sendInterval) return;
  sendInterval = setInterval(() => {
    if (ws && ws.readyState === 1 && imuActive) {
      const msg = { q: lastQuat };
      if (lat !== null) { msg.lat = lat; msg.lon = lon; }
      ws.send(JSON.stringify(msg));
    }
  }, 50); // 20 Hz
}

function stopSending() {
  if (sendInterval) { clearInterval(sendInterval); sendInterval = null; }
}

// --- DeviceOrientation ---
function onOrientation(ev) {
  if (ev.alpha === null) return;
  imuActive = true;
  setDot('imu-dot', 'ok');
  lastQuat = eulerToQuaternion(ev.alpha, ev.beta, ev.gamma);
  document.getElementById('euler-val').textContent =
    ev.alpha.toFixed(0) + '/' + ev.beta.toFixed(0) + '/' + ev.gamma.toFixed(0);
}

// --- Geolocation ---
function setGPS(newLat, newLon, source) {
  lat = newLat; lon = newLon;
  setDot('gps-dot', 'ok');
  document.getElementById('gps-val').textContent = lat.toFixed(4) + ', ' + lon.toFixed(4) + ' (' + source + ')';
}

function startGPS() {
  if (!navigator.geolocation || location.protocol === 'http:') {
    setDot('gps-dot', 'err');
    document.getElementById('gps-val').textContent = 'needs HTTPS — use manual';
    return;
  }
  navigator.geolocation.watchPosition(
    (pos) => { setGPS(pos.coords.latitude, pos.coords.longitude, 'gps'); },
    () => { setDot('gps-dot', 'err'); document.getElementById('gps-val').textContent = 'denied — use manual'; },
    { enableHighAccuracy: false, maximumAge: 30000 }
  );
}

document.getElementById('btn-set-gps').addEventListener('click', () => {
  const mLat = parseFloat(document.getElementById('manual-lat').value);
  const mLon = parseFloat(document.getElementById('manual-lon').value);
  if (!isNaN(mLat) && !isNaN(mLon)) {
    setGPS(mLat, mLon, 'manual');
  }
});

// --- Permission button ---
document.getElementById('perm-btn').addEventListener('click', async () => {
  // iOS requires explicit permission request
  if (typeof DeviceOrientationEvent !== 'undefined' && typeof DeviceOrientationEvent.requestPermission === 'function') {
    try {
      const state = await DeviceOrientationEvent.requestPermission();
      if (state !== 'granted') { alert('Orientation permission denied'); return; }
    } catch(e) { alert('Permission error: ' + e); return; }
  }
  window.addEventListener('deviceorientation', onOrientation);
  startGPS();
  document.getElementById('perm-btn').classList.add('hidden');
  document.getElementById('main-ui').classList.remove('hidden');
  connectWS();
});

// --- Controls ---
document.getElementById('btn-solve').addEventListener('click', () => {
  if (ws && ws.readyState === 1) ws.send(JSON.stringify({cmd: 'solve'}));
});

['tracking', 'pa-sample', 'pa-fix'].forEach(mode => {
  document.getElementById('btn-' + mode).addEventListener('click', (ev) => {
    document.querySelectorAll('.controls button').forEach(b => b.classList.remove('active'));
    ev.target.classList.add('active');
    const modeMap = {'tracking':'tracking', 'pa-sample':'pa_sampling', 'pa-fix':'pa_fix'};
    if (ws && ws.readyState === 1) ws.send(JSON.stringify({cmd: 'mode', value: modeMap[mode]}));
  });
});

const fovSlider = document.getElementById('fov-slider');
fovSlider.addEventListener('input', () => {
  document.getElementById('fov-val').textContent = fovSlider.value + '°';
});
fovSlider.addEventListener('change', () => {
  if (ws && ws.readyState === 1) ws.send(JSON.stringify({cmd: 'fov', value: parseFloat(fovSlider.value)}));
});
</script>
</body>
</html>
"""


def get_local_ip():
    """Get the local IP address reachable by devices on the same network."""
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"


class IMUServer:
    """aiohttp server: serves phone UI, receives IMU quaternions via WebSocket."""

    def __init__(self, on_orientation_update, on_command=None, port=8080):
        self.on_orientation_update = on_orientation_update
        self.on_command = on_command
        self.port = port
        self._app = None
        self._runner = None
        self._ws_clients = set()

    async def start(self):
        """Start the aiohttp server."""
        self._app = web.Application()
        self._app.router.add_get("/", self.handle_index)
        self._app.router.add_get("/ws", self.handle_websocket)

        self._runner = web.AppRunner(self._app)
        await self._runner.setup()
        site = web.TCPSite(self._runner, "0.0.0.0", self.port)
        await site.start()

    async def stop(self):
        """Stop the server and close all WebSocket connections."""
        for ws in list(self._ws_clients):
            await ws.close()
        self._ws_clients.clear()
        if self._runner:
            await self._runner.cleanup()

    async def handle_index(self, request):
        """Serve the embedded phone HTML."""
        return web.Response(text=PHONE_HTML, content_type="text/html")

    async def handle_websocket(self, request):
        """Handle WebSocket connections from the phone."""
        ws = web.WebSocketResponse()
        await ws.prepare(request)
        self._ws_clients.add(ws)
        log.info("Phone connected via WebSocket")

        try:
            async for msg in ws:
                if msg.type == WSMsgType.TEXT:
                    try:
                        data = json.loads(msg.data)
                    except json.JSONDecodeError:
                        continue

                    if "cmd" in data:
                        # Command message
                        if self.on_command:
                            await self.on_command(data["cmd"], data.get("value"))
                    elif "q" in data:
                        # Orientation update
                        q = data["q"]
                        lat = data.get("lat")
                        lon = data.get("lon")
                        if isinstance(q, list) and len(q) == 4:
                            await self.on_orientation_update(q, lat, lon)
                elif msg.type == WSMsgType.ERROR:
                    log.warning("WebSocket error: %s", ws.exception())
        finally:
            self._ws_clients.discard(ws)
            log.info("Phone disconnected")

        return ws

    async def broadcast_state(self, fb_b64, ra_deg, dec_deg, roll_deg,
                              fov_deg=None, alt_deg=None, az_deg=None,
                              solve_status=None):
        """Push framebuffer + coords + solve status to all connected phone clients."""
        msg = {
            "fb": fb_b64,
            "ra": ra_deg,
            "dec": dec_deg,
            "roll": roll_deg,
        }
        if fov_deg is not None:
            msg["fov"] = fov_deg
        if alt_deg is not None:
            msg["alt"] = alt_deg
        if az_deg is not None:
            msg["az"] = az_deg
        if solve_status is not None:
            msg["solve"] = solve_status

        payload = json.dumps(msg)
        dead = set()
        for ws in self._ws_clients:
            try:
                await ws.send_str(payload)
            except Exception:
                dead.add(ws)
        self._ws_clients -= dead

    @property
    def client_count(self):
        return len(self._ws_clients)

    @property
    def url(self):
        ip = get_local_ip()
        return f"http://{ip}:{self.port}"
