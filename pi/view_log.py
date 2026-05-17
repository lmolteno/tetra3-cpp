#!/usr/bin/env python3
"""
Live viewer for pi_tracker JSON log via SSH.
Usage: python3 view_log.py [host]  (default: 192.168.31.237)
"""

import sys
import subprocess
import json
import base64
import math

HOST = sys.argv[1] if len(sys.argv) > 1 else "192.168.31.237"
LOG  = "/var/log/pi_tracker.log"

W, H = 128, 64

# ── ANSI helpers ──────────────────────────────────────────────────────────────
ESC   = "\033"
RESET = f"{ESC}[0m"
BOLD  = f"{ESC}[1m"
DIM   = f"{ESC}[2m"
WHITE = f"{ESC}[97m"
CYAN  = f"{ESC}[96m"
GREEN = f"{ESC}[92m"
YELLOW= f"{ESC}[93m"
RED   = f"{ESC}[91m"

def clear():  print(f"{ESC}[2J{ESC}[H", end="")
def home():   print(f"{ESC}[H", end="")
def erase_down(): print(f"{ESC}[J", end="")


# ── Framebuffer decode ────────────────────────────────────────────────────────
def decode_fb(b64: str) -> list[list[bool]]:
    """Decode base64 SH1106 page-major framebuffer → pixels[y][x]."""
    data = base64.b64decode(b64)
    px = [[False] * W for _ in range(H)]
    for page in range(H // 8):
        for x in range(W):
            byte = data[x + page * W]
            for bit in range(8):
                px[page * 8 + bit][x] = bool(byte & (1 << bit))
    return px


def render_fb(px: list[list[bool]]) -> list[str]:
    """Render 128×64 pixel grid as 128×32 half-block terminal lines."""
    lines = []
    for row in range(0, H, 2):
        chars = []
        for x in range(W):
            top = px[row][x]
            bot = px[row + 1][x] if row + 1 < H else False
            if top and bot:
                chars.append("█")
            elif top:
                chars.append("▀")
            elif bot:
                chars.append("▄")
            else:
                chars.append(" ")
        lines.append("".join(chars))
    return lines


# ── Detail formatting ─────────────────────────────────────────────────────────
def deg(r): return math.degrees(r)

def fmt_timing(t: dict) -> str:
    keys = ["capture", "wait", "bin", "detect", "solve_wall", "result_age"]
    parts = [f"{k}={t[k]:.0f}" for k in keys if k in t]
    return "  ".join(parts) + " ms"


def build_details(obj: dict) -> list[str]:
    lines = []
    mode   = obj.get("mode", "?")
    solved = obj.get("solved", False)

    # ── status bar ────────────────────────────────────────────────────────────
    status_col = GREEN if solved else RED
    status_str = "SOLVED" if solved else "no fix"
    exp_ms  = obj.get("exposure_us", 0) / 1000
    gain    = obj.get("gain", "?")
    lens    = obj.get("lens_pos")
    lens_s  = f"  lens={lens:.0f}" if lens is not None else ""
    lines.append(
        f"  {BOLD}mode={CYAN}{mode}{RESET}  "
        f"{status_col}{BOLD}{status_str}{RESET}  "
        f"exp={exp_ms:.0f}ms  gain={gain}{lens_s}"
    )

    # ── solve result ──────────────────────────────────────────────────────────
    if solved:
        ra   = deg(obj["ra_rad"])
        dec  = deg(obj["dec_rad"])
        roll = deg(obj["roll_rad"])
        fov  = deg(obj["fov_rad"])
        nm   = obj.get("num_matches", "?")
        rmse = obj.get("rmse", 0)
        lines.append(
            f"  {WHITE}RA={ra:7.3f}°  Dec={dec:+7.3f}°  "
            f"Roll={roll:+6.2f}°  FOV={fov:.3f}°  "
            f"matches={nm}  rmse={rmse:.3f}{RESET}"
        )

    # ── IMU pointing ──────────────────────────────────────────────────────────
    imu = obj.get("imu")
    if imu:
        ra   = deg(imu["ra_rad"])
        dec  = deg(imu["dec_rad"])
        roll = deg(imu["roll_rad"])
        alt  = deg(imu.get("alt_rad", 0))
        az   = deg(imu.get("az_rad", 0))
        if az < 0: az += 360.0
        cal  = imu.get("calibrated", False)
        tag  = f"{GREEN}cal{RESET}" if cal else f"{YELLOW}raw{RESET}"
        col  = CYAN if cal else YELLOW
        lines.append(
            f"  {col}IMU [{tag}{col}]  "
            f"alt={alt:+6.2f}°  az={az:6.2f}°    "
            f"RA={ra:7.3f}°  Dec={dec:+7.3f}°  Roll={roll:+6.2f}°{RESET}"
        )

    # ── polar align ───────────────────────────────────────────────────────────
    pa = obj.get("pa")
    if pa:
        ns = pa.get("num_samples", 0)
        lines.append(f"  {YELLOW}PA samples={ns}{RESET}", )
        if pa.get("pole_valid"):
            pra  = deg(pa["pole_ra_rad"])
            pdec = deg(pa["pole_dec_rad"])
            arc  = pa.get("arc_deg", 0)
            lines.append(
                f"  {YELLOW}pole RA={pra:.3f}° Dec={pdec:.3f}°  "
                f"offset={arc*60:.1f} arcmin{RESET}"
            )
        if "alt_arcmin" in pa:
            alt = pa["alt_arcmin"]
            az  = pa["az_arcmin"]
            tot = pa["total_arcmin"]
            lines.append(
                f"  {YELLOW}alt={alt:+.1f}'  az={az:+.1f}'  "
                f"total={tot:.1f}'{RESET}"
            )

    # ── timing ────────────────────────────────────────────────────────────────
    if "timing_ms" in obj:
        lines.append(f"  {DIM}{fmt_timing(obj['timing_ms'])}{RESET}")

    return lines


# ── Main render ───────────────────────────────────────────────────────────────
BORDER_TOP    = "┌" + "─" * W + "┐"
BORDER_BOTTOM = "└" + "─" * W + "┘"

def render(obj: dict):
    home()

    # OLED display box
    print(DIM + BORDER_TOP + RESET)
    if "fb" in obj:
        fb_lines = render_fb(decode_fb(obj["fb"]))
        for l in fb_lines:
            print(f"{DIM}│{RESET}{WHITE}{l}{RESET}{DIM}│{RESET}")
    else:
        for _ in range(H // 2):
            print(f"{DIM}│{' ' * W}│{RESET}")
    print(DIM + BORDER_BOTTOM + RESET)

    # Details below
    for l in build_details(obj):
        print(l)

    erase_down()
    sys.stdout.flush()


# ── Entry point ───────────────────────────────────────────────────────────────
def main():
    clear()
    print(f"Connecting to {HOST}:{LOG} …")

    cmd = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
           HOST, f"tail -f -n 1 {LOG}"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    try:
        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            render(obj)
    except KeyboardInterrupt:
        pass
    finally:
        proc.terminate()
        print(RESET)


if __name__ == "__main__":
    main()
