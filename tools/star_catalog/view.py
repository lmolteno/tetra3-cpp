#!/usr/bin/env python3
"""Live viewer for pi_tracker NDJSON output.

Renders the 128x64 SH1106 framebuffer to your terminal using Unicode
half-block characters (so 1 char = 2 vertical pixels) and prints the solve
state below it.

Usage:
    ssh pi './pi_tracker ...' | tools/star_catalog/view.py
    ./pi_tracker ... | tools/star_catalog/view.py
"""

import base64
import json
import math
import sys


# SH1106 page-major: each byte holds 8 vertical pixels of one column,
# bit 0 = top of page. Page n covers y in [n*8, n*8+7].
def decode_fb(b64):
    # Tolerate copy-paste mangling: pad up to a multiple of 4.
    pad = (-len(b64)) % 4
    try:
        raw = base64.b64decode(b64 + "=" * pad)
    except Exception:
        return None
    if len(raw) != 1024:
        return None
    rows = [[False] * 128 for _ in range(64)]
    for y in range(64):
        page = y // 8
        bit = 1 << (y % 8)
        for x in range(128):
            rows[y][x] = bool(raw[x + page * 128] & bit)
    return rows


HALF_TOP = "▀"  # ▀
HALF_BOT = "▄"  # ▄
FULL = "█"      # █


def render_fb(rows):
    out = []
    for y0 in range(0, 64, 2):
        line = []
        for x in range(128):
            t = rows[y0][x]
            b = rows[y0 + 1][x]
            line.append(FULL if (t and b) else HALF_TOP if t else HALF_BOT if b else " ")
        out.append("".join(line))
    return "\n".join(out)


def render_status(o):
    parts = [f'mode={o.get("mode", "?")}']
    if o.get("solved"):
        ra = math.degrees(o["ra_rad"])
        dec = math.degrees(o["dec_rad"])
        fov = math.degrees(o["fov_rad"])
        roll = math.degrees(o["roll_rad"])
        parts.append("\033[32msolved\033[0m")
        parts.append(f"ra={ra:7.3f}°")
        parts.append(f"dec={dec:+7.3f}°")
        parts.append(f"fov={fov:5.2f}°")
        parts.append(f"roll={roll:6.1f}°")
        parts.append(f"matches={o['num_matches']}")
        parts.append(f"rmse={o['rmse']:.1f}\"")
        parts.append(f"solve={o['solve_time_ms']:.1f}ms")
    elif o.get("solved") is False:
        parts.append("\033[31munsolved\033[0m")
    bits = []
    if "lens_pos" in o:
        bits.append(f"lens={o['lens_pos']}")
    if "exposure_us" in o:
        bits.append(f"exp={o['exposure_us']/1000:.0f}ms")
    if "gain" in o:
        bits.append(f"gain={o['gain']}")
    if bits:
        parts.append("[" + " ".join(bits) + "]")
    if pa := o.get("pa"):
        ps = []
        if "num_samples" in pa:
            ps.append(f"n={pa['num_samples']}")
        if pa.get("pole_valid"):
            ps.append(f"arc={pa.get('arc_deg', 0):.1f}°")
        if "total_arcmin" in pa:
            ps.append(f"err={pa['total_arcmin']:.1f}'")
        if ps:
            parts.append("(pa " + " ".join(ps) + ")")
    if "frame_path" in o:
        parts.append(o["frame_path"])
    out = ["  ".join(parts)]
    if t := o.get("timing_ms"):
        order = ["capture", "wait", "bin", "archive", "load",
                 "detect", "bg", "mask", "cluster", "match",
                 "solve_wall", "result_age"]
        bits = [f"{k}={t[k]:.0f}" for k in order if k in t]
        if bits:
            out.append("  timing_ms: " + " ".join(bits))
    return "\n".join(out)


def main():
    sys.stdout.write("\033[2J\033[?25l")  # clear screen, hide cursor
    try:
        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue
            try:
                o = json.loads(line)
            except json.JSONDecodeError:
                continue

            sys.stdout.write("\033[H")  # cursor home
            if "fb" in o:
                rows = decode_fb(o["fb"])
                if rows:
                    sys.stdout.write(render_fb(rows) + "\n")
            sys.stdout.write("\n" + render_status(o) + "\n")
            sys.stdout.write("\033[J")  # clear to end of screen
            sys.stdout.flush()
    except KeyboardInterrupt:
        pass
    finally:
        sys.stdout.write("\033[?25h")  # show cursor


if __name__ == "__main__":
    main()
