#!/usr/bin/env python3
"""Live viewer for pi_tracker NDJSON output.

Renders the 128x64 SH1106 framebuffer to your terminal as Unicode braille
(2x4 pixels per char, so 64 cols × 16 rows — fits in any terminal and
gets the OLED's 2:1 aspect ratio right).

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


# Braille pattern: 2x4 pixels per char, base codepoint U+2800.
# Bit position by (x_in_block, y_in_block):
#     (0,0)=0   (1,0)=3
#     (0,1)=1   (1,1)=4
#     (0,2)=2   (1,2)=5
#     (0,3)=6   (1,3)=7
_BRAILLE_BIT = [
    [0, 3],
    [1, 4],
    [2, 5],
    [6, 7],
]


def render_fb_lines(rows):
    out = []
    for y0 in range(0, 64, 4):
        line = []
        for x0 in range(0, 128, 2):
            bits = 0
            for dy in range(4):
                y = y0 + dy
                if y >= 64:
                    continue
                for dx in range(2):
                    x = x0 + dx
                    if x >= 128:
                        continue
                    if rows[y][x]:
                        bits |= 1 << _BRAILLE_BIT[dy][dx]
            line.append(chr(0x2800 + bits))
        out.append("".join(line))
    return out


def render_status_lines(o):
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
        order = ["capture", "wait", "bin", "frame_age", "archive", "load",
                 "detect", "bg", "mask", "cluster", "match",
                 "solve_wall", "result_age"]
        bits = [f"{k}={t[k]:.0f}" for k in order if k in t]
        if bits:
            out.append("  timing: " + " ".join(bits))
    return out


# Alternate screen buffer — the viewer redraws in place without polluting
# the user's terminal scrollback. We use absolute cursor positioning per row
# (no trailing \n's) so a small terminal can't accidentally scroll the alt
# screen between ticks, and clear-to-EOL on each line so a shorter line
# doesn't leave residue from the previous tick.
ENTER_ALT   = "\033[?1049h"
EXIT_ALT    = "\033[?1049l"
HIDE_CURSOR = "\033[?25l"
SHOW_CURSOR = "\033[?25h"
CLEAR_EOL   = "\033[K"
CLEAR_BELOW = "\033[J"


def cursor_to(row, col=1):  # 1-indexed
    return f"\033[{row};{col}H"


def main():
    sys.stdout.write(ENTER_ALT + HIDE_CURSOR + cursor_to(1) + "\033[2J")
    sys.stdout.flush()
    try:
        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue
            try:
                o = json.loads(line)
            except json.JSONDecodeError:
                continue

            buf = []
            row = 1
            if "fb" in o:
                rows = decode_fb(o["fb"])
                if rows:
                    for fb_line in render_fb_lines(rows):
                        buf.append(cursor_to(row))
                        buf.append(fb_line)
                        buf.append(CLEAR_EOL)
                        row += 1
            row += 1  # blank separator
            for status_line in render_status_lines(o):
                buf.append(cursor_to(row))
                buf.append(status_line)
                buf.append(CLEAR_EOL)
                row += 1
            # Clear any rows below the last status line that may have content
            # from a previous tick (e.g. timing line that's now absent).
            buf.append(cursor_to(row))
            buf.append(CLEAR_BELOW)
            sys.stdout.write("".join(buf))
            sys.stdout.flush()
    except KeyboardInterrupt:
        pass
    finally:
        sys.stdout.write(SHOW_CURSOR + EXIT_ALT)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
