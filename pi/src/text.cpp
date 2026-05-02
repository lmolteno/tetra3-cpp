#include "text.h"
#include "font_6x8.h"

namespace text {

static void draw_char(ICanvas &canvas, int x, int y, char c, bool large) {
    int idx;
    if (static_cast<unsigned char>(c) == 0xB0) {
        idx = 95;
    } else {
        idx = c - 32;
    }
    if (idx < 0 || idx >= 96) return;

    const uint8_t *glyph = font_6x8[idx];
    int sx = large ? 2 : 1;
    int sy = large ? 2 : 1;

    for (int col = 0; col < 6; col++) {
        uint8_t column_data = glyph[col];
        for (int row = 0; row < 8; row++) {
            if (column_data & (1 << row)) {
                for (int dy = 0; dy < sy; dy++) {
                    for (int dx = 0; dx < sx; dx++) {
                        canvas.draw_pixel(x + col * sx + dx, y + row * sy + dy);
                    }
                }
            }
        }
    }
}

void draw(ICanvas &canvas, int x, int y, const std::string &s, bool large) {
    int char_w = large ? 12 : 6;
    for (size_t i = 0; i < s.size(); i++) {
        draw_char(canvas, x + static_cast<int>(i) * char_w, y, s[i], large);
    }
}

}
