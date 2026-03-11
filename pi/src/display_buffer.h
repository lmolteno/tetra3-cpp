#pragma once

#include "display.h"
#include "font_6x8.h"
#include <array>
#include <cstring>

// Headless display with readable framebuffer. Same pixel format as SdlDisplay.
class BufferDisplay : public IDisplay {
public:
    void clear() override { fb_.fill(0); }

    void draw_pixel(int x, int y) override {
        if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return;
        fb_[x + (y / 8) * WIDTH] |= (1 << (y % 8));
    }

    void draw_text(int x, int y, const std::string &text, bool large = false) override {
        int char_w = large ? 12 : 6;
        for (size_t i = 0; i < text.size(); i++) {
            draw_char(x + static_cast<int>(i) * char_w, y, text[i], large);
        }
    }

    void update() override {} // nothing to do

    const uint8_t *get_framebuffer() const override { return fb_.data(); }

private:
    std::array<uint8_t, FB_SIZE> fb_{};

    void draw_char(int x, int y, char c, bool large) {
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
                            draw_pixel(x + col * sx + dx, y + row * sy + dy);
                        }
                    }
                }
            }
        }
    }
};
