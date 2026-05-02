#pragma once

#include "render/canvas.h"
#include <array>

// Headless 128x64 1-bit canvas. Same byte layout as the SH1106 (page-major:
// each byte is 8 vertical pixels), so the framebuffer can be blitted directly
// to the OLED by an external bridge process.
class BufferCanvas : public ICanvas {
public:
    static constexpr int W = 128;
    static constexpr int H = 64;
    static constexpr std::size_t FB_SIZE = W * H / 8;  // 1024

    void clear() override { fb_.fill(0); }

    void draw_pixel(int x, int y) override {
        if (x < 0 || x >= W || y < 0 || y >= H) return;
        fb_[x + (y / 8) * W] |= (1 << (y % 8));
    }

    void update() override {}

    int width() const override { return W; }
    int height() const override { return H; }
    const uint8_t *framebuffer() const override { return fb_.data(); }
    std::size_t framebuffer_size() const override { return FB_SIZE; }

private:
    std::array<uint8_t, FB_SIZE> fb_{};
};
