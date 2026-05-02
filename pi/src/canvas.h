#pragma once

#include <cstdint>
#include <cstddef>

// Pure pixel canvas. No knowledge of what's drawn or where pixels go after
// `update()` — implementations may push to hardware, a window, or just hold
// a framebuffer in memory.
class ICanvas {
public:
    virtual ~ICanvas() = default;

    virtual void clear() = 0;
    virtual void draw_pixel(int x, int y) = 0;
    virtual void update() = 0;

    virtual int width() const = 0;
    virtual int height() const = 0;

    // 1-bit framebuffer access for canvases that have one. Page-major
    // (SH1106-style): byte at offset (x + (y/8)*width), bit (1 << (y%8)).
    virtual const uint8_t *framebuffer() const { return nullptr; }
    virtual std::size_t framebuffer_size() const { return 0; }

    void draw_hline(int x0, int x1, int y) {
        for (int x = x0; x <= x1; x++) draw_pixel(x, y);
    }
};
