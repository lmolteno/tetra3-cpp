#ifdef HAS_SDL_DISPLAY

#include "display_sdl.h"
#include "font_6x8.h"
#include <iostream>
#include <cstring>

SdlDisplay::SdlDisplay(int scale) : scale(scale) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL init failed: " << SDL_GetError() << std::endl;
        return;
    }

    window = SDL_CreateWindow("OLED Emulator",
                              SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                              WIDTH * scale, HEIGHT * scale,
                              SDL_WINDOW_SHOWN);
    if (!window) {
        std::cerr << "Window creation failed: " << SDL_GetError() << std::endl;
        return;
    }

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        std::cerr << "Renderer creation failed: " << SDL_GetError() << std::endl;
        return;
    }

    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB332,
                                SDL_TEXTUREACCESS_STREAMING,
                                WIDTH, HEIGHT);

    clear();
}

SdlDisplay::~SdlDisplay() {
    if (texture) SDL_DestroyTexture(texture);
    if (renderer) SDL_DestroyRenderer(renderer);
    if (window) SDL_DestroyWindow(window);
    SDL_Quit();
}

void SdlDisplay::set_fb_pixel(int x, int y) {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return;
    int byte_idx = x + (y / 8) * WIDTH;
    framebuffer[byte_idx] |= (1 << (y % 8));
}

bool SdlDisplay::get_fb_pixel(int x, int y) const {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return false;
    int byte_idx = x + (y / 8) * WIDTH;
    return (framebuffer[byte_idx] >> (y % 8)) & 1;
}

void SdlDisplay::clear() {
    framebuffer.fill(0);
}

void SdlDisplay::draw_char(int x, int y, char c, bool large) {
    int idx;
    if (static_cast<unsigned char>(c) == 0xB0) {
        idx = 95; // degree symbol at slot 127 (127-32)
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
                        set_fb_pixel(x + col * sx + dx, y + row * sy + dy);
                    }
                }
            }
        }
    }
}

void SdlDisplay::draw_text(int x, int y, const std::string &text, bool large) {
    int char_w = large ? 12 : 6;
    for (size_t i = 0; i < text.size(); i++) {
        draw_char(x + i * char_w, y, text[i], large);
    }
}

void SdlDisplay::draw_pixel(int x, int y) {
    set_fb_pixel(x, y);
}

void SdlDisplay::update() {
    if (!renderer || !texture) return;

    // Pump events to keep window responsive
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_QUIT) {
            // Could set a flag here
        }
    }

    // Convert 1-bit framebuffer to RGB332 texture
    uint8_t pixels[WIDTH * HEIGHT];
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            pixels[y * WIDTH + x] = get_fb_pixel(x, y) ? 0xFF : 0x00;
        }
    }

    SDL_UpdateTexture(texture, nullptr, pixels, WIDTH);
    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    SDL_RenderPresent(renderer);
}

#endif // HAS_SDL_DISPLAY
