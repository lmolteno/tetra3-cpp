#pragma once

#ifdef HAS_SDL_DISPLAY

#include "display.h"
#include <SDL2/SDL.h>
#include <array>

class SdlDisplay : public IDisplay {
public:
    SdlDisplay(int scale = 4);
    ~SdlDisplay();

    void clear() override;
    void draw_text(int x, int y, const std::string &text, bool large = false) override;
    void draw_pixel(int x, int y) override;
    void update() override;
    const uint8_t *get_framebuffer() const override { return framebuffer.data(); }

private:
    int scale;
    SDL_Window *window = nullptr;
    SDL_Renderer *renderer = nullptr;
    SDL_Texture *texture = nullptr;
    std::array<uint8_t, FB_SIZE> framebuffer{};

    void set_fb_pixel(int x, int y);
    bool get_fb_pixel(int x, int y) const;
    void draw_char(int x, int y, char c, bool large);
};

#endif // HAS_SDL_DISPLAY
