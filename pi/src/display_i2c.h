#pragma once

#include "display.h"
#include <array>

class I2cSh1106Display : public IDisplay {
public:
    I2cSh1106Display(const char *device = "/dev/i2c-1", uint8_t address = 0x3C);
    ~I2cSh1106Display();

    void clear() override;
    void draw_text(int x, int y, const std::string &text, bool large = false) override;
    void draw_pixel(int x, int y) override;
    void update() override;

private:
    int fd = -1;
    uint8_t addr;
    std::array<uint8_t, FB_SIZE> framebuffer{};
    bool initialized = false;

    // Same font as SDL version (declared extern, defined in display_i2c.cpp)
    void set_fb_pixel(int x, int y);
    void draw_char(int x, int y, char c, bool large);
    bool send_command(uint8_t cmd);
    bool send_data(const uint8_t *data, size_t len);
};
