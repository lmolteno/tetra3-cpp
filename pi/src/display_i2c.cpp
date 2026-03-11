#include "display_i2c.h"
#include "font_6x8.h"
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/i2c-dev.h>

I2cSh1106Display::I2cSh1106Display(const char *device, uint8_t address) : addr(address) {
    fd = open(device, O_RDWR);
    if (fd < 0) {
        std::cerr << "Failed to open I2C device: " << device << std::endl;
        return;
    }

    if (ioctl(fd, I2C_SLAVE, addr) < 0) {
        std::cerr << "Failed to set I2C address 0x" << std::hex << (int)addr << std::endl;
        close(fd);
        fd = -1;
        return;
    }

    // SH1106 init sequence
    const uint8_t init_cmds[] = {
        0xAE,       // Display OFF
        0xD5, 0x80, // Set display clock divide
        0xA8, 0x3F, // Set multiplex ratio (64)
        0xD3, 0x00, // Set display offset
        0x40,       // Set start line
        0xAD, 0x8B, // Set charge pump (internal VCC)
        0xA1,       // Segment remap (mirror X)
        0xC8,       // COM output scan direction (mirror Y)
        0xDA, 0x12, // Set COM pins config
        0x81, 0xFF, // Set contrast (max)
        0xD9, 0x1F, // Set pre-charge period
        0xDB, 0x40, // Set VCOMH deselect level
        0x33,       // Set VPP to 9V
        0xA6,       // Normal display (not inverted)
        0xAF,       // Display ON
    };

    for (uint8_t cmd : init_cmds) {
        if (!send_command(cmd)) {
            std::cerr << "SH1106 init failed at command 0x" << std::hex << (int)cmd << std::endl;
            close(fd);
            fd = -1;
            return;
        }
    }

    initialized = true;
    clear();
    update();
    std::cout << "SH1106 OLED initialized on " << device << " at 0x" << std::hex << (int)addr << std::dec << std::endl;
}

I2cSh1106Display::~I2cSh1106Display() {
    if (fd >= 0) {
        // Display off
        send_command(0xAE);
        close(fd);
    }
}

bool I2cSh1106Display::send_command(uint8_t cmd) {
    uint8_t buf[2] = {0x00, cmd}; // Co=0, D/C#=0 (command)
    return write(fd, buf, 2) == 2;
}

bool I2cSh1106Display::send_data(const uint8_t *data, size_t len) {
    // Send in chunks with 0x40 prefix (data mode)
    std::vector<uint8_t> buf(len + 1);
    buf[0] = 0x40; // Co=0, D/C#=1 (data)
    std::memcpy(&buf[1], data, len);
    return write(fd, buf.data(), buf.size()) == static_cast<ssize_t>(buf.size());
}

void I2cSh1106Display::set_fb_pixel(int x, int y) {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return;
    int byte_idx = x + (y / 8) * WIDTH;
    framebuffer[byte_idx] |= (1 << (y % 8));
}

void I2cSh1106Display::clear() {
    framebuffer.fill(0);
}

void I2cSh1106Display::draw_char(int x, int y, char c, bool large) {
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

void I2cSh1106Display::draw_text(int x, int y, const std::string &text, bool large) {
    int char_w = large ? 12 : 6;
    for (size_t i = 0; i < text.size(); i++) {
        draw_char(x + i * char_w, y, text[i], large);
    }
}

void I2cSh1106Display::draw_pixel(int x, int y) {
    set_fb_pixel(x, y);
}

void I2cSh1106Display::update() {
    if (!initialized) return;

    // SH1106 uses page addressing with 2-pixel column offset
    for (int page = 0; page < 8; page++) {
        send_command(0xB0 + page);           // Set page address
        send_command(0x02);                   // Set lower column address (offset 2)
        send_command(0x10);                   // Set higher column address
        send_data(&framebuffer[page * WIDTH], WIDTH);
    }
}
