#pragma once

#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>
// I2C_SLAVE is 0x0703 on all Linux architectures; define it here so we don't
// need linux/i2c-dev.h (not installed in Alpine musl build environments).
#ifndef I2C_SLAVE
#define I2C_SLAVE 0x0703
#endif

// Thin I2C driver for the SSD1306 128x64 OLED.
//
// Expects the same page-major framebuffer layout as BufferCanvas:
//   byte at (x + page*128), bit (1 << bit_within_page).
// That's exactly what SSD1306 horizontal addressing mode ingests, so
// flush() is a single write of the raw 1024-byte buffer.
class I2cSsd1306 {
public:
    static constexpr int W = 128;
    static constexpr int H = 64;
    static constexpr std::size_t FB_SIZE = W * H / 8;  // 1024 bytes, 8 pages

    // bus: Linux i2c device node, e.g. "/dev/i2c-1"
    // addr: 7-bit I2C address (0x3C or 0x3D depending on SA0 pin)
    explicit I2cSsd1306(const char *bus = "/dev/i2c-1", uint8_t addr = 0x3C) {
        fd_ = ::open(bus, O_RDWR);
        if (fd_ < 0)
            throw std::runtime_error(std::string("SSD1306: open ") + bus + " failed");
        if (::ioctl(fd_, I2C_SLAVE, addr) < 0)
            throw std::runtime_error("SSD1306: I2C_SLAVE ioctl failed");
        init_display();
    }

    ~I2cSsd1306() { if (fd_ >= 0) ::close(fd_); }

    I2cSsd1306(const I2cSsd1306 &) = delete;
    I2cSsd1306 &operator=(const I2cSsd1306 &) = delete;

    // Push a FB_SIZE-byte page-major framebuffer (BufferCanvas::framebuffer())
    // to the display in one bulk I2C write.
    void flush(const uint8_t *fb) {
        // Reset the column+page address window to the full 128x64 each frame.
        // In horizontal addressing mode the controller auto-wraps, but this
        // ensures a dropped write doesn't desync the pointer.
        static constexpr uint8_t addr_cmds[] = {
            0x00,
            0x21, 0x00, 0x7F,  // column address: 0 to 127
            0x22, 0x00, 0x07,  // page address:   0 to 7
        };
        write_raw(addr_cmds, sizeof(addr_cmds));

        // Control byte 0x40 = "data follows", then raw framebuffer.
        uint8_t buf[1 + FB_SIZE];
        buf[0] = 0x40;
        std::memcpy(buf + 1, fb, FB_SIZE);
        write_raw(buf, sizeof(buf));
    }

private:
    int fd_ = -1;

    void write_raw(const uint8_t *data, std::size_t len) {
        if (::write(fd_, data, len) != static_cast<ssize_t>(len))
            throw std::runtime_error("SSD1306: I2C write failed");
    }

    void cmd(std::initializer_list<uint8_t> bytes) {
        // 0x00 = control byte: Co=0 (all following are commands), D/C#=0
        uint8_t buf[32];
        buf[0] = 0x00;
        std::size_t i = 1;
        for (uint8_t b : bytes) buf[i++] = b;
        write_raw(buf, i);
    }

    void init_display() {
        cmd({0xAE});              // display off
        cmd({0xD5, 0x80});       // clock divider / oscillator frequency
        cmd({0xA8, 0x3F});       // multiplex ratio = 63 (64 rows)
        cmd({0xD3, 0x00});       // display offset = 0
        cmd({0x40});             // display start line = 0
        cmd({0x8D, 0x14});       // charge pump: enable (internal VCC)
        cmd({0x20, 0x00});       // memory addressing: horizontal — enables bulk flush
        cmd({0xA1});             // segment remap: col 127 -> SEG0 (flip X)
        cmd({0xC8});             // COM scan: row N-1 -> COM0 (flip Y)
        cmd({0xDA, 0x12});       // COM pins: alt config, no remap (for 128x64)
        cmd({0x81, 0xCF});       // contrast
        cmd({0xD9, 0xF1});       // pre-charge period
        cmd({0xDB, 0x40});       // VCOMH deselect level
        cmd({0xA4});             // output follows RAM (not all-on)
        cmd({0xA6});             // normal (non-inverted) display
        cmd({0xAF});             // display on
    }
};
