#pragma once

#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <unistd.h>
#include <poll.h>

// Sysfs-based GPIO input. Active-low: button shorts the pin to GND;
// pin is HIGH at rest via an external or SoC pull-up.
//
// NOTE: the Pi's BCM GPIO has configurable bias but sysfs can't set it.
// Either wire a physical pull-up, or add `gpio=<N>=pu` to /boot/config.txt
// (or run `raspi-gpio set <N> pu` before launching).
//
// Requires root or membership in the 'gpio' group to export sysfs pins.
class GpioButton {
public:
    explicit GpioButton(int pin, int debounce_ms = 50)
        : pin_(pin), debounce_ms_(debounce_ms)
    {
        sysfs_write("/sys/class/gpio/export", std::to_string(pin));

        std::string base = "/sys/class/gpio/gpio" + std::to_string(pin);
        sysfs_write(base + "/direction", "in");
        sysfs_write(base + "/edge",      "falling");  // active-low

        std::string vpath = base + "/value";
        fd_ = ::open(vpath.c_str(), O_RDONLY | O_NONBLOCK);
        if (fd_ < 0)
            throw std::runtime_error("GpioButton: open " + vpath + " failed");

        // Discard any pending interrupt on open.
        char buf[4];
        (void)::read(fd_, buf, sizeof(buf));
    }

    ~GpioButton() {
        if (fd_ >= 0) ::close(fd_);
        sysfs_write("/sys/class/gpio/unexport", std::to_string(pin_));
    }

    GpioButton(const GpioButton &) = delete;
    GpioButton &operator=(const GpioButton &) = delete;

    // Non-blocking (timeout_ms=0) or blocking poll.
    // Returns true if a debounced falling edge was detected.
    bool poll_press(int timeout_ms = 0) {
        struct pollfd pfd{fd_, POLLPRI | POLLERR, 0};
        if (::poll(&pfd, 1, timeout_ms) <= 0) return false;

        char buf[4]{};
        ::lseek(fd_, 0, SEEK_SET);
        (void)::read(fd_, buf, sizeof(buf));
        if (buf[0] != '0') return false;  // '0' = LOW = pressed

        // Debounce: wait and confirm still pressed
        ::usleep(static_cast<useconds_t>(debounce_ms_) * 1000);
        ::lseek(fd_, 0, SEEK_SET);
        (void)::read(fd_, buf, sizeof(buf));
        return (buf[0] == '0');
    }

private:
    int pin_;
    int debounce_ms_;
    int fd_ = -1;

    static void sysfs_write(const std::string &path, const std::string &val) {
        int fd = ::open(path.c_str(), O_WRONLY);
        if (fd < 0) return;  // silently ignore (e.g. already exported)
        (void)::write(fd, val.c_str(), val.size());
        ::close(fd);
    }
};
