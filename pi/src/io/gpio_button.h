#pragma once

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <poll.h>

// Sysfs-based GPIO input. Active-low: button shorts the pin to GND;
// pin is HIGH at rest via an external or SoC pull-up.
//
// Takes a BCM pin number. Automatically resolves the sysfs base offset,
// which moved from 0 to 512 on Pi kernels >= 6.6.
//
// NOTE: the Pi's BCM GPIO has configurable bias but sysfs can't set it.
// Either wire a physical pull-up, or add `gpio=<N>=pu` to /boot/usercfg.txt
// (or run `raspi-gpio set <N> pu` before launching).
//
// Requires root or membership in the 'gpio' group to export sysfs pins.
class GpioButton {
public:
    explicit GpioButton(int bcm_pin, int debounce_ms = 50)
        : pin_(bcm_pin + gpio_chip_base()), debounce_ms_(debounce_ms)
    {
        sysfs_write("/sys/class/gpio/export", std::to_string(pin_));

        std::string dir = "/sys/class/gpio/gpio" + std::to_string(pin_);
        sysfs_write(dir + "/direction", "in");
        sysfs_write(dir + "/edge",      "falling");  // active-low

        std::string vpath = dir + "/value";
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
    int pin_;        // resolved sysfs pin number (bcm + chip base)
    int debounce_ms_;
    int fd_ = -1;

    // Scan /sys/class/gpio for a gpiochip labelled as the BCM pinctrl
    // driver and return its base offset.  Kernels < 6.6 use base 0;
    // kernels >= 6.6 use 512 (or similar).  Falls back to the largest
    // base found if no label match, which works fine on a Pi 2B with a
    // single GPIO controller.
    static int gpio_chip_base() {
        DIR *d = ::opendir("/sys/class/gpio");
        if (!d) return 0;

        int best_base = 0;
        struct dirent *e;
        while ((e = ::readdir(d))) {
            if (::strncmp(e->d_name, "gpiochip", 8) != 0) continue;
            std::string dir = "/sys/class/gpio/" + std::string(e->d_name);

            char buf[32]{};
            int fd = ::open((dir + "/base").c_str(), O_RDONLY);
            if (fd < 0) continue;
            (void)::read(fd, buf, sizeof(buf) - 1);
            ::close(fd);
            int base = ::atoi(buf);

            // Prefer a chip whose label names the BCM pinctrl driver.
            char lbl[64]{};
            fd = ::open((dir + "/label").c_str(), O_RDONLY);
            if (fd >= 0) {
                (void)::read(fd, lbl, sizeof(lbl) - 1);
                ::close(fd);
                if (::strstr(lbl, "pinctrl") || ::strstr(lbl, "bcm")) {
                    ::closedir(d);
                    return base;
                }
            }

            if (base > best_base) best_base = base;
        }
        ::closedir(d);
        return best_base;
    }

    static void sysfs_write(const std::string &path, const std::string &val) {
        int fd = ::open(path.c_str(), O_WRONLY);
        if (fd < 0) return;  // silently ignore (e.g. already exported)
        (void)::write(fd, val.c_str(), val.size());
        ::close(fd);
    }
};
