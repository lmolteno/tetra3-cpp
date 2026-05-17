#pragma once

#include <cstdint>
#include <string>

namespace imu {

struct Reading {
    float accel[3];   // m/s²,  body frame
    float gyro[3];    // rad/s, body frame
    float mag[3];     // µT,    body frame (already remapped from AK8963 axes)
    bool  mag_valid;  // true if a fresh mag sample arrived this read
};

// Bare-bones MPU-9250 I2C driver.
//
// Opens its own fd to the bus — the kernel serialises per-bus access, so this
// is safe to run on a thread while the OLED is being driven from the main
// loop on the same /dev/i2c-N.
//
// AK8963 magnetometer is reached via bypass mode (same bus, address 0x0C).
// Accel/gyro/mag samples are remapped into a single right-handed body frame
// matching the MPU-9250 chip frame (AK8963's X/Y are swapped and its Z is
// inverted relative to the MPU, per datasheet figure 9-1).
class Mpu9250 {
public:
    // bus       e.g. "/dev/i2c-1"
    // mpu_addr  0x68 (AD0 low) or 0x69 (AD0 high)
    Mpu9250(const std::string &bus, std::uint8_t mpu_addr = 0x68);
    ~Mpu9250();

    Mpu9250(const Mpu9250 &)            = delete;
    Mpu9250 &operator=(const Mpu9250 &) = delete;

    // Configure MPU + AK8963 (DLPF, scales, bypass enable, continuous mag).
    // Returns false on bus open failure, WHO_AM_I mismatch, or transfer
    // failure. After a failure the object is non-functional.
    bool init();

    // Pull one fresh accel+gyro sample. The magnetometer runs at a slower
    // rate (100 Hz continuous mode 2) so it's only valid sometimes —
    // out.mag_valid reports that.
    bool read(Reading &out);

    bool ok() const { return fd_ >= 0; }

private:
    int           fd_              = -1;
    std::uint8_t  mpu_addr_        = 0x68;
    // Factory per-axis sensitivity correction read out of the AK8963 fuse ROM.
    float         mag_asa_[3]      = {1.0f, 1.0f, 1.0f};
    // Sensor scales chosen in init(); kept as members so read() doesn't hard-
    // code them.
    float         accel_lsb_per_g_  = 16384.0f;  // ±2g
    float         gyro_lsb_per_dps_ = 131.0f;    // ±250 dps

    bool set_slave(std::uint8_t addr);
    bool write_reg(std::uint8_t addr, std::uint8_t reg, std::uint8_t val);
    bool read_bytes(std::uint8_t addr, std::uint8_t reg,
                    std::uint8_t *buf, int len);
};

}  // namespace imu
