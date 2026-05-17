#include "imu/mpu9250.h"

#include <chrono>
#include <cmath>
#include <thread>

#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>

// I2C_SLAVE is 0x0703 on all Linux architectures; define it here so we don't
// need linux/i2c-dev.h (not installed in Alpine musl build environments).
#ifndef I2C_SLAVE
#define I2C_SLAVE 0x0703
#endif

namespace imu {

namespace {

constexpr std::uint8_t AK_ADDR = 0x0C;  // AK8963 magnetometer (via MPU bypass)

// MPU-9250 registers (datasheet §3)
constexpr std::uint8_t REG_SMPLRT_DIV   = 0x19;
constexpr std::uint8_t REG_CONFIG       = 0x1A;
constexpr std::uint8_t REG_GYRO_CONFIG  = 0x1B;
constexpr std::uint8_t REG_ACCEL_CONFIG = 0x1C;
constexpr std::uint8_t REG_INT_PIN_CFG  = 0x37;
constexpr std::uint8_t REG_ACCEL_XOUT_H = 0x3B;
constexpr std::uint8_t REG_PWR_MGMT_1   = 0x6B;
constexpr std::uint8_t REG_WHO_AM_I     = 0x75;

// AK8963 registers
constexpr std::uint8_t AK_WIA   = 0x00;
constexpr std::uint8_t AK_ST1   = 0x02;
constexpr std::uint8_t AK_HXL   = 0x03;
constexpr std::uint8_t AK_CNTL1 = 0x0A;
constexpr std::uint8_t AK_ASAX  = 0x10;

inline std::int16_t s16(std::uint8_t hi, std::uint8_t lo) {
    return static_cast<std::int16_t>((std::uint16_t(hi) << 8) | lo);
}

}  // namespace

Mpu9250::Mpu9250(const std::string &bus, std::uint8_t mpu_addr)
    : mpu_addr_(mpu_addr) {
    fd_ = ::open(bus.c_str(), O_RDWR);
}

Mpu9250::~Mpu9250() {
    if (fd_ >= 0) ::close(fd_);
}

bool Mpu9250::set_slave(std::uint8_t addr) {
    return ::ioctl(fd_, I2C_SLAVE, addr) >= 0;
}

bool Mpu9250::write_reg(std::uint8_t addr, std::uint8_t reg, std::uint8_t val) {
    if (!set_slave(addr)) return false;
    std::uint8_t buf[2] = { reg, val };
    return ::write(fd_, buf, 2) == 2;
}

bool Mpu9250::read_bytes(std::uint8_t addr, std::uint8_t reg,
                         std::uint8_t *buf, int len) {
    if (!set_slave(addr)) return false;
    if (::write(fd_, &reg, 1) != 1) return false;
    return ::read(fd_, buf, len) == len;
}

bool Mpu9250::init() {
    if (fd_ < 0) return false;

    using namespace std::chrono_literals;

    std::uint8_t who = 0;
    if (!read_bytes(mpu_addr_, REG_WHO_AM_I, &who, 1)) return false;
    // MPU-9250 = 0x71, some clones (MPU-9255) = 0x73. Be tolerant.
    if (who != 0x71 && who != 0x73 && who != 0x70) return false;

    if (!write_reg(mpu_addr_, REG_PWR_MGMT_1, 0x80)) return false;  // reset
    std::this_thread::sleep_for(100ms);
    if (!write_reg(mpu_addr_, REG_PWR_MGMT_1, 0x01)) return false;  // PLL w/ X gyro
    std::this_thread::sleep_for(50ms);

    write_reg(mpu_addr_, REG_CONFIG,       0x03);  // DLPF: 41 Hz
    write_reg(mpu_addr_, REG_SMPLRT_DIV,   0x04);  // 1 kHz / (1+4) = 200 Hz
    write_reg(mpu_addr_, REG_GYRO_CONFIG,  0x00);  // ±250 dps  → 131 LSB/dps
    write_reg(mpu_addr_, REG_ACCEL_CONFIG, 0x00);  // ±2 g       → 16384 LSB/g
    write_reg(mpu_addr_, REG_INT_PIN_CFG,  0x02);  // BYPASS_EN → AK8963 visible

    std::this_thread::sleep_for(10ms);

    std::uint8_t ak_who = 0;
    if (!read_bytes(AK_ADDR, AK_WIA, &ak_who, 1)) return false;
    if (ak_who != 0x48) return false;

    write_reg(AK_ADDR, AK_CNTL1, 0x00);    // power down
    std::this_thread::sleep_for(10ms);
    write_reg(AK_ADDR, AK_CNTL1, 0x0F);    // fuse ROM access
    std::this_thread::sleep_for(10ms);

    std::uint8_t asa[3] = {0};
    if (read_bytes(AK_ADDR, AK_ASAX, asa, 3)) {
        for (int i = 0; i < 3; i++)
            mag_asa_[i] = (static_cast<float>(asa[i]) - 128.0f) / 256.0f + 1.0f;
    }

    write_reg(AK_ADDR, AK_CNTL1, 0x00);    // power down
    std::this_thread::sleep_for(10ms);
    write_reg(AK_ADDR, AK_CNTL1, 0x16);    // continuous 2 (100 Hz) + 16-bit
    std::this_thread::sleep_for(10ms);

    return true;
}

bool Mpu9250::read(Reading &out) {
    out.mag_valid = false;
    if (fd_ < 0) return false;

    std::uint8_t buf[14];
    if (!read_bytes(mpu_addr_, REG_ACCEL_XOUT_H, buf, 14)) return false;

    std::int16_t ax = s16(buf[0],  buf[1]);
    std::int16_t ay = s16(buf[2],  buf[3]);
    std::int16_t az = s16(buf[4],  buf[5]);
    // buf[6..7] is temperature.
    std::int16_t gx = s16(buf[8],  buf[9]);
    std::int16_t gy = s16(buf[10], buf[11]);
    std::int16_t gz = s16(buf[12], buf[13]);

    constexpr float G       = 9.80665f;
    constexpr float DEG2RAD = static_cast<float>(M_PI) / 180.0f;
    out.accel[0] = ax / accel_lsb_per_g_  * G;
    out.accel[1] = ay / accel_lsb_per_g_  * G;
    out.accel[2] = az / accel_lsb_per_g_  * G;
    out.gyro[0]  = gx / gyro_lsb_per_dps_ * DEG2RAD;
    out.gyro[1]  = gy / gyro_lsb_per_dps_ * DEG2RAD;
    out.gyro[2]  = gz / gyro_lsb_per_dps_ * DEG2RAD;

    // Magnetometer is asynchronous (continuous mode 2 = 100 Hz). We poll ST1
    // to find out whether a new sample is ready; if not, we just leave
    // mag_valid=false and the caller carries the last value.
    std::uint8_t st1 = 0;
    if (!read_bytes(AK_ADDR, AK_ST1, &st1, 1)) return true;  // accel/gyro ok
    if (!(st1 & 0x01)) return true;  // no fresh mag sample

    std::uint8_t mb[7];
    if (!read_bytes(AK_ADDR, AK_HXL, mb, 7)) return true;
    std::uint8_t st2 = mb[6];
    if (st2 & 0x08) return true;     // HOFL: magnetic overflow → drop

    // AK8963 outputs little-endian.
    std::int16_t mx = s16(mb[1], mb[0]);
    std::int16_t my = s16(mb[3], mb[2]);
    std::int16_t mz = s16(mb[5], mb[4]);

    // 16-bit mode: 0.15 µT/LSB. Apply per-axis factory sensitivity.
    constexpr float MAG_LSB = 0.15f;
    float mxc = mx * MAG_LSB * mag_asa_[0];
    float myc = my * MAG_LSB * mag_asa_[1];
    float mzc = mz * MAG_LSB * mag_asa_[2];

    // AK8963 axes are rotated relative to MPU-9250 body axes (datasheet
    // fig. 9-1): mpu_x = ak_y, mpu_y = ak_x, mpu_z = -ak_z.
    out.mag[0] =  myc;
    out.mag[1] =  mxc;
    out.mag[2] = -mzc;
    out.mag_valid = true;
    return true;
}

}  // namespace imu
