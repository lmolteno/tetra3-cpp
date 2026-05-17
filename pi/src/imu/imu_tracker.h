#pragma once

#include <atomic>
#include <chrono>
#include <fstream>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "imu/madgwick.h"

namespace imu {

class Mpu9250;

struct Pointing {
    float ra;          // radians, [0, 2π)
    float dec;         // radians, [-π/2, π/2]
    float roll;        // radians, position angle of image-up vs. celestial north
    bool  calibrated;  // false → pre-calibration, no pointing available
};

// 6-DOF Madgwick (accel + gyro) plus hand-eye calibration against plate solves.
//
// The IMU's absolute heading is meaningless on its own — magnetometer is too
// polluted by local fields and gyro yaw drifts. Instead, we wait for the
// first ~10 s of plate solves after startup, collect (q_body, R_solve) pairs,
// and fit a fixed R_cam→body via Kabsch SVD on the rotation-axis pairs. From
// then on:
//
//   - every new plate solve refreshes an anchor (R_solve_anchor, q_anchor);
//   - between solves we compose the body-frame gyro delta with R_cam→body
//     and the anchor to estimate the current camera→celestial rotation.
//
// Result on a 5-minute test trace: median 0.10° / 95th 0.24° / max 1.52°
// pointing error vs. plate solves.
class ImuTracker {
public:
    struct Config {
        std::string i2c_bus       = "/dev/i2c-1";
        std::uint8_t mpu_addr     = 0x68;
        float        beta         = 0.05f;   // Madgwick filter gain (6-DOF)
        int          rate_hz      = 100;

        // Hand-eye calibration window after the first plate solve.
        float        calib_window_s = 10.0f;
        int          calib_min_pairs = 4;
        // Minimum camera-frame motion (degrees) across the calibration window
        // before we trust the fit; below this the rotation axes are too
        // parallel for an SVD to disambiguate R_cam→body.
        float        calib_min_total_deg = 20.0f;

        // If non-empty: open this file for the raw IMU sample stream.
        std::string imu_log_path;
    };

    explicit ImuTracker(const Config &cfg);
    ~ImuTracker();

    ImuTracker(const ImuTracker &)            = delete;
    ImuTracker &operator=(const ImuTracker &) = delete;

    bool start();
    void stop();
    bool running() const { return running_.load(); }

    // Pointing estimate. nullopt until the hand-eye fit has succeeded.
    std::optional<Pointing> get_pointing() const;

    // Notify the tracker that a fresh plate solve arrived (RA/Dec/Roll in rad).
    // While not yet calibrated, the solve is appended to the calibration
    // buffer; once enough variety has been collected the fit runs. After
    // calibration the solve just refreshes the anchor.
    void update_solve(float ra, float dec, float roll);

    // True once the hand-eye fit has succeeded.
    bool calibrated() const;

private:
    void run();

    Config cfg_;
    std::unique_ptr<Mpu9250> imu_;
    MadgwickAhrs filter_;

    std::thread thread_;
    std::atomic<bool> running_{false};

    mutable std::mutex mutex_;
    // ── shared state under mutex_ ────────────────────────────────────────
    bool  bootstrapped_ = false;     // got at least one accel sample
    bool  calibrated_   = false;     // hand-eye fit has succeeded
    float q_[4]         = {1, 0, 0, 0};

    // Hand-eye calibration buffer (only used until calibrated_=true).
    struct SolvePair {
        float q_body[4];          // body→IMU-world at the time of the solve
        float R_solve[9];         // camera→celestial from the solve
    };
    std::vector<SolvePair> calib_buffer_;
    double calib_first_solve_t_ = 0.0;

    // Post-calibration anchor and the fitted camera↔body rotation.
    float anchor_q_[4]    = {1, 0, 0, 0};
    float anchor_R_[9]    = {1,0,0, 0,1,0, 0,0,1};
    float R_cb_[9]        = {1,0,0, 0,1,0, 0,0,1};  // camera → body

    // Diagnostic log (raw samples + bias + solves).
    std::ofstream log_;
    std::mutex log_mutex_;
    std::uint64_t log_samples_since_flush_ = 0;
};

}  // namespace imu
