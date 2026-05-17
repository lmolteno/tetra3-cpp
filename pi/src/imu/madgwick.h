#pragma once

namespace imu {

// Madgwick AHRS filter — quaternion sensor fusion of accel + gyro + mag.
// Based on Sebastian Madgwick's 2010 paper / reference implementation.
//
// Frame convention (matches Madgwick's reference): q is the orientation of
// the sensor body relative to a North-East-Down world frame, applied as
//   v_world = R(q) * v_body
// where R(q) is the standard Hamilton quaternion→matrix.  q[0..3] = [w,x,y,z].
class MadgwickAhrs {
public:
    // beta = filter gain. Higher → more weight on accel/mag (faster correction
    // of gyro drift, more accel/mag noise leakage). 0.1 is a reasonable default
    // for a near-stationary star tracker.
    explicit MadgwickAhrs(float beta = 0.1f) : beta_(beta) {}

    // gyro:  rad/s, accel: any consistent unit (normalised internally),
    // mag:   any consistent unit (normalised internally), dt: seconds.
    // mag_valid=false → fall back to the 6-DOF (accel+gyro only) variant.
    void update(const float gyro[3], const float accel[3],
                const float mag[3], bool mag_valid, float dt);

    // TRIAD-style bootstrap from a single accel + mag sample. Useful at
    // start-up so the filter doesn't drift in from q=[1,0,0,0] over the first
    // few seconds. Safe to call repeatedly when stationary.
    void set_from_accel_mag(const float accel[3], const float mag[3]);

    void quaternion(float q_out[4]) const {
        q_out[0] = q_[0]; q_out[1] = q_[1];
        q_out[2] = q_[2]; q_out[3] = q_[3];
    }

private:
    float q_[4] = {1.0f, 0.0f, 0.0f, 0.0f};
    float beta_;
};

}  // namespace imu
