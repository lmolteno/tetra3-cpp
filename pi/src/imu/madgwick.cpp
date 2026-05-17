#include "imu/madgwick.h"

#include <cmath>

namespace imu {

namespace {

inline float inv_sqrt(float x) { return 1.0f / std::sqrt(x); }

}  // namespace

void MadgwickAhrs::update(const float gyro[3], const float accel[3],
                          const float mag[3], bool mag_valid, float dt) {
    float q0 = q_[0], q1 = q_[1], q2 = q_[2], q3 = q_[3];
    float gx = gyro[0], gy = gyro[1], gz = gyro[2];

    // Quaternion derivative from gyro alone.
    float qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
    float qDot2 = 0.5f * ( q0 * gx + q2 * gz - q3 * gy);
    float qDot3 = 0.5f * ( q0 * gy - q1 * gz + q3 * gx);
    float qDot4 = 0.5f * ( q0 * gz + q1 * gy - q2 * gx);

    float ax = accel[0], ay = accel[1], az = accel[2];
    float anorm2 = ax * ax + ay * ay + az * az;
    if (anorm2 > 1e-6f) {
        float ainv = inv_sqrt(anorm2);
        ax *= ainv; ay *= ainv; az *= ainv;

        if (mag_valid) {
            float mx = mag[0], my = mag[1], mz = mag[2];
            float mnorm2 = mx * mx + my * my + mz * mz;
            if (mnorm2 > 1e-6f) {
                float minv = inv_sqrt(mnorm2);
                mx *= minv; my *= minv; mz *= minv;

                // ── Madgwick MARG gradient (9-DOF) ────────────────────────
                float _2q0mx = 2.0f * q0 * mx;
                float _2q0my = 2.0f * q0 * my;
                float _2q0mz = 2.0f * q0 * mz;
                float _2q1mx = 2.0f * q1 * mx;
                float _2q0 = 2.0f * q0;
                float _2q1 = 2.0f * q1;
                float _2q2 = 2.0f * q2;
                float _2q3 = 2.0f * q3;
                float _2q0q2 = 2.0f * q0 * q2;
                float _2q2q3 = 2.0f * q2 * q3;
                float q0q0 = q0 * q0, q0q1 = q0 * q1, q0q2 = q0 * q2, q0q3 = q0 * q3;
                float q1q1 = q1 * q1, q1q2 = q1 * q2, q1q3 = q1 * q3;
                float q2q2 = q2 * q2, q2q3 = q2 * q3;
                float q3q3 = q3 * q3;

                // Reference direction of Earth's magnetic field after rotation
                // into world frame, then projected to (bx, 0, bz) by removing
                // the east component.
                float hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 +
                           _2q1 * my * q2 + _2q1 * mz * q3 -
                           mx * q2q2 - mx * q3q3;
                float hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 +
                           _2q1mx * q2 - my * q1q1 + my * q2q2 +
                           _2q2 * mz * q3 - my * q3q3;
                float _2bx = std::sqrt(hx * hx + hy * hy);
                float _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 +
                             _2q1mx * q3 - mz * q1q1 +
                             _2q2 * my * q3 - mz * q2q2 + mz * q3q3;
                float _4bx = 2.0f * _2bx;
                float _4bz = 2.0f * _2bz;

                float s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax)
                         +  _2q1 * (2.0f * q0q1 + _2q2q3 - ay)
                         -  _2bz * q2 *
                             (_2bx * (0.5f - q2q2 - q3q3) +
                              _2bz * (q1q3 - q0q2) - mx)
                         + (-_2bx * q3 + _2bz * q1) *
                             (_2bx * (q1q2 - q0q3) +
                              _2bz * (q0q1 + q2q3) - my)
                         +  _2bx * q2 *
                             (_2bx * (q0q2 + q1q3) +
                              _2bz * (0.5f - q1q1 - q2q2) - mz);
                float s1 =  _2q3 * (2.0f * q1q3 - _2q0q2 - ax)
                         +  _2q0 * (2.0f * q0q1 + _2q2q3 - ay)
                         -  4.0f * q1 * (1.0f - 2.0f * q1q1 - 2.0f * q2q2 - az)
                         +  _2bz * q3 *
                             (_2bx * (0.5f - q2q2 - q3q3) +
                              _2bz * (q1q3 - q0q2) - mx)
                         + (_2bx * q2 + _2bz * q0) *
                             (_2bx * (q1q2 - q0q3) +
                              _2bz * (q0q1 + q2q3) - my)
                         + (_2bx * q3 - _4bz * q1) *
                             (_2bx * (q0q2 + q1q3) +
                              _2bz * (0.5f - q1q1 - q2q2) - mz);
                float s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax)
                         +  _2q3 * (2.0f * q0q1 + _2q2q3 - ay)
                         -  4.0f * q2 * (1.0f - 2.0f * q1q1 - 2.0f * q2q2 - az)
                         + (-_4bx * q2 - _2bz * q0) *
                             (_2bx * (0.5f - q2q2 - q3q3) +
                              _2bz * (q1q3 - q0q2) - mx)
                         + (_2bx * q1 + _2bz * q3) *
                             (_2bx * (q1q2 - q0q3) +
                              _2bz * (q0q1 + q2q3) - my)
                         + (_2bx * q0 - _4bz * q2) *
                             (_2bx * (q0q2 + q1q3) +
                              _2bz * (0.5f - q1q1 - q2q2) - mz);
                float s3 =  _2q1 * (2.0f * q1q3 - _2q0q2 - ax)
                         +  _2q2 * (2.0f * q0q1 + _2q2q3 - ay)
                         + (-_4bx * q3 + _2bz * q1) *
                             (_2bx * (0.5f - q2q2 - q3q3) +
                              _2bz * (q1q3 - q0q2) - mx)
                         + (-_2bx * q0 + _2bz * q2) *
                             (_2bx * (q1q2 - q0q3) +
                              _2bz * (q0q1 + q2q3) - my)
                         +  _2bx * q1 *
                             (_2bx * (q0q2 + q1q3) +
                              _2bz * (0.5f - q1q1 - q2q2) - mz);

                float sinv = inv_sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
                qDot1 -= beta_ * s0 * sinv;
                qDot2 -= beta_ * s1 * sinv;
                qDot3 -= beta_ * s2 * sinv;
                qDot4 -= beta_ * s3 * sinv;
            }
        } else {
            // ── 6-DOF (accel only) ─────────────────────────────────────────
            float _2q0 = 2.0f * q0;
            float _2q1 = 2.0f * q1;
            float _2q2 = 2.0f * q2;
            float _2q3 = 2.0f * q3;
            float _4q0 = 4.0f * q0;
            float _4q1 = 4.0f * q1;
            float _4q2 = 4.0f * q2;
            float _8q1 = 8.0f * q1;
            float _8q2 = 8.0f * q2;
            float q0q0 = q0 * q0;
            float q1q1 = q1 * q1;
            float q2q2 = q2 * q2;
            float q3q3 = q3 * q3;

            float s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
            float s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1
                     - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2
                     + _4q1 * az;
            float s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3
                     - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2
                     + _4q2 * az;
            float s3 = 4.0f * q1q1 * q3 - _2q1 * ax
                     + 4.0f * q2q2 * q3 - _2q2 * ay;
            float sinv = inv_sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
            qDot1 -= beta_ * s0 * sinv;
            qDot2 -= beta_ * s1 * sinv;
            qDot3 -= beta_ * s2 * sinv;
            qDot4 -= beta_ * s3 * sinv;
        }
    }

    q0 += qDot1 * dt;
    q1 += qDot2 * dt;
    q2 += qDot3 * dt;
    q3 += qDot4 * dt;

    float qinv = inv_sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q_[0] = q0 * qinv;
    q_[1] = q1 * qinv;
    q_[2] = q2 * qinv;
    q_[3] = q3 * qinv;
}

void MadgwickAhrs::set_from_accel_mag(const float accel[3], const float mag[3]) {
    // Build a body→world (NED) rotation from accel ("up" in body) and mag
    // (horizontal "north" after tilt-compensation), then convert to a
    // quaternion. Madgwick's convention: v_world = R(q) * v_body.

    float ax = accel[0], ay = accel[1], az = accel[2];
    float an = std::sqrt(ax*ax + ay*ay + az*az);
    if (an < 1e-6f) return;
    // Accel reads the support force opposing gravity, so its direction is the
    // body-frame "up". World "up" = -down = [0, 0, -1] in NED.
    float ux = ax/an, uy = ay/an, uz = az/an;

    float mx = mag[0], my = mag[1], mz = mag[2];
    float mn = std::sqrt(mx*mx + my*my + mz*mz);
    if (mn < 1e-6f) return;
    mx /= mn; my /= mn; mz /= mn;

    // Project mag onto the horizontal plane (orthogonal to up).
    float d = mx*ux + my*uy + mz*uz;
    float nx = mx - d*ux, ny = my - d*uy, nz = mz - d*uz;
    float nn = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (nn < 1e-6f) return;
    nx /= nn; ny /= nn; nz /= nn;        // "north" expressed in body coords

    // Down in body coords:
    float dx = -ux, dy = -uy, dz = -uz;
    // East = down × north (NED right-handed)
    float ex = dy*nz - dz*ny;
    float ey = dz*nx - dx*nz;
    float ez = dx*ny - dy*nx;

    // R(q) sends v_body to v_world. So columns of R = world basis vectors
    // expressed in… wait, no — R maps body to world, so columns of R are the
    // images of body basis vectors in the world frame. Equivalently, rows of R
    // are world basis vectors expressed in body coords: row 0 = north_body,
    // row 1 = east_body, row 2 = down_body.
    float R00 = nx, R01 = ny, R02 = nz;
    float R10 = ex, R11 = ey, R12 = ez;
    float R20 = dx, R21 = dy, R22 = dz;

    // Shepperd / Markley quaternion extraction.
    float trace = R00 + R11 + R22;
    float qw, qx, qy, qz;
    if (trace > 0.0f) {
        float s = 0.5f / std::sqrt(trace + 1.0f);
        qw = 0.25f / s;
        qx = (R21 - R12) * s;
        qy = (R02 - R20) * s;
        qz = (R10 - R01) * s;
    } else if (R00 > R11 && R00 > R22) {
        float s = 2.0f * std::sqrt(1.0f + R00 - R11 - R22);
        qw = (R21 - R12) / s;
        qx = 0.25f * s;
        qy = (R01 + R10) / s;
        qz = (R02 + R20) / s;
    } else if (R11 > R22) {
        float s = 2.0f * std::sqrt(1.0f + R11 - R00 - R22);
        qw = (R02 - R20) / s;
        qx = (R01 + R10) / s;
        qy = 0.25f * s;
        qz = (R12 + R21) / s;
    } else {
        float s = 2.0f * std::sqrt(1.0f + R22 - R00 - R11);
        qw = (R10 - R01) / s;
        qx = (R02 + R20) / s;
        qy = (R12 + R21) / s;
        qz = 0.25f * s;
    }
    float qn = 1.0f / std::sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
    q_[0] = qw * qn;
    q_[1] = qx * qn;
    q_[2] = qy * qn;
    q_[3] = qz * qn;
}

}  // namespace imu
