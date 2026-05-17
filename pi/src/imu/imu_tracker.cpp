#include "imu/imu_tracker.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

#include "imu/mpu9250.h"
#include "imu/sky_math.h"

namespace imu {

namespace {

double unix_now() {
    using namespace std::chrono;
    return duration<double>(
        system_clock::now().time_since_epoch()).count();
}

// ── Small quaternion / matrix helpers (3x3 row-major) ────────────────────────

void qmul(const float a[4], const float b[4], float out[4]) {
    out[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    out[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    out[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    out[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

void qinv(const float a[4], float out[4]) {
    out[0] =  a[0]; out[1] = -a[1]; out[2] = -a[2]; out[3] = -a[3];
}

// Build (axis, angle) from a rotation matrix. Returns angle in radians;
// axis has unit length. For very small angles axis is set to (1,0,0).
void R_to_axis_angle(const float R[9], float axis[3], float &angle) {
    float tr = R[0] + R[4] + R[8];
    float c = (tr - 1.0f) * 0.5f;
    if (c > 1.0f) c = 1.0f;
    if (c < -1.0f) c = -1.0f;
    angle = std::acos(c);
    float s = std::sin(angle);
    if (std::fabs(s) < 1e-5f) {
        axis[0] = 1; axis[1] = 0; axis[2] = 0;
        return;
    }
    axis[0] = (R[7] - R[5]) / (2 * s);
    axis[1] = (R[2] - R[6]) / (2 * s);
    axis[2] = (R[3] - R[1]) / (2 * s);
    float n = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (n > 0) { axis[0]/=n; axis[1]/=n; axis[2]/=n; }
}

// Kabsch fit: given paired unit vectors {a_i}, {b_i}, find R such that
// R * b_i ≈ a_i. Returns R as 3x3 row-major (proper rotation, det = +1).
// Falls back to identity if rank-deficient.
bool kabsch(const std::vector<std::array<float,3>> &A,
            const std::vector<std::array<float,3>> &B,
            float R[9]) {
    if (A.size() < 2 || A.size() != B.size()) return false;
    // H = sum b_i * a_i^T (3x3).
    double H[9] = {0};
    for (size_t i = 0; i < A.size(); i++) {
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                H[r*3 + c] += B[i][r] * A[i][c];
    }

    // SVD of H via Jacobi eigen-decomposition of H^T H. Tiny matrix, fine.
    double HtH[9];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            double s = 0;
            for (int k = 0; k < 3; k++) s += H[k*3 + r] * H[k*3 + c];
            HtH[r*3 + c] = s;
        }

    double V[9] = {1,0,0, 0,1,0, 0,0,1};
    double A_[9]; for (int i = 0; i < 9; i++) A_[i] = HtH[i];
    for (int iter = 0; iter < 30; iter++) {
        // Pick largest off-diagonal.
        int p = 0, q = 1; double mx = std::fabs(A_[1]);
        for (int r = 0; r < 3; r++) for (int c = r+1; c < 3; c++) {
            double v = std::fabs(A_[r*3 + c]);
            if (v > mx) { mx = v; p = r; q = c; }
        }
        if (mx < 1e-12) break;
        double app = A_[p*3 + p], aqq = A_[q*3 + q], apq = A_[p*3 + q];
        double phi = (aqq - app) / (2 * apq);
        double t = phi >= 0 ? 1.0/(phi + std::sqrt(1+phi*phi))
                            : 1.0/(phi - std::sqrt(1+phi*phi));
        double c = 1.0/std::sqrt(1+t*t), s = t*c;
        for (int i = 0; i < 3; i++) {
            double aip = A_[i*3 + p], aiq = A_[i*3 + q];
            A_[i*3 + p] = c*aip - s*aiq;
            A_[i*3 + q] = s*aip + c*aiq;
        }
        for (int j = 0; j < 3; j++) {
            double apj = A_[p*3 + j], aqj = A_[q*3 + j];
            A_[p*3 + j] = c*apj - s*aqj;
            A_[q*3 + j] = s*apj + c*aqj;
            double vip = V[j*3 + p], viq = V[j*3 + q];
            V[j*3 + p] = c*vip - s*viq;
            V[j*3 + q] = s*vip + c*viq;
        }
    }
    // Eigenvalues on diagonal of A_, eigenvectors in columns of V.
    double sv[3] = {std::sqrt(std::max(0.0, A_[0])),
                    std::sqrt(std::max(0.0, A_[4])),
                    std::sqrt(std::max(0.0, A_[8]))};
    // U = H * V * diag(1/sv).
    double U[9];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            double s = 0;
            for (int k = 0; k < 3; k++) s += H[r*3 + k] * V[k*3 + c];
            U[r*3 + c] = sv[c] > 1e-9 ? s / sv[c] : 0.0;
        }

    // R = V * U^T (with reflection correction).
    double VUt[9];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            double s = 0;
            for (int k = 0; k < 3; k++) s += V[r*3 + k] * U[c*3 + k];
            VUt[r*3 + c] = s;
        }
    double det =
          VUt[0]*(VUt[4]*VUt[8] - VUt[5]*VUt[7])
        - VUt[1]*(VUt[3]*VUt[8] - VUt[5]*VUt[6])
        + VUt[2]*(VUt[3]*VUt[7] - VUt[4]*VUt[6]);
    if (det < 0) {
        // Flip the smallest singular value's sign.
        int kmin = 0; for (int k = 1; k < 3; k++) if (sv[k] < sv[kmin]) kmin = k;
        for (int r = 0; r < 3; r++) U[r*3 + kmin] = -U[r*3 + kmin];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) {
                double s = 0;
                for (int k = 0; k < 3; k++) s += V[r*3 + k] * U[c*3 + k];
                VUt[r*3 + c] = s;
            }
    }
    for (int i = 0; i < 9; i++) R[i] = static_cast<float>(VUt[i]);
    return true;
}

}  // namespace

ImuTracker::ImuTracker(const Config &cfg)
    : cfg_(cfg), filter_(cfg.beta) {}

ImuTracker::~ImuTracker() { stop(); }

bool ImuTracker::calibrated() const {
    std::lock_guard<std::mutex> lk(mutex_);
    return calibrated_;
}

bool ImuTracker::start() {
    imu_ = std::make_unique<Mpu9250>(cfg_.i2c_bus, cfg_.mpu_addr);
    if (!imu_->ok() || !imu_->init()) {
        imu_.reset();
        return false;
    }
    if (!cfg_.imu_log_path.empty()) {
        log_.open(cfg_.imu_log_path, std::ios::out | std::ios::trunc);
        if (log_.is_open()) {
            log_ << "# pi_tracker IMU sample log\n"
                 << "# I,t_unix,ax,ay,az,gx,gy,gz,mx,my,mz,mag_valid,qw,qx,qy,qz\n"
                 << "#   raw a (m/s² body), raw g (rad/s body, NOT bias-subtracted),\n"
                 << "#   raw m (µT body), q is body→world from 6-DOF Madgwick\n"
                 << "# S,t_unix,ra_rad,dec_rad,roll_rad\n"
                 << "# B,t_unix,bx,by,bz   (gyro bias, rad/s body, subtracted before filter)\n"
                 << "# K,t_unix,R00..R22,calib_total_deg  (hand-eye fit on calibration)\n";
            log_.flush();
        } else {
            std::cerr << "IMU log: failed to open " << cfg_.imu_log_path << "\n";
        }
    }
    running_.store(true);
    thread_ = std::thread([this] { run(); });
    return true;
}

void ImuTracker::stop() {
    if (running_.exchange(false) && thread_.joinable()) thread_.join();
    imu_.reset();
}

void ImuTracker::run() {
    using clock = std::chrono::steady_clock;
    const auto period = std::chrono::nanoseconds(
        static_cast<std::int64_t>(1e9 / cfg_.rate_hz));
    auto next = clock::now();

    // ── 2 s of stationary gyro bias calibration ────────────────────────────
    float gyro_bias[3] = {0, 0, 0};
    {
        constexpr int N_BIAS = 200;
        int got = 0;
        Reading r;
        while (got < N_BIAS && running_.load()) {
            next += period;
            if (imu_->read(r)) {
                gyro_bias[0] += r.gyro[0];
                gyro_bias[1] += r.gyro[1];
                gyro_bias[2] += r.gyro[2];
                got++;
            }
            std::this_thread::sleep_until(next);
        }
        if (got > 0) {
            gyro_bias[0] /= got;
            gyro_bias[1] /= got;
            gyro_bias[2] /= got;
        }
        constexpr float RAD2DEG = 180.0f / 3.14159265358979f;
        std::cerr << "IMU gyro bias (deg/s): "
                  << gyro_bias[0] * RAD2DEG << " "
                  << gyro_bias[1] * RAD2DEG << " "
                  << gyro_bias[2] * RAD2DEG
                  << " (avg of " << got << " samples)\n";
        if (log_.is_open()) {
            std::lock_guard<std::mutex> lk(log_mutex_);
            log_ << "B," << std::fixed << unix_now()
                 << "," << gyro_bias[0]
                 << "," << gyro_bias[1]
                 << "," << gyro_bias[2] << "\n";
            log_.flush();
        }
    }

    auto last_t = clock::now();
    while (running_.load()) {
        next += period;

        Reading r;
        if (imu_->read(r)) {
            auto now = clock::now();
            float dt = std::chrono::duration<float>(now - last_t).count();
            last_t = now;
            if (dt <= 0 || dt > 0.5f) dt = 1.0f / cfg_.rate_hz;

            float gyro[3] = {
                r.gyro[0] - gyro_bias[0],
                r.gyro[1] - gyro_bias[1],
                r.gyro[2] - gyro_bias[2],
            };

            // 6-DOF Madgwick: accel keeps gravity locked, gyro integrates the
            // rest. We never feed mag — its absolute heading is unreliable
            // and we don't need an absolute world frame; hand-eye does the
            // alignment via plate solves.
            float dummy_mag[3] = {0, 0, 0};
            if (!bootstrapped_) {
                filter_.set_from_accel_mag(r.accel, dummy_mag);
                std::lock_guard<std::mutex> lk(mutex_);
                filter_.quaternion(q_);
                bootstrapped_ = true;
            } else {
                filter_.update(gyro, r.accel, dummy_mag, /*mag_valid=*/false, dt);
                std::lock_guard<std::mutex> lk(mutex_);
                filter_.quaternion(q_);
            }

            if (log_.is_open()) {
                float qsnap[4];
                filter_.quaternion(qsnap);
                std::lock_guard<std::mutex> lk(log_mutex_);
                log_ << "I," << std::fixed << unix_now()
                     << "," << r.accel[0] << "," << r.accel[1] << "," << r.accel[2]
                     << "," << r.gyro[0]  << "," << r.gyro[1]  << "," << r.gyro[2]
                     << "," << r.mag[0]   << "," << r.mag[1]   << "," << r.mag[2]
                     << "," << (r.mag_valid ? 1 : 0)
                     << "," << qsnap[0] << "," << qsnap[1] << "," << qsnap[2] << "," << qsnap[3]
                     << "\n";
                if (++log_samples_since_flush_ >= 100) {
                    log_.flush();
                    log_samples_since_flush_ = 0;
                }
            }
        }

        std::this_thread::sleep_until(next);
    }
}

std::optional<Pointing> ImuTracker::get_pointing() const {
    std::lock_guard<std::mutex> lk(mutex_);
    if (!calibrated_) return std::nullopt;

    // body delta in body frame: q_delta_body = anchor_q^{-1} * q_current
    float q_anchor_inv[4]; qinv(anchor_q_, q_anchor_inv);
    float q_delta[4]; qmul(q_anchor_inv, q_, q_delta);
    // Normalise (defensive)
    float qn = 1.0f / std::sqrt(q_delta[0]*q_delta[0] + q_delta[1]*q_delta[1]
                              + q_delta[2]*q_delta[2] + q_delta[3]*q_delta[3]);
    q_delta[0]*=qn; q_delta[1]*=qn; q_delta[2]*=qn; q_delta[3]*=qn;

    float R_body_delta[9];
    sky::quat_to_matrix(q_delta, R_body_delta);

    // Camera-frame delta: R_cam_delta = R_cb^T * R_body_delta * R_cb
    float R_cb_T[9]; sky::transpose3(R_cb_, R_cb_T);
    float tmp[9]; sky::matmul3(R_cb_T, R_body_delta, tmp);
    float R_cam_delta[9]; sky::matmul3(tmp, R_cb_, R_cam_delta);

    // R_cam_to_cel(t) = anchor_R * R_cam_delta
    float R_cam_to_cel[9]; sky::matmul3(anchor_R_, R_cam_delta, R_cam_to_cel);

    Pointing p;
    sky::matrix_to_radec_roll(R_cam_to_cel, p.ra, p.dec, p.roll);
    p.calibrated = true;
    return p;
}

void ImuTracker::update_solve(float ra, float dec, float roll) {
    double t = unix_now();
    if (log_.is_open()) {
        std::lock_guard<std::mutex> lk(log_mutex_);
        log_ << "S," << std::fixed << t
             << "," << ra << "," << dec << "," << roll << "\n";
        log_.flush();
    }

    std::lock_guard<std::mutex> lk(mutex_);
    if (!bootstrapped_) return;

    float R_solve[9];
    sky::radec_roll_to_matrix(ra, dec, roll, R_solve);

    if (calibrated_) {
        // Just refresh the anchor.
        anchor_q_[0]=q_[0]; anchor_q_[1]=q_[1]; anchor_q_[2]=q_[2]; anchor_q_[3]=q_[3];
        for (int i = 0; i < 9; i++) anchor_R_[i] = R_solve[i];
        return;
    }

    // Pre-calibration: append to the buffer.
    SolvePair sp;
    sp.q_body[0]=q_[0]; sp.q_body[1]=q_[1]; sp.q_body[2]=q_[2]; sp.q_body[3]=q_[3];
    for (int i = 0; i < 9; i++) sp.R_solve[i] = R_solve[i];

    if (calib_buffer_.empty()) calib_first_solve_t_ = t;
    calib_buffer_.push_back(sp);

    // Try to fit if enough time has elapsed and we have enough pairs.
    double elapsed = t - calib_first_solve_t_;
    if (elapsed < cfg_.calib_window_s ||
        (int)calib_buffer_.size() < cfg_.calib_min_pairs + 1) return;

    // Build (body_axis_i, camera_axis_i) pairs from each solve vs. solve_0.
    const SolvePair &p0 = calib_buffer_.front();
    float q0_inv[4]; qinv(p0.q_body, q0_inv);
    float R0_solve_T[9]; sky::transpose3(p0.R_solve, R0_solve_T);

    std::vector<std::array<float,3>> body_axes, cam_axes;
    float total_cam_deg = 0.0f;
    constexpr float DEG = 180.0f / 3.14159265358979f;
    for (size_t i = 1; i < calib_buffer_.size(); i++) {
        const SolvePair &pi = calib_buffer_[i];
        float qd[4]; qmul(q0_inv, pi.q_body, qd);
        float R_body[9]; sky::quat_to_matrix(qd, R_body);
        float R_cam[9]; sky::matmul3(R0_solve_T, pi.R_solve, R_cam);
        float ax_b[3], ax_c[3], ang_b, ang_c;
        R_to_axis_angle(R_body, ax_b, ang_b);
        R_to_axis_angle(R_cam,  ax_c, ang_c);
        // Skip pairs with tiny camera-frame rotation — they're noise.
        if (ang_c < 5.0f / DEG) continue;
        body_axes.push_back({ax_b[0], ax_b[1], ax_b[2]});
        cam_axes.push_back({ax_c[0], ax_c[1], ax_c[2]});
        total_cam_deg += ang_c * DEG;
    }

    if ((int)body_axes.size() < cfg_.calib_min_pairs ||
        total_cam_deg < cfg_.calib_min_total_deg) {
        // Not yet enough motion; keep buffering.
        return;
    }

    // Solve R_cb * cam_axis ≈ body_axis.
    float R_cb_fit[9];
    if (!kabsch(body_axes, cam_axes, R_cb_fit)) return;

    for (int i = 0; i < 9; i++) R_cb_[i] = R_cb_fit[i];
    // Anchor to the most recent solve in the buffer.
    const SolvePair &latest = calib_buffer_.back();
    for (int i = 0; i < 4; i++) anchor_q_[i] = latest.q_body[i];
    for (int i = 0; i < 9; i++) anchor_R_[i] = latest.R_solve[i];
    calibrated_ = true;
    calib_buffer_.clear();

    float ax[3], ang;
    R_to_axis_angle(R_cb_, ax, ang);
    std::cerr << "IMU calibrated after " << elapsed << "s ("
              << body_axes.size() << " axis pairs, "
              << total_cam_deg << "° total). R_cam->body ≡ "
              << (ang * DEG) << "° about ["
              << ax[0] << "," << ax[1] << "," << ax[2] << "]\n";

    if (log_.is_open()) {
        std::lock_guard<std::mutex> lkl(log_mutex_);
        log_ << "K," << std::fixed << t;
        for (int i = 0; i < 9; i++) log_ << "," << R_cb_[i];
        log_ << "," << total_cam_deg << "\n";
        log_.flush();
    }
}

}  // namespace imu
