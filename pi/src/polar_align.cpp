#include "polar_align.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

PolarAligner::PolarAligner(float observer_lat_deg, float observer_lon_deg)
    : observer_lat_rad_(observer_lat_deg * static_cast<float>(M_PI) / 180.0f),
      observer_lon_rad_(observer_lon_deg * static_cast<float>(M_PI) / 180.0f) {}

void PolarAligner::clear_samples() { samples_.clear(); }

void PolarAligner::add_sample(float ra_rad, float dec_rad, double timestamp_s) {
    samples_.push_back({ra_rad, dec_rad, timestamp_s});
}

std::vector<PASample> PolarAligner::derotate_samples() const {
    if (samples_.empty()) return {};

    std::vector<PASample> out = samples_;
    double t0 = out[0].timestamp;

    for (auto &s : out) {
        double dt = s.timestamp - t0;
        // Earth rotation only shifts RA; dec is unchanged
        s.ra -= static_cast<float>(SIDEREAL_RATE * dt);
    }
    return out;
}

PoleEstimate PolarAligner::fit_circle(const std::vector<PASample> &pts) const {
    int n = static_cast<int>(pts.size());
    if (n < 3) return {};

    // Convert all points to unit vectors
    std::vector<double> vx(n), vy(n), vz(n);
    for (int i = 0; i < n; i++) {
        double cd = std::cos(pts[i].dec);
        vx[i] = cd * std::cos(pts[i].ra);
        vy[i] = cd * std::sin(pts[i].ra);
        vz[i] = std::sin(pts[i].dec);
    }

    // Spherical circle fit: all points should satisfy n·p_i = c (constant)
    // where n is the pole (unit) vector and c = cos(angular_radius).
    //
    // Rewrite as: n·p_i = n·p_j for all i,j, i.e. n·(p_i - p_j) = 0.
    // This means n is perpendicular to all difference vectors.
    // With n points we get n-1 constraints. Solve via SVD-like approach:
    //
    // Build matrix D where each row is (p_i - p_0), find n as the
    // direction that minimizes |D·n|^2 subject to |n|=1.
    // This is the left null-space, found via the smallest eigenvector of D^T·D.

    // Build D^T·D (3x3 symmetric matrix)
    double ATA[3][3] = {};
    for (int i = 1; i < n; i++) {
        double dx = vx[i] - vx[0];
        double dy = vy[i] - vy[0];
        double dz = vz[i] - vz[0];
        ATA[0][0] += dx * dx; ATA[0][1] += dx * dy; ATA[0][2] += dx * dz;
        ATA[1][1] += dy * dy; ATA[1][2] += dy * dz;
        ATA[2][2] += dz * dz;
    }
    ATA[1][0] = ATA[0][1];
    ATA[2][0] = ATA[0][2];
    ATA[2][1] = ATA[1][2];

    // Find the eigenvector with the smallest eigenvalue of ATA (3x3).
    // Use inverse power iteration starting from the mean direction
    // (which is a good initial guess for the pole).
    double mx = 0, my = 0, mz = 0;
    for (int i = 0; i < n; i++) { mx += vx[i]; my += vy[i]; mz += vz[i]; }
    double mnorm = std::sqrt(mx * mx + my * my + mz * mz);
    if (mnorm < 1e-10) return {};
    double ex = mx / mnorm, ey = my / mnorm, ez = mz / mnorm;

    // Inverse iteration: solve ATA * y = x, normalize, repeat.
    // Add a tiny shift to avoid singularity: (ATA + eps*I) * y = x.
    for (int iter = 0; iter < 50; iter++) {
        // Solve (ATA + eps*I) * y = e  via Cramer's rule
        double eps = 1e-12;
        double a00 = ATA[0][0] + eps, a01 = ATA[0][1], a02 = ATA[0][2];
        double a10 = ATA[1][0], a11 = ATA[1][1] + eps, a12 = ATA[1][2];
        double a20 = ATA[2][0], a21 = ATA[2][1], a22 = ATA[2][2] + eps;

        double det = a00 * (a11 * a22 - a12 * a21)
                   - a01 * (a10 * a22 - a12 * a20)
                   + a02 * (a10 * a21 - a11 * a20);
        if (std::fabs(det) < 1e-30) break;

        double yx = (ex * (a11 * a22 - a12 * a21)
                   - a01 * (ey * a22 - a12 * ez)
                   + a02 * (ey * a21 - a11 * ez)) / det;
        double yy = (a00 * (ey * a22 - a12 * ez)
                   - ex * (a10 * a22 - a12 * a20)
                   + a02 * (a10 * ez - ey * a20)) / det;
        double yz = (a00 * (a11 * ez - ey * a21)
                   - a01 * (a10 * ez - ey * a20)
                   + ex * (a10 * a21 - a11 * a20)) / det;

        double ynorm = std::sqrt(yx * yx + yy * yy + yz * yz);
        if (ynorm < 1e-15) break;
        ex = yx / ynorm; ey = yy / ynorm; ez = yz / ynorm;
    }

    // Choose the pole direction closest to the mean (not antipodal)
    if (ex * mx + ey * my + ez * mz < 0) {
        ex = -ex; ey = -ey; ez = -ez;
    }

    // Convert back to RA/Dec
    PoleEstimate est;
    est.dec = std::asin(std::clamp(static_cast<float>(ez), -1.0f, 1.0f));
    est.ra = std::atan2(static_cast<float>(ey), static_cast<float>(ex));
    if (est.ra < 0) est.ra += 2.0f * static_cast<float>(M_PI);
    est.valid = true;

    // Compute arc span of samples (max angular separation between any two)
    float max_ang = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            float dot = static_cast<float>(vx[i] * vx[j] + vy[i] * vy[j] + vz[i] * vz[j]);
            dot = std::clamp(dot, -1.0f, 1.0f);
            float ang = std::acos(dot);
            if (ang > max_ang) max_ang = ang;
        }
    }
    est.arc_deg = max_ang * 180.0f / static_cast<float>(M_PI);

    return est;
}

PoleEstimate PolarAligner::estimate_pole() const {
    if (samples_.size() < 3) return {};
    auto derotated = derotate_samples();
    return fit_circle(derotated);
}

double PolarAligner::lst_rad(double unix_time) const {
    // Julian date from Unix timestamp
    double jd = unix_time / 86400.0 + 2440587.5;
    double T = (jd - 2451545.0) / 36525.0;
    // GMST in degrees
    double gmst = std::fmod(280.46061837 + 360.98564736629 * (jd - 2451545.0)
                            + 0.000387933 * T * T, 360.0);
    if (gmst < 0) gmst += 360.0;
    // LST = GMST + observer longitude
    double lst_deg = gmst + observer_lon_rad_ * 180.0 / M_PI;
    return lst_deg * M_PI / 180.0;
}

void PolarAligner::eq_to_altaz(float ra, float dec, double lst,
                                float &alt, float &az) const {
    float ha = static_cast<float>(lst) - ra;
    float sin_alt = std::sin(observer_lat_rad_) * std::sin(dec)
                  + std::cos(observer_lat_rad_) * std::cos(dec) * std::cos(ha);
    alt = std::asin(std::clamp(sin_alt, -1.0f, 1.0f));

    float cos_alt = std::cos(alt);
    if (cos_alt < 1e-6f) {
        az = 0;
        return;
    }

    float sin_az = -std::cos(dec) * std::sin(ha) / cos_alt;
    float cos_az = (std::sin(dec) - std::sin(observer_lat_rad_) * sin_alt)
                 / (std::cos(observer_lat_rad_) * cos_alt);
    az = std::atan2(sin_az, cos_az);
}

AltAzOffset PolarAligner::pole_error(const PoleEstimate &pole) const {
    if (!pole.valid) return {};

    // Get current time for LST
    auto now = std::chrono::system_clock::now();
    double unix_time = std::chrono::duration<double>(
        now.time_since_epoch()).count();
    double lst = lst_rad(unix_time);

    // True celestial pole
    float true_dec = (observer_lat_rad_ >= 0)
                   ? static_cast<float>(M_PI / 2)
                   : static_cast<float>(-M_PI / 2);
    float true_ra = 0.0f; // RA is arbitrary for the pole, but use 0

    // Alt/az of true pole (should be alt=lat, az=0 for north or az=180 for south)
    float true_alt, true_az;
    eq_to_altaz(true_ra, true_dec, lst, true_alt, true_az);

    // Alt/az of estimated mount pole
    float mount_alt, mount_az;
    eq_to_altaz(pole.ra, pole.dec, lst, mount_alt, mount_az);

    AltAzOffset off;
    off.alt_arcmin = (mount_alt - true_alt) * 180.0f / static_cast<float>(M_PI) * 60.0f;

    // Az difference, handle wrap
    float daz = mount_az - true_az;
    if (daz > static_cast<float>(M_PI)) daz -= 2.0f * static_cast<float>(M_PI);
    if (daz < -static_cast<float>(M_PI)) daz += 2.0f * static_cast<float>(M_PI);
    off.az_arcmin = daz * 180.0f / static_cast<float>(M_PI) * 60.0f;

    off.total_arcmin = std::sqrt(off.alt_arcmin * off.alt_arcmin
                               + off.az_arcmin * off.az_arcmin);
    return off;
}
