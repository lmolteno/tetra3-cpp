#pragma once

#include <vector>
#include <cmath>
#include <cstdint>

// Observer location -- override at compile time or later via GNSS
#ifndef OBSERVER_LAT_DEG
#define OBSERVER_LAT_DEG -33.86 // Sydney
#endif
#ifndef OBSERVER_LON_DEG
#define OBSERVER_LON_DEG 151.21
#endif

enum class AppMode {
    Tracking,    // Normal plate-solve display
    PASampling,  // Polar alignment: collecting rotation samples
    PAFix        // Polar alignment: showing alt/az correction
};

struct PASample {
    float ra;          // radians, solved
    float dec;         // radians, solved
    double timestamp;  // seconds, monotonic (steady_clock)
};

struct PoleEstimate {
    float ra  = 0;     // radians
    float dec = 0;     // radians
    bool valid = false;
    float arc_deg = 0; // arc span of samples in degrees
};

struct AltAzOffset {
    float alt_arcmin = 0; // positive = mount pole too high
    float az_arcmin  = 0; // positive = mount pole too far east
    float total_arcmin = 0;
};

class PolarAligner {
public:
    PolarAligner(float observer_lat_deg, float observer_lon_deg);

    void clear_samples();
    void add_sample(float ra_rad, float dec_rad, double timestamp_s);
    int num_samples() const { return static_cast<int>(samples_.size()); }

    // Estimate the mount's polar axis from collected samples.
    // Needs >= 3 samples with sufficient arc span.
    PoleEstimate estimate_pole() const;

    // Compute alt/az offset between estimated pole and true celestial pole.
    AltAzOffset pole_error(const PoleEstimate &pole) const;

private:
    static constexpr double SIDEREAL_RATE = 7.2921159e-5; // rad/s

    float observer_lat_rad_;
    float observer_lon_rad_;
    std::vector<PASample> samples_;

    // De-rotate samples to the epoch of the first sample
    // (removes Earth's sidereal rotation from RA)
    std::vector<PASample> derotate_samples() const;

    // Fit a small circle through 3+ unit vectors, return its center (pole)
    PoleEstimate fit_circle(const std::vector<PASample> &derotated) const;

    // Compute Local Sidereal Time in radians from a Unix timestamp
    double lst_rad(double unix_time) const;

    // Convert equatorial (ra, dec) to alt/az for the observer at a given LST
    void eq_to_altaz(float ra, float dec, double lst,
                     float &alt, float &az) const;
};

// Live mount-pole tracker for PAFix mode.
//
// When PAFix starts, snapshot the fitted mount pole and the camera's current
// (ra, dec). Each subsequent solve gives a new camera direction. Under the
// assumption that the operator only turns the alt/az knobs (i.e. moves the
// whole mount rigidly, without touching RA/Dec axes), the same rigid 3D
// rotation that moved the camera also moved the mount pole. We recover that
// rotation as the unique no-twist rotation taking the predicted-no-move
// camera direction (just sidereally drifted) to the actual current camera
// direction, and apply it to the snapshotted pole.
//
// Operates entirely in celestial coordinates — no lat/lon needed for the
// update. Numerical alt/az still goes through PolarAligner::pole_error
// (which does need lat/lon).
class PAFixTracker {
public:
    static constexpr double SIDEREAL_RATE = 7.2921159e-5;  // rad/s

    void start(float pole_ra, float pole_dec,
               float cam_ra,  float cam_dec,
               double t_seconds);
    void stop();
    bool active() const { return active_; }

    // Live pole estimate at time t_now given the current camera direction.
    // Falls back to the snapshotted pole (just sidereally drifted) if the
    // camera is within numerical noise of the predicted-no-move position.
    PoleEstimate live_pole(float cam_ra_now, float cam_dec_now,
                           double t_seconds) const;

private:
    float  pole0_ra_   = 0;
    float  pole0_dec_  = 0;
    float  cam0_ra_    = 0;
    float  cam0_dec_   = 0;
    double t0_         = 0;
    bool   active_     = false;
};
