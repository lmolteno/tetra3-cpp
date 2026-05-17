#pragma once

#include "app/polar_align.h"
#include <tetra3/types.h>

#include <cstdint>
#include <cstddef>
#include <iosfwd>
#include <string>

// One-line NDJSON publisher. Writes the current daemon state plus a base64
// snapshot of the framebuffer. Consumers (web UI, I2C bridge) read it line by
// line from the daemon's stdout.
namespace json_sink {

// Per-tick latency breakdown. Negative values are omitted from output. All
// times in milliseconds. capture / archive are this-tick wall clocks; load /
// detect / match come from solve_cli; solve_wall is pi_tracker's
// submit→result-arrived wall clock; result_age is how stale the displayed
// result is at publish time.
struct Timing {
    float capture     = -1.0f;   // main-thread mailbox poll cost (was wait+bin pre-decouple)
    float wait        = -1.0f;   // capture thread: time between completed requests
    float bin         = -1.0f;   // capture thread: 2x2 Bayer averaging
    float frame_age   = -1.0f;   // age of the displayed frame at publish time
    float archive     = -1.0f;
    float load        = -1.0f;
    float detect      = -1.0f;
    float bg          = -1.0f;   // detect sub-stage: tiled BG estimation
    float mask        = -1.0f;   // detect sub-stage: bilinear threshold mask
    float cluster     = -1.0f;   // detect sub-stage: flood fill + centroiding
    float match       = -1.0f;
    float solve_wall  = -1.0f;
    float result_age  = -1.0f;
};

struct State {
    AppMode mode = AppMode::Tracking;
    SolveResult result{};

    // Polar-align state. Ignored unless mode is PASampling or PAFix.
    int pa_samples = 0;
    PoleEstimate pole{};
    AltAzOffset offset{};

    // Optional fields. Sentinel values omit them from output.
    float    lens_pos    = -1.0f;     // negative  -> unknown
    uint32_t exposure_us = 0;         // 0         -> unknown
    float    gain        = -1.0f;     // negative  -> unknown
    std::string frame_path;            // empty     -> no archived frame

    Timing timing;

    // IMU pointing. imu_valid=false → omitted from output.
    bool  imu_valid      = false;
    float imu_ra         = 0.0f;
    float imu_dec        = 0.0f;
    float imu_roll       = 0.0f;
    bool  imu_calibrated = false;

    const std::uint8_t *fb = nullptr;
    std::size_t fb_size = 0;
};

void publish(std::ostream &out, const State &s);

}
