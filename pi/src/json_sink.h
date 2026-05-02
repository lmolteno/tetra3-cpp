#pragma once

#include "polar_align.h"
#include <tetra3/types.h>

#include <cstdint>
#include <cstddef>
#include <iosfwd>
#include <string>

// One-line NDJSON publisher. Writes the current daemon state plus a base64
// snapshot of the framebuffer. Consumers (web UI, I2C bridge) read it line by
// line from the daemon's stdout.
namespace json_sink {

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

    const std::uint8_t *fb = nullptr;
    std::size_t fb_size = 0;
};

void publish(std::ostream &out, const State &s);

}
