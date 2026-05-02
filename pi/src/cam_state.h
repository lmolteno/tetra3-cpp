#pragma once

#include <cstdint>
#include <string>

// Persistent manual camera settings: lens position, exposure time, analogue
// gain. Stored as one text file per field under <state_dir> (atomic via
// tmp+rename) so a tweak to one value doesn't risk corrupting the others.
namespace cam_state {

constexpr float    DEFAULT_LENS_POS    = 0.0f;       // diopters; 0 = infinity
constexpr uint32_t DEFAULT_EXPOSURE_US = 100'000;    // 100 ms
constexpr float    DEFAULT_GAIN        = 1.0f;       // analogue gain multiplier

struct Settings {
    float    lens_pos    = DEFAULT_LENS_POS;
    uint32_t exposure_us = DEFAULT_EXPOSURE_US;
    float    gain        = DEFAULT_GAIN;
};

// Best-effort: missing or unreadable files fall back to defaults.
Settings load(const std::string &state_dir);

bool save_lens_pos   (const std::string &state_dir, float v);
bool save_exposure_us(const std::string &state_dir, uint32_t v);
bool save_gain       (const std::string &state_dir, float v);

}
