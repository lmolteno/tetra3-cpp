#pragma once

#include "app/polar_align.h"

#include <optional>
#include <string>
#include <vector>

// Reads NDJSON commands from stdin without blocking the main loop. Each
// poll() drains whatever has arrived since the last call and returns the
// parsed commands.
//
// Recognized keys per line:
//   {"mode": "tracking" | "pa_sampling" | "pa_fix"}
//   {"lens_pos":    <float>}     diopters
//   {"exposure_us": <int>}       microseconds
//   {"gain":        <float>}     analogue-gain multiplier
// Multiple keys per line are allowed and treated as a single command.
class ControlInput {
public:
    struct Command {
        std::optional<AppMode>  mode;
        std::optional<float>    lens_pos;
        std::optional<uint32_t> exposure_us;
        std::optional<float>    gain;
    };

    ControlInput();

    // Returns commands parsed since last call. Empty vector if nothing.
    std::vector<Command> poll();

private:
    std::string buf_;
    bool eof_ = false;

    static std::optional<Command> parse_line(const std::string &line);
};
