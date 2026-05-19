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
//   {"solve": {                  inject a fake solve result — used by the
//      "ra_rad":  <float>,        PA simulator / accuracy test tools to drive
//      "dec_rad": <float>,        the daemon without running solve_cli. Pair
//      "roll_rad":<float>,        with --no-solver. ra/dec/roll/fov required;
//      "fov_rad": <float>,        rmse and matches default to 0/4. Set
//      "rmse":    <float>,        "solved":false to inject a no-solve frame.
//      "matches": <int>,
//      "solved":  <bool>}}
// Multiple keys per line are allowed and treated as a single command.
class ControlInput {
public:
    struct InjectedSolve {
        bool   solved = true;
        float  ra_rad = 0, dec_rad = 0, roll_rad = 0, fov_rad = 0;
        float  rmse   = 0;
        int    num_matches = 0;
    };

    struct Command {
        std::optional<AppMode>       mode;
        std::optional<float>         lens_pos;
        std::optional<uint32_t>      exposure_us;
        std::optional<float>         gain;
        std::optional<InjectedSolve> solve;
    };

    ControlInput();

    // Returns commands parsed since last call. Empty vector if nothing.
    std::vector<Command> poll();

private:
    std::string buf_;
    bool eof_ = false;

    static std::optional<Command> parse_line(const std::string &line);
};
