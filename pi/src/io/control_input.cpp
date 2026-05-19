#include "io/control_input.h"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>

ControlInput::ControlInput() {
    int flags = fcntl(STDIN_FILENO, F_GETFL, 0);
    if (flags >= 0) {
        fcntl(STDIN_FILENO, F_SETFL, flags | O_NONBLOCK);
    }
}

std::vector<ControlInput::Command> ControlInput::poll() {
    std::vector<Command> out;
    if (eof_) return out;

    char chunk[1024];
    while (true) {
        ssize_t n = read(STDIN_FILENO, chunk, sizeof(chunk));
        if (n > 0) {
            buf_.append(chunk, static_cast<std::size_t>(n));
            continue;
        }
        if (n == 0) {
            eof_ = true;
            break;
        }
        if (errno == EAGAIN || errno == EWOULDBLOCK) break;
        if (errno == EINTR) continue;
        eof_ = true;
        break;
    }

    std::size_t start = 0;
    for (std::size_t i = 0; i < buf_.size(); i++) {
        if (buf_[i] != '\n') continue;
        std::string line = buf_.substr(start, i - start);
        start = i + 1;
        if (line.empty()) continue;
        if (auto cmd = parse_line(line)) out.push_back(*cmd);
    }
    if (start > 0) buf_.erase(0, start);

    return out;
}

namespace {

std::string find_string(const std::string &s, const char *key) {
    std::string needle = std::string("\"") + key + "\"";
    auto pos = s.find(needle);
    if (pos == std::string::npos) return {};
    pos = s.find(':', pos);
    if (pos == std::string::npos) return {};
    auto q1 = s.find('"', pos + 1);
    if (q1 == std::string::npos) return {};
    auto q2 = s.find('"', q1 + 1);
    if (q2 == std::string::npos) return {};
    return s.substr(q1 + 1, q2 - q1 - 1);
}

std::optional<float> find_float(const std::string &s, const char *key) {
    std::string needle = std::string("\"") + key + "\"";
    auto pos = s.find(needle);
    if (pos == std::string::npos) return std::nullopt;
    pos = s.find(':', pos);
    if (pos == std::string::npos) return std::nullopt;
    char *end = nullptr;
    const char *start = s.c_str() + pos + 1;
    while (*start == ' ' || *start == '\t') start++;
    if (*start == '"' || *start == '{' || *start == '[') return std::nullopt;
    float v = std::strtof(start, &end);
    if (end == start) return std::nullopt;
    return v;
}

}

std::optional<ControlInput::Command> ControlInput::parse_line(const std::string &line) {
    Command cmd;

    std::string mode = find_string(line, "mode");
    if (!mode.empty()) {
        if      (mode == "tracking")    cmd.mode = AppMode::Tracking;
        else if (mode == "pa_sampling") cmd.mode = AppMode::PASampling;
        else if (mode == "pa_fix")      cmd.mode = AppMode::PAFix;
    }

    cmd.lens_pos = find_float(line, "lens_pos");
    cmd.gain     = find_float(line, "gain");
    if (auto v = find_float(line, "exposure_us")) {
        if (*v >= 0) cmd.exposure_us = static_cast<uint32_t>(*v);
    }

    // Look for a "solve" key. We only inject a solve if at least one of the
    // numeric fields is present — otherwise the inner search just sees the
    // outer "mode" key and produces a bogus inject.
    if (line.find("\"solve\"") != std::string::npos) {
        InjectedSolve s;
        bool any = false;
        if (auto v = find_float(line, "ra_rad"))   { s.ra_rad   = *v; any = true; }
        if (auto v = find_float(line, "dec_rad"))  { s.dec_rad  = *v; any = true; }
        if (auto v = find_float(line, "roll_rad")) { s.roll_rad = *v; any = true; }
        if (auto v = find_float(line, "fov_rad"))  { s.fov_rad  = *v; any = true; }
        if (auto v = find_float(line, "rmse"))     { s.rmse     = *v; }
        if (auto v = find_float(line, "matches"))  { s.num_matches = static_cast<int>(*v); }
        // "solved":false explicitly injects a no-solve frame.
        if (line.find("\"solved\":false") != std::string::npos) {
            s.solved = false;
            any = true;
        }
        if (any) cmd.solve = s;
    }

    if (!cmd.mode && !cmd.lens_pos && !cmd.exposure_us
        && !cmd.gain && !cmd.solve) return std::nullopt;
    return cmd;
}
