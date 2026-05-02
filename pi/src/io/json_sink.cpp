#include "io/json_sink.h"
#include <ostream>

namespace json_sink {

namespace {

const char B64[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

void write_base64(std::ostream &out, const std::uint8_t *data, std::size_t len) {
    for (std::size_t i = 0; i < len; i += 3) {
        std::uint32_t n = static_cast<std::uint32_t>(data[i]) << 16;
        if (i + 1 < len) n |= static_cast<std::uint32_t>(data[i + 1]) << 8;
        if (i + 2 < len) n |= data[i + 2];
        out << B64[(n >> 18) & 63];
        out << B64[(n >> 12) & 63];
        out << ((i + 1 < len) ? B64[(n >> 6) & 63] : '=');
        out << ((i + 2 < len) ? B64[n & 63] : '=');
    }
}

const char *mode_str(AppMode m) {
    switch (m) {
    case AppMode::Tracking:   return "tracking";
    case AppMode::PASampling: return "pa_sampling";
    case AppMode::PAFix:      return "pa_fix";
    }
    return "tracking";
}

}

void publish(std::ostream &out, const State &s) {
    out << "{\"mode\":\"" << mode_str(s.mode) << "\"";

    out << ",\"solved\":" << (s.result.solved ? "true" : "false");
    if (s.result.solved) {
        out << ",\"ra_rad\":"        << s.result.ra
            << ",\"dec_rad\":"       << s.result.dec
            << ",\"roll_rad\":"      << s.result.roll
            << ",\"fov_rad\":"       << s.result.fov
            << ",\"rmse\":"          << s.result.rmse
            << ",\"num_matches\":"   << s.result.num_matches
            << ",\"solve_time_ms\":" << s.result.solve_time_ms;
    }

    if (s.mode == AppMode::PASampling || s.mode == AppMode::PAFix) {
        out << ",\"pa\":{"
            << "\"num_samples\":" << s.pa_samples
            << ",\"pole_valid\":" << (s.pole.valid ? "true" : "false");
        if (s.pole.valid) {
            out << ",\"pole_ra_rad\":"  << s.pole.ra
                << ",\"pole_dec_rad\":" << s.pole.dec
                << ",\"arc_deg\":"      << s.pole.arc_deg;
        }
        if (s.mode == AppMode::PAFix) {
            out << ",\"alt_arcmin\":"   << s.offset.alt_arcmin
                << ",\"az_arcmin\":"    << s.offset.az_arcmin
                << ",\"total_arcmin\":" << s.offset.total_arcmin;
        }
        out << "}";
    }

    if (s.lens_pos >= 0.0f) {
        out << ",\"lens_pos\":" << s.lens_pos;
    }
    if (s.exposure_us > 0) {
        out << ",\"exposure_us\":" << s.exposure_us;
    }
    if (s.gain >= 0.0f) {
        out << ",\"gain\":" << s.gain;
    }

    if (!s.frame_path.empty()) {
        out << ",\"frame_path\":\"" << s.frame_path << "\"";
    }

    if (s.fb && s.fb_size > 0) {
        out << ",\"fb\":\"";
        write_base64(out, s.fb, s.fb_size);
        out << "\"";
    }

    out << "}\n";
    out.flush();
}

}
