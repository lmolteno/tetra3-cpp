#include "app/cam_state.h"

#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

namespace cam_state {

namespace {

std::string path_for(const std::string &dir, const std::string &name) {
    return dir + "/" + name;
}

template <typename T>
bool save_value(const std::string &dir, const std::string &name, T v) {
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        std::cerr << "cam_state: create_directories(" << dir << ") failed: "
                  << ec.message() << "\n";
        return false;
    }

    std::string final_path = path_for(dir, name);
    std::string tmp_path   = final_path + ".tmp";

    {
        std::ofstream out(tmp_path, std::ios::trunc);
        if (!out) {
            std::cerr << "cam_state: open(" << tmp_path << ") failed\n";
            return false;
        }
        out << v << "\n";
        if (!out) {
            std::cerr << "cam_state: write to " << tmp_path << " failed\n";
            return false;
        }
    }

    fs::rename(tmp_path, final_path, ec);
    if (ec) {
        std::cerr << "cam_state: rename failed: " << ec.message() << "\n";
        return false;
    }
    return true;
}

template <typename T>
T load_value(const std::string &dir, const std::string &name, T def) {
    std::ifstream in(path_for(dir, name));
    if (!in) return def;
    T v;
    in >> v;
    if (!in) return def;
    return v;
}

}

Settings load(const std::string &state_dir) {
    Settings s;
    s.lens_pos    = load_value<float>   (state_dir, "lens_position", DEFAULT_LENS_POS);
    s.exposure_us = load_value<uint32_t>(state_dir, "exposure_us",   DEFAULT_EXPOSURE_US);
    s.gain        = load_value<float>   (state_dir, "gain",          DEFAULT_GAIN);
    return s;
}

bool save_lens_pos   (const std::string &dir, float v)    { return save_value(dir, "lens_position", v); }
bool save_exposure_us(const std::string &dir, uint32_t v) { return save_value(dir, "exposure_us",   v); }
bool save_gain       (const std::string &dir, float v)    { return save_value(dir, "gain",          v); }

}
