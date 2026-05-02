#include "capture/frame_writer.h"

#include <cstdio>
#include <filesystem>
#include <iostream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace fs = std::filesystem;

FrameWriter::FrameWriter(std::string dir, int slots)
    : dir_(std::move(dir)), slots_(slots) {}

bool FrameWriter::ensure_dir() const {
    std::error_code ec;
    fs::create_directories(dir_, ec);
    if (ec) {
        std::cerr << "FrameWriter: create_directories(" << dir_ << ") failed: "
                  << ec.message() << "\n";
        return false;
    }
    return true;
}

std::string FrameWriter::slot_path(int idx) const {
    char name[16];
    std::snprintf(name, sizeof(name), "%02d.png", idx);
    return dir_ + "/" + name;
}

std::string FrameWriter::write(const Frame &frame) {
    if (frame.data.empty()) return {};

    int idx = next_;
    next_ = (next_ + 1) % slots_;

    // stb_image_write only supports 8-bit PNG output. Downscale by the
    // frame's actual bit depth so 8-bit values use the full 0..255 range
    // regardless of whether the camera produced 10/12/16-bit data.
    int shift = frame.bit_depth - 8;
    if (shift < 0) shift = 0;

    std::vector<uint8_t> px(static_cast<size_t>(frame.width) * frame.height);
    for (size_t i = 0; i < px.size(); i++) {
        uint16_t v = frame.data[i] >> shift;
        px[i] = static_cast<uint8_t>(v > 255 ? 255 : v);
    }

    std::string path = slot_path(idx);
    int stride = frame.width;  // 1 byte per pixel, 1 channel
    if (!stbi_write_png(path.c_str(), frame.width, frame.height, 1,
                        px.data(), stride)) {
        std::cerr << "FrameWriter: stbi_write_png(" << path << ") failed\n";
        return {};
    }
    return path;
}
