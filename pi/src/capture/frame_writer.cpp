#include "capture/frame_writer.h"

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <vector>

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
    std::snprintf(name, sizeof(name), "%02d.pgm", idx);
    return dir_ + "/" + name;
}

std::string FrameWriter::write(const Frame &frame) {
    if (frame.data.empty()) return {};

    int idx = next_;
    next_ = (next_ + 1) % slots_;

    // 16-bit binary PGM (P5). Spec is big-endian, so byte-swap the uint16s
    // into a temp buffer then issue one fwrite — about as fast as a memcpy on
    // tmpfs. solve_cli picks this up via the PGM-16 fast path in image_loader
    // (no zlib, no codec, full bit depth preserved).
    std::vector<std::uint8_t> bytes(frame.data.size() * 2);
    for (std::size_t i = 0; i < frame.data.size(); i++) {
        bytes[2 * i]     = static_cast<std::uint8_t>(frame.data[i] >> 8);
        bytes[2 * i + 1] = static_cast<std::uint8_t>(frame.data[i] & 0xFF);
    }

    std::string path = slot_path(idx);
    FILE *f = std::fopen(path.c_str(), "wb");
    if (!f) {
        std::cerr << "FrameWriter: fopen(" << path << ") failed\n";
        return {};
    }
    if (std::fprintf(f, "P5\n%d %d\n65535\n", frame.width, frame.height) < 0
        || std::fwrite(bytes.data(), 1, bytes.size(), f) != bytes.size()) {
        std::cerr << "FrameWriter: write to " << path << " failed\n";
        std::fclose(f);
        return {};
    }
    std::fclose(f);
    return path;
}
