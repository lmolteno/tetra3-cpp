#include "capture/image_loader.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#define STBI_ONLY_JPEG
#define STBI_ONLY_PNM
#define STBI_ONLY_TGA
#define STBI_ONLY_BMP
#define STBI_NO_HDR
#define STBI_NO_LINEAR
#include <stb_image.h>

ImageFileSource::ImageFileSource(const std::string &path, bool verbose,
                                  bool watch)
    : file_path(path), verbose(verbose), watch_(watch) {}

static int64_t file_mtime_ns(const std::string &path) {
    struct stat st{};
    if (::stat(path.c_str(), &st) != 0) return 0;
#if defined(__APPLE__)
    return static_cast<int64_t>(st.st_mtimespec.tv_sec) * 1000000000LL
         + st.st_mtimespec.tv_nsec;
#else
    return static_cast<int64_t>(st.st_mtim.tv_sec) * 1000000000LL
         + st.st_mtim.tv_nsec;
#endif
}

namespace {

// stb_image's PNM loader only handles 8-bit PGM. We see 16-bit PGM frequently
// (synthetic images, raw dumps), so we have a small parser of our own.
bool load_pgm_16(const std::string &path, std::vector<uint16_t> &out,
                 int &w, int &h, int &bit_depth) {
    FILE *f = std::fopen(path.c_str(), "rb");
    if (!f) return false;

    auto skip_ws = [&]() {
        for (;;) {
            int c = std::fgetc(f);
            if (c == EOF) return;
            if (c == '#') {
                while (c != '\n' && c != EOF) c = std::fgetc(f);
                continue;
            }
            if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
            std::ungetc(c, f);
            return;
        }
    };
    auto read_int = [&]() -> int {
        skip_ws();
        int v = 0;
        bool any = false;
        for (;;) {
            int c = std::fgetc(f);
            if (c < '0' || c > '9') {
                if (c != EOF) std::ungetc(c, f);
                break;
            }
            v = v * 10 + (c - '0');
            any = true;
        }
        return any ? v : -1;
    };

    char magic[2];
    if (std::fread(magic, 1, 2, f) != 2 || magic[0] != 'P' || magic[1] != '5') {
        std::fclose(f);
        return false;
    }
    int width = read_int();
    int height = read_int();
    int maxval = read_int();
    if (width <= 0 || height <= 0 || maxval <= 0) {
        std::fclose(f);
        return false;
    }
    // Single whitespace byte after maxval, per spec.
    std::fgetc(f);

    out.resize(static_cast<size_t>(width) * height);
    if (maxval <= 255) {
        std::vector<uint8_t> tmp(out.size());
        if (std::fread(tmp.data(), 1, tmp.size(), f) != tmp.size()) {
            std::fclose(f);
            return false;
        }
        for (size_t i = 0; i < tmp.size(); i++) {
            out[i] = static_cast<uint16_t>(tmp[i]) * 257;  // expand to full 16-bit
        }
        bit_depth = 8;
    } else {
        // 16-bit big-endian per PGM spec.
        std::vector<uint8_t> tmp(out.size() * 2);
        if (std::fread(tmp.data(), 1, tmp.size(), f) != tmp.size()) {
            std::fclose(f);
            return false;
        }
        for (size_t i = 0; i < out.size(); i++) {
            out[i] = (static_cast<uint16_t>(tmp[2*i]) << 8) | tmp[2*i + 1];
        }
        bit_depth = 16;
    }
    w = width;
    h = height;
    std::fclose(f);
    return true;
}

}

bool ImageFileSource::do_load() {
    // Try our PGM-16 fast path first (covers the common synthetic / raw-dump
    // case).
    int w = 0, h = 0, bit_depth = 0;
    std::vector<uint16_t> pgm;
    if (load_pgm_16(file_path, pgm, w, h, bit_depth)) {
        frame.width = w;
        frame.height = h;
        frame.bit_depth = bit_depth;
        frame.data = std::move(pgm);
        loaded = true;
        if (verbose) {
            std::cout << "Loaded image: " << file_path << " (" << w << "x" << h
                      << ", " << bit_depth << "-bit PGM)" << std::endl;
        }
        return true;
    }

    // Fall back to stb_image: 16-bit-aware load for PNG, otherwise 8-bit. We
    // force a single grayscale channel; stb does the colour conversion.
    int channels = 0;
    if (stbi_is_16_bit(file_path.c_str())) {
        stbi_us *px = stbi_load_16(file_path.c_str(), &w, &h, &channels, 1);
        if (!px) {
            std::cerr << "Failed to load image: " << file_path
                      << " (" << stbi_failure_reason() << ")" << std::endl;
            return false;
        }
        frame.data.assign(px, px + static_cast<size_t>(w) * h);
        stbi_image_free(px);
        bit_depth = 16;
    } else {
        stbi_uc *px = stbi_load(file_path.c_str(), &w, &h, &channels, 1);
        if (!px) {
            std::cerr << "Failed to load image: " << file_path
                      << " (" << stbi_failure_reason() << ")" << std::endl;
            return false;
        }
        frame.data.resize(static_cast<size_t>(w) * h);
        for (size_t i = 0; i < frame.data.size(); i++) {
            frame.data[i] = static_cast<uint16_t>(px[i]) * 257;  // 8-bit -> 16-bit
        }
        stbi_image_free(px);
        bit_depth = 8;
    }

    frame.width = w;
    frame.height = h;
    frame.bit_depth = bit_depth;
    loaded = true;
    if (verbose) {
        std::cout << "Loaded image: " << file_path << " (" << w << "x" << h
                  << ", " << bit_depth << "-bit)" << std::endl;
    }
    return true;
}

bool ImageFileSource::initialize() {
    if (!do_load()) return false;
    last_mtime_ns_ = file_mtime_ns(file_path);
    delivered_ = false;
    return true;
}

std::optional<Frame> ImageFileSource::capture() {
    if (!loaded) return std::nullopt;

    if (!watch_) {
        // Original behaviour: always return the loaded frame.
        return frame;
    }

    // Watch mode: deliver the loaded frame once, then only re-deliver when
    // the file's mtime advances (the simulator/test writes a new image).
    if (!delivered_) {
        delivered_ = true;
        return frame;
    }

    int64_t mt = file_mtime_ns(file_path);
    if (mt == 0 || mt == last_mtime_ns_) {
        // No update yet. Return nullopt so the daemon's "no fresh frame"
        // path runs — render the last solve, don't resubmit.
        return std::nullopt;
    }

    if (!do_load()) {
        // Partial write or transient error — try again next tick.
        return std::nullopt;
    }
    last_mtime_ns_ = mt;
    return frame;
}
