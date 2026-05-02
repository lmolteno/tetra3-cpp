#pragma once

#include "frame_source.h"
#include <string>

// Writes incoming frames as 16-bit PNGs into a directory, rotating through a
// fixed slot count (default 10). Path is meant to live in tmpfs so this
// doesn't burn the SD card.
class FrameWriter {
public:
    explicit FrameWriter(std::string dir, int slots = 10);

    // Returns true if the directory exists or was created.
    bool ensure_dir() const;

    // Write `frame` into the next slot and return its absolute path. Empty
    // string on failure.
    std::string write(const Frame &frame);

private:
    std::string dir_;
    int slots_;
    int next_ = 0;

    std::string slot_path(int idx) const;
};
