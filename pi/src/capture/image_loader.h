#pragma once

#include "capture/frame_source.h"
#include <cstdint>
#include <string>

class ImageFileSource : public IFrameSource {
public:
    // watch = true: capture() polls the file's mtime each tick. The first
    // capture returns the initially loaded frame; subsequent captures return
    // a freshly reloaded frame whenever the file is rewritten, and nullopt
    // when nothing has changed since the last delivery. Default (watch=false)
    // returns the loaded frame indefinitely — the original test-image
    // behaviour.
    explicit ImageFileSource(const std::string &path, bool verbose = false,
                             bool watch = false);
    bool initialize() override;
    std::optional<Frame> capture() override;

private:
    bool do_load();

    std::string file_path;
    bool verbose = false;
    bool watch_ = false;
    bool loaded = false;
    bool delivered_ = false;     // true once we've returned the current frame
    int64_t last_mtime_ns_ = 0;  // mtime in nanoseconds; 0 = unset
    Frame frame;
};
