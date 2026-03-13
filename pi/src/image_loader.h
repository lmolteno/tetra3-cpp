#pragma once

#include "frame_source.h"
#include <string>

class ImageFileSource : public IFrameSource {
public:
    explicit ImageFileSource(const std::string &path, bool verbose = false);
    bool initialize() override;
    std::optional<Frame> capture() override;

private:
    std::string file_path;
    bool verbose = false;
    bool loaded = false;
    Frame frame;
};
