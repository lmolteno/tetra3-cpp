#pragma once

#include "frame_source.h"
#include <string>

class ImageFileSource : public IFrameSource {
public:
    explicit ImageFileSource(const std::string &path);
    bool initialize() override;
    std::optional<Frame> capture() override;

private:
    std::string file_path;
    bool loaded = false;
    Frame frame;
};
