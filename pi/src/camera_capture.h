#pragma once

#ifdef HAS_LIBCAMERA

#include "frame_source.h"
#include <memory>

class LibcameraCapture : public IFrameSource {
public:
    LibcameraCapture(int width = 2028, int height = 1520);
    ~LibcameraCapture();
    bool initialize() override;
    std::optional<Frame> capture() override;

private:
    int req_width, req_height;
    struct Impl;
    std::unique_ptr<Impl> impl;
};

#endif // HAS_LIBCAMERA
