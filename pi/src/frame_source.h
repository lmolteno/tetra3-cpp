#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

struct Frame {
    std::vector<uint16_t> data; // 16-bit grayscale
    int width, height, bit_depth;
};

class IFrameSource {
public:
    virtual ~IFrameSource() = default;
    virtual bool initialize() = 0;
    virtual std::optional<Frame> capture() = 0;
};
