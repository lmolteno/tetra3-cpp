#pragma once

#ifdef HAS_LIBCAMERA

#include "capture/frame_source.h"
#include <atomic>
#include <memory>

class LibcameraCapture : public IFrameSource {
public:
    LibcameraCapture(int width = 2028, int height = 1520);
    ~LibcameraCapture();

    bool initialize() override;
    std::optional<Frame> capture() override;

    // Manual AF / exposure / gain only — no AE, no AWB. Set any of these
    // before initialize() to seed the starting values; afterwards updates are
    // applied on the next request.
    void set_lens_position(float diopters);
    void set_exposure_us(uint32_t us);
    void set_gain(float multiplier);

    float    lens_position() const { return lens_position_.load(); }
    uint32_t exposure_us()   const { return exposure_us_.load();   }
    float    gain()          const { return gain_.load();          }

    // Per-call breakdown of the most recent capture(): wait_ms is the time
    // blocked on libcamera (i.e. exposure + readout + ISP), bin_ms is the
    // 2x2 Bayer averaging in software. -1 if not yet measured.
    float last_wait_ms() const { return last_wait_ms_.load(); }
    float last_bin_ms()  const { return last_bin_ms_.load();  }

private:
    int req_width, req_height;
    std::atomic<float>    lens_position_{0.0f};
    std::atomic<uint32_t> exposure_us_{100'000};   // 100 ms
    std::atomic<float>    gain_{1.0f};
    std::atomic<float>    last_wait_ms_{-1.0f};
    std::atomic<float>    last_bin_ms_{-1.0f};

    struct Impl;
    std::unique_ptr<Impl> impl;
};

#endif // HAS_LIBCAMERA
