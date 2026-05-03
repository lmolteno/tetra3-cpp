#pragma once

#ifdef HAS_LIBCAMERA

#include "capture/frame_source.h"
#include <atomic>
#include <chrono>
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

    // Per-frame breakdown of the most recent frame produced by the internal
    // capture thread: wait_ms is the time the thread sat blocked on
    // libcamera (i.e. exposure + readout + ISP, modulo our queue depth),
    // bin_ms is the 2x2 Bayer averaging in software. -1 if not yet measured.
    float last_wait_ms() const { return last_wait_ms_.load(); }
    float last_bin_ms()  const { return last_bin_ms_.load();  }

private:
    int req_width, req_height;
    std::atomic<float>    lens_position_{0.0f};
    std::atomic<uint32_t> exposure_us_{100'000};   // 100 ms
    std::atomic<float>    gain_{1.0f};
    std::atomic<float>    last_wait_ms_{-1.0f};
    std::atomic<float>    last_bin_ms_{-1.0f};

    // Used by capture_thread_main() to compute wait_ms between request pops.
    std::chrono::steady_clock::time_point last_request_pop_{};

    void capture_thread_main();

    struct Impl;
    std::unique_ptr<Impl> impl;
};

#endif // HAS_LIBCAMERA
