#ifdef HAS_LIBCAMERA

#include "capture/camera_capture.h"
#include <array>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
#include <sys/mman.h>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif
#include <libcamera/libcamera.h>
#include <libcamera/control_ids.h>
#include <libcamera/formats.h>
#include <libcamera/base/object.h>

using namespace libcamera;

// Inherits from libcamera::Object so it can be the receiver of the
// requestCompleted signal in libcamera >= 0.2 (the API now requires an Object
// receiver instead of accepting bare lambdas).
struct LibcameraCapture::Impl : public Object {
    std::unique_ptr<CameraManager> cm;
    std::shared_ptr<Camera> camera;
    std::unique_ptr<CameraConfiguration> config;
    std::unique_ptr<FrameBufferAllocator> allocator;
    std::vector<std::unique_ptr<Request>> requests;
    Stream *stream = nullptr;

    // FIFO of completed requests pushed by libcamera's worker thread, popped
    // by the internal capture thread. We keep all allocated requests queued
    // continuously (re-queued the moment they come back) so the camera always
    // has work to do — no missed sensor cycles between consumer reads.
    std::mutex mtx;
    std::condition_variable cv;
    std::deque<Request *> completed;

    // Last control values we actually wrote into a Request. Compared against
    // the atomics on every re-queue so we only set controls when something
    // actually changed — Pi's IPA treats any control-on-request as a change
    // and lags new values by ~2 frames before they take effect.
    float    applied_lens        = std::nanf("");
    uint32_t applied_exposure_us = 0;
    float    applied_gain        = std::nanf("");

    // Internal capture thread: continuously consumes completed requests, runs
    // 2x2 bin, and stashes the most recent Frame in the mailbox. capture() is
    // a non-blocking poll on this mailbox.
    std::thread thread;
    std::atomic<bool> running{false};

    std::mutex mailbox_mu;
    std::optional<Frame> latest;
    std::uint64_t latest_seq = 0;
    std::uint64_t consumed_seq = 0;

    void on_request_completed(Request *request) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            completed.push_back(request);
        }
        cv.notify_one();
    }
};

LibcameraCapture::LibcameraCapture(int width, int height)
    : req_width(width), req_height(height), impl(std::make_unique<Impl>()) {}

LibcameraCapture::~LibcameraCapture() {
    // Stop the capture thread first so it isn't trying to use camera/buffers
    // while libcamera is tearing them down.
    if (impl->running.load()) {
        impl->running.store(false);
        impl->cv.notify_all();
        if (impl->thread.joinable()) impl->thread.join();
    }
    if (impl->camera) {
        impl->camera->stop();
        impl->camera->release();
    }
    if (impl->cm) {
        impl->cm->stop();
    }
}

bool LibcameraCapture::initialize() {
    impl->cm = std::make_unique<CameraManager>();
    int ret = impl->cm->start();
    if (ret) {
        std::cerr << "Failed to start CameraManager: " << ret << std::endl;
        return false;
    }

    auto cameras = impl->cm->cameras();
    if (cameras.empty()) {
        std::cerr << "No cameras found" << std::endl;
        return false;
    }

    impl->camera = cameras[0];
    ret = impl->camera->acquire();
    if (ret) {
        std::cerr << "Failed to acquire camera" << std::endl;
        return false;
    }

    // Raw Bayer stream. The default libcamera picks for the Pi is the CSI-2
    // packed 10-bit format (5 bytes per 4 pixels), which is awkward to read.
    // Force the unpacked variant so each pixel sits in one little-endian
    // uint16 with the value in the low 10 bits.
    impl->config = impl->camera->generateConfiguration({StreamRole::Raw});
    if (!impl->config) {
        std::cerr << "Failed to generate camera configuration" << std::endl;
        return false;
    }

    StreamConfiguration &stream_config = impl->config->at(0);
    stream_config.size.width = req_width * 2;   // Raw is full Bayer, we'll downsample
    stream_config.size.height = req_height * 2;
    stream_config.pixelFormat = formats::SBGGR10;  // IMX708/477/219 are all BGGR

    CameraConfiguration::Status status = impl->config->validate();
    if (status == CameraConfiguration::Invalid) {
        std::cerr << "Camera configuration invalid" << std::endl;
        return false;
    }
    if (stream_config.pixelFormat != formats::SBGGR10) {
        std::cerr << "Warning: requested SBGGR10 but got "
                  << stream_config.pixelFormat.toString()
                  << " — the 2x2 Bayer averaging assumes 1 px per uint16\n";
    }

    ret = impl->camera->configure(impl->config.get());
    if (ret) {
        std::cerr << "Failed to configure camera: " << ret << std::endl;
        return false;
    }

    impl->stream = impl->config->at(0).stream();

    // Allocate buffers
    impl->allocator = std::make_unique<FrameBufferAllocator>(impl->camera);
    ret = impl->allocator->allocate(impl->stream);
    if (ret < 0) {
        std::cerr << "Failed to allocate buffers" << std::endl;
        return false;
    }

    // Create requests
    for (const auto &buffer : impl->allocator->buffers(impl->stream)) {
        auto request = impl->camera->createRequest();
        if (!request) {
            std::cerr << "Failed to create request" << std::endl;
            return false;
        }
        request->addBuffer(impl->stream, buffer.get());
        impl->requests.push_back(std::move(request));
    }

    // Connect with ConnectionTypeDirect so the slot runs on libcamera's worker
    // thread (no event loop required on the main thread).
    impl->camera->requestCompleted.connect(
        impl.get(), &Impl::on_request_completed, ConnectionTypeDirect);

    // Fully manual: AE off (so ExposureTime + AnalogueGain stick), AF off
    // (LensPosition is the source of truth). FrameDurationLimits pinned to the
    // exposure time so the pipeline doesn't try to clamp it. Initial values
    // are set here and tracked in applied_* so apply_pending_controls() knows
    // what's already in effect and only writes deltas.
    const uint32_t initial_exposure = exposure_us_.load();
    const float    initial_gain     = gain_.load();
    const float    initial_lens     = lens_position_.load();
    const std::array<int64_t, 2> fdur = { initial_exposure, initial_exposure };

    ControlList start_controls(impl->camera->controls());
    start_controls.set(controls::AeEnable, false);
    start_controls.set(controls::AfMode,        controls::AfModeManual);
    start_controls.set(controls::ExposureTime,  static_cast<int32_t>(initial_exposure));
    start_controls.set(controls::AnalogueGain,  initial_gain);
    start_controls.set(controls::LensPosition,  initial_lens);
    start_controls.set(controls::FrameDurationLimits,
                       Span<const int64_t, 2>(fdur));

    ret = impl->camera->start(&start_controls);
    if (ret) {
        std::cerr << "Failed to start camera: " << ret << std::endl;
        return false;
    }

    impl->applied_lens        = initial_lens;
    impl->applied_exposure_us = initial_exposure;
    impl->applied_gain        = initial_gain;

    // Pre-fill the pipeline. Every allocated request is queued up front so
    // the camera always has frames in flight. As each one completes, the
    // capture thread re-queues it.
    for (auto &req : impl->requests) {
        impl->camera->queueRequest(req.get());
    }

    // Spawn the internal capture thread. It runs wait+bin continuously and
    // stashes the most recent Frame in the mailbox; capture() polls that
    // mailbox without blocking on libcamera.
    last_request_pop_ = std::chrono::steady_clock::now();
    impl->running.store(true);
    impl->thread = std::thread(&LibcameraCapture::capture_thread_main, this);

    std::cerr << "Camera initialized: " << stream_config.toString()
              << " (lens=" << initial_lens
              << " exposure_us=" << initial_exposure
              << " gain=" << initial_gain
              << ", " << impl->requests.size() << " requests pre-queued)" << std::endl;
    return true;
}

void LibcameraCapture::set_lens_position(float diopters) { lens_position_.store(diopters); }
void LibcameraCapture::set_exposure_us(uint32_t us)      { exposure_us_.store(us); }
void LibcameraCapture::set_gain(float multiplier)        { gain_.store(multiplier); }

namespace {
// 2x2 reduce a row range. Each thread gets its own [y0, y1) strip; reads
// come from a shared mmap'd buffer, writes go to non-overlapping ranges of
// out_buf, so no synchronisation is needed inside the parallel region.
void bin_rows(const uint8_t *bytes, unsigned int stride, int out_width,
              uint16_t *out_buf, int y0, int y1) {
    for (int y = y0; y < y1; y++) {
        int by = y * 2;
        const uint16_t *row0 = reinterpret_cast<const uint16_t *>(bytes + by       * stride);
        const uint16_t *row1 = reinterpret_cast<const uint16_t *>(bytes + (by + 1) * stride);
        uint16_t *out = out_buf + static_cast<size_t>(y) * out_width;
        int x = 0;
#if defined(__ARM_NEON)
        for (; x + 8 <= out_width; x += 8) {
            uint16x8x2_t r0 = vld2q_u16(row0 + x * 2);
            uint16x8x2_t r1 = vld2q_u16(row1 + x * 2);
            uint16x8_t s0  = vaddq_u16(r0.val[0], r0.val[1]);
            uint16x8_t s1  = vaddq_u16(r1.val[0], r1.val[1]);
            uint16x8_t avg = vshrq_n_u16(vaddq_u16(s0, s1), 2);
            vst1q_u16(out + x, avg);
        }
#endif
        for (; x < out_width; x++) {
            int bx = x * 2;
            uint32_t sum = row0[bx] + row0[bx + 1] + row1[bx] + row1[bx + 1];
            out[x] = static_cast<uint16_t>(sum >> 2);
        }
    }
}
}

void LibcameraCapture::capture_thread_main() {
    using clock = std::chrono::steady_clock;
    const StreamConfiguration &sc = impl->config->at(0);
    const int raw_width  = sc.size.width;
    const int raw_height = sc.size.height;
    const unsigned int stride = sc.stride;
    const int out_width  = raw_width  / 2;
    const int out_height = raw_height / 2;

    auto requeue = [&](Request *r) {
        r->reuse(Request::ReuseBuffers);
        auto &rc = r->controls();

        float    new_lens   = lens_position_.load();
        uint32_t new_exp_us = exposure_us_.load();
        float    new_gain   = gain_.load();

        if (new_lens != impl->applied_lens) {
            rc.set(controls::LensPosition, new_lens);
            impl->applied_lens = new_lens;
        }
        if (new_exp_us != impl->applied_exposure_us) {
            rc.set(controls::ExposureTime, static_cast<int32_t>(new_exp_us));
            const std::array<int64_t, 2> fdur = { new_exp_us, new_exp_us };
            rc.set(controls::FrameDurationLimits, Span<const int64_t, 2>(fdur));
            impl->applied_exposure_us = new_exp_us;
        }
        if (new_gain != impl->applied_gain) {
            rc.set(controls::AnalogueGain, new_gain);
            impl->applied_gain = new_gain;
        }
        impl->camera->queueRequest(r);
    };

    while (impl->running.load()) {
        Request *req = nullptr;
        {
            std::unique_lock<std::mutex> lock(impl->mtx);
            // Time-bounded wait so we re-check `running` and can shut down.
            impl->cv.wait_for(lock, std::chrono::milliseconds(100), [&]{
                return !impl->completed.empty() || !impl->running.load();
            });
            if (!impl->running.load()) return;
            if (impl->completed.empty()) continue;
            auto t_wait_ms =
                std::chrono::duration<float, std::milli>(clock::now() - last_request_pop_).count();
            last_wait_ms_.store(t_wait_ms);
            last_request_pop_ = clock::now();
            req = impl->completed.front();
            impl->completed.pop_front();
        }

        if (req->status() != Request::RequestComplete) {
            std::cerr << "Request failed (status=" << static_cast<int>(req->status()) << ")" << std::endl;
            requeue(req);
            continue;
        }

        const auto &buffers = req->buffers();
        auto it = buffers.find(impl->stream);
        if (it == buffers.end()) { requeue(req); continue; }
        FrameBuffer *buffer = it->second;
        const auto &planes = buffer->planes();
        if (planes.empty()) { requeue(req); continue; }

        void *data = mmap(nullptr, planes[0].length, PROT_READ, MAP_SHARED,
                          planes[0].fd.get(), planes[0].offset);
        if (data == MAP_FAILED) {
            std::cerr << "Failed to mmap buffer" << std::endl;
            requeue(req);
            continue;
        }

        Frame frame;
        frame.width = out_width;
        frame.height = out_height;
        frame.bit_depth = 10;
        frame.data.resize(static_cast<size_t>(out_width) * out_height);

        auto t_bin = clock::now();
        const uint8_t *bytes   = static_cast<const uint8_t *>(data);
        uint16_t      *out_buf = frame.data.data();

        // Pi 2 (and every newer Pi) has 4 cores. Spawn N-1 workers and run the
        // last strip on this thread. std::thread create+join is ~few hundred
        // us — negligible against the ms-scale work.
        unsigned hwc = std::thread::hardware_concurrency();
        int n_threads = (hwc == 0) ? 1 : std::min(4u, hwc);
        if (n_threads <= 1 || out_height < 32) {
            bin_rows(bytes, stride, out_width, out_buf, 0, out_height);
        } else {
            const int rows_per = (out_height + n_threads - 1) / n_threads;
            std::vector<std::thread> workers;
            workers.reserve(n_threads - 1);
            for (int t = 0; t < n_threads - 1; t++) {
                int y0 = t * rows_per;
                int y1 = std::min((t + 1) * rows_per, out_height);
                if (y0 >= y1) break;
                workers.emplace_back(bin_rows, bytes, stride, out_width, out_buf, y0, y1);
            }
            bin_rows(bytes, stride, out_width, out_buf,
                     (n_threads - 1) * rows_per, out_height);
            for (auto &w : workers) w.join();
        }
        last_bin_ms_.store(
            std::chrono::duration<float, std::milli>(clock::now() - t_bin).count());

        munmap(data, planes[0].length);

        // Drop the previous unread frame if the consumer hasn't picked it up;
        // we only ever keep the latest. This is the structural fix that lets
        // capture() be a fast non-blocking poll while the camera runs at its
        // own rate.
        {
            std::lock_guard<std::mutex> lock(impl->mailbox_mu);
            impl->latest = std::move(frame);
            impl->latest_seq++;
        }

        requeue(req);
    }
}

std::optional<Frame> LibcameraCapture::capture() {
    std::lock_guard<std::mutex> lock(impl->mailbox_mu);
    if (!impl->latest || impl->latest_seq == impl->consumed_seq) {
        return std::nullopt;
    }
    impl->consumed_seq = impl->latest_seq;
    Frame f = std::move(*impl->latest);
    impl->latest.reset();
    return f;
}

#endif // HAS_LIBCAMERA
