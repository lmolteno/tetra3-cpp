#ifdef HAS_LIBCAMERA

#include "camera_capture.h"
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <sys/mman.h>
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

    // Completed request signalled from libcamera's worker thread.
    std::mutex mtx;
    std::condition_variable cv;
    Request *completed_request = nullptr;
    bool request_completed = false;

    void on_request_completed(Request *request) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            completed_request = request;
            request_completed = true;
        }
        cv.notify_all();
    }
};

LibcameraCapture::LibcameraCapture(int width, int height)
    : req_width(width), req_height(height), impl(std::make_unique<Impl>()) {}

LibcameraCapture::~LibcameraCapture() {
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
    // (LensPosition is the source of truth). Seeded from the configured values.
    ControlList start_controls(impl->camera->controls());
    start_controls.set(controls::AeEnable, false);
    start_controls.set(controls::ExposureTime,  static_cast<int32_t>(exposure_us_.load()));
    start_controls.set(controls::AnalogueGain,  gain_.load());
    start_controls.set(controls::AfMode, controls::AfModeManual);
    start_controls.set(controls::LensPosition, lens_position_.load());

    ret = impl->camera->start(&start_controls);
    if (ret) {
        std::cerr << "Failed to start camera: " << ret << std::endl;
        return false;
    }

    std::cout << "Camera initialized: " << stream_config.toString()
              << " (lens=" << lens_position_.load()
              << " exposure_us=" << exposure_us_.load()
              << " gain=" << gain_.load() << ")" << std::endl;
    return true;
}

void LibcameraCapture::set_lens_position(float diopters) { lens_position_.store(diopters); }
void LibcameraCapture::set_exposure_us(uint32_t us)      { exposure_us_.store(us); }
void LibcameraCapture::set_gain(float multiplier)        { gain_.store(multiplier); }

std::optional<Frame> LibcameraCapture::capture() {
    if (!impl->camera || impl->requests.empty()) return std::nullopt;

    // Apply the latest manual settings on this request so web UI tweaks take
    // effect on the next exposure.
    auto &rc = impl->requests[0]->controls();
    rc.set(controls::LensPosition, lens_position_.load());
    rc.set(controls::ExposureTime, static_cast<int32_t>(exposure_us_.load()));
    rc.set(controls::AnalogueGain, gain_.load());

    {
        std::lock_guard<std::mutex> lock(impl->mtx);
        impl->request_completed = false;
        impl->completed_request = nullptr;
    }
    impl->camera->queueRequest(impl->requests[0].get());

    // Wait for the request_completed signal callback to fire.
    {
        std::unique_lock<std::mutex> lock(impl->mtx);
        impl->cv.wait(lock, [&]{ return impl->request_completed; });
    }

    if (impl->completed_request->status() != Request::RequestComplete) {
        std::cerr << "Request failed" << std::endl;
        return std::nullopt;
    }

    // Get the buffer
    const auto &buffers = impl->completed_request->buffers();
    auto it = buffers.find(impl->stream);
    if (it == buffers.end()) return std::nullopt;

    FrameBuffer *buffer = it->second;
    const auto &planes = buffer->planes();
    if (planes.empty()) return std::nullopt;

    // Map the buffer
    void *data = mmap(nullptr, planes[0].length, PROT_READ, MAP_SHARED,
                      planes[0].fd.get(), planes[0].offset);
    if (data == MAP_FAILED) {
        std::cerr << "Failed to mmap buffer" << std::endl;
        return std::nullopt;
    }

    // 2x2-bin the Bayer mosaic into grayscale. With SBGGR10 each pixel is one
    // uint16; libcamera reports the row stride in bytes (may be padded), so we
    // index rows by stride and reinterpret each row as uint16_t*.
    const StreamConfiguration &sc = impl->config->at(0);
    int raw_width  = sc.size.width;
    int raw_height = sc.size.height;
    unsigned int stride = sc.stride;
    int out_width  = raw_width  / 2;
    int out_height = raw_height / 2;

    Frame frame;
    frame.width = out_width;
    frame.height = out_height;
    frame.bit_depth = 10;
    frame.data.resize(static_cast<size_t>(out_width) * out_height);

    const uint8_t *bytes = static_cast<const uint8_t *>(data);
    for (int y = 0; y < out_height; y++) {
        int by = y * 2;
        const uint16_t *row0 = reinterpret_cast<const uint16_t *>(bytes + by       * stride);
        const uint16_t *row1 = reinterpret_cast<const uint16_t *>(bytes + (by + 1) * stride);
        uint16_t *out = frame.data.data() + static_cast<size_t>(y) * out_width;
        for (int x = 0; x < out_width; x++) {
            int bx = x * 2;
            uint32_t sum = row0[bx] + row0[bx + 1] + row1[bx] + row1[bx + 1];
            out[x] = static_cast<uint16_t>(sum >> 2);
        }
    }

    munmap(data, planes[0].length);

    // Reuse the request
    impl->completed_request->reuse(Request::ReuseBuffers);

    return frame;
}

#endif // HAS_LIBCAMERA
