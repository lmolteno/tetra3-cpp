#ifdef HAS_LIBCAMERA

#include "camera_capture.h"
#include <iostream>
#include <libcamera/libcamera.h>

using namespace libcamera;

struct LibcameraCapture::Impl {
    std::unique_ptr<CameraManager> cm;
    std::shared_ptr<Camera> camera;
    std::unique_ptr<CameraConfiguration> config;
    std::unique_ptr<FrameBufferAllocator> allocator;
    std::vector<std::unique_ptr<Request>> requests;
    Stream *stream = nullptr;

    // Completed request storage
    Request *completed_request = nullptr;
    bool request_completed = false;
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

    // Configure for raw Bayer (SRGGB12 for IMX477)
    impl->config = impl->camera->generateConfiguration({StreamRole::Raw});
    if (!impl->config) {
        std::cerr << "Failed to generate camera configuration" << std::endl;
        return false;
    }

    StreamConfiguration &stream_config = impl->config->at(0);
    stream_config.size.width = req_width * 2;  // Raw is full Bayer, we'll downsample
    stream_config.size.height = req_height * 2;

    CameraConfiguration::Status status = impl->config->validate();
    if (status == CameraConfiguration::Invalid) {
        std::cerr << "Camera configuration invalid" << std::endl;
        return false;
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

    // Connect signal for completed requests
    impl->camera->requestCompleted.connect(
        [this](Request *request) {
            impl->completed_request = request;
            impl->request_completed = true;
        });

    ret = impl->camera->start();
    if (ret) {
        std::cerr << "Failed to start camera: " << ret << std::endl;
        return false;
    }

    std::cout << "Camera initialized: " << stream_config.toString() << std::endl;
    return true;
}

std::optional<Frame> LibcameraCapture::capture() {
    if (!impl->camera || impl->requests.empty()) return std::nullopt;

    // Queue a request
    impl->request_completed = false;
    impl->camera->queueRequest(impl->requests[0].get());

    // Wait for completion (simple polling)
    while (!impl->request_completed) {
        impl->cm->processEvents();
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

    // Convert Bayer to grayscale by averaging 2x2 superpixels
    int raw_width = impl->config->at(0).size.width;
    int raw_height = impl->config->at(0).size.height;
    int out_width = raw_width / 2;
    int out_height = raw_height / 2;

    Frame frame;
    frame.width = out_width;
    frame.height = out_height;
    frame.bit_depth = 12;
    frame.data.resize(out_width * out_height);

    // SRGGB12 packed: each pixel is 12 bits
    // For simplicity, treat as 16-bit (most libcamera raw formats are unpacked to 16-bit)
    const uint16_t *raw = static_cast<const uint16_t *>(data);

    for (int y = 0; y < out_height; y++) {
        for (int x = 0; x < out_width; x++) {
            int by = y * 2;
            int bx = x * 2;
            // Average 2x2 Bayer superpixel
            uint32_t sum = raw[by * raw_width + bx] +
                           raw[by * raw_width + bx + 1] +
                           raw[(by + 1) * raw_width + bx] +
                           raw[(by + 1) * raw_width + bx + 1];
            frame.data[y * out_width + x] = static_cast<uint16_t>(sum / 4);
        }
    }

    munmap(data, planes[0].length);

    // Reuse the request
    impl->completed_request->reuse(Request::ReuseBuffers);

    return frame;
}

#endif // HAS_LIBCAMERA
