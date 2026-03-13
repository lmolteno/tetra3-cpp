#include "image_loader.h"
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

ImageFileSource::ImageFileSource(const std::string &path, bool verbose)
    : file_path(path), verbose(verbose) {}

bool ImageFileSource::initialize() {
    cv::Mat img = cv::imread(file_path, cv::IMREAD_UNCHANGED | cv::IMREAD_ANYDEPTH);
    if (img.empty()) {
        std::cerr << "Failed to load image: " << file_path << std::endl;
        return false;
    }

    // Convert to grayscale if needed
    cv::Mat gray;
    if (img.channels() > 1) {
        cv::cvtColor(img, gray, cv::COLOR_BGR2GRAY);
    } else {
        gray = img;
    }

    // Determine bit depth and convert to 16-bit
    cv::Mat gray16;
    if (gray.depth() == CV_16U) {
        gray16 = gray;
        frame.bit_depth = 16;
    } else if (gray.depth() == CV_8U) {
        gray.convertTo(gray16, CV_16U, 257.0); // scale 0-255 to 0-65535
        frame.bit_depth = 8;
    } else if (gray.depth() == CV_32F || gray.depth() == CV_64F) {
        double minVal, maxVal;
        cv::minMaxLoc(gray, &minVal, &maxVal);
        gray.convertTo(gray16, CV_16U, 65535.0 / maxVal);
        frame.bit_depth = 16;
    } else {
        gray.convertTo(gray16, CV_16U);
        frame.bit_depth = 16;
    }

    frame.width = gray16.cols;
    frame.height = gray16.rows;
    frame.data.resize(frame.width * frame.height);

    // Copy to contiguous buffer
    for (int y = 0; y < frame.height; y++) {
        const uint16_t *row = gray16.ptr<uint16_t>(y);
        std::copy(row, row + frame.width, frame.data.data() + y * frame.width);
    }

    loaded = true;
    if (verbose) {
        std::cout << "Loaded image: " << file_path << " (" << frame.width << "x" << frame.height
                  << ", " << frame.bit_depth << "-bit)" << std::endl;
    }
    return true;
}

std::optional<Frame> ImageFileSource::capture() {
    if (!loaded) return std::nullopt;
    return frame; // Return the same frame each time
}
