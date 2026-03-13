#pragma once

#include <tetra3/types.h>
#include "frame_source.h"
#include <vector>
#include <cstdint>

struct DetectedStar {
    Centroid centroid; // sub-pixel (y, x)
    double snr, peak_value, flux;
    int cluster_size;
};

struct StarDetectorConfig {
    int sigma_clip_iterations = 3;
    float sigma_clip_factor = 3.0f;
    float detection_sigma = 5.0f;
    int min_cluster_size = 3;
    int max_cluster_size = 500;
    int max_stars = 200;
    bool verbose = false;
};

class StarDetector {
public:
    explicit StarDetector(const StarDetectorConfig &config = {});

    std::vector<DetectedStar> detect(const Frame &frame);

private:
    StarDetectorConfig config;

    struct BackgroundStats {
        double mean, stddev;
    };

    BackgroundStats estimate_background(const uint16_t *data, int width, int height);

    struct Point { int y, x; };

    int flood_fill(const std::vector<uint8_t> &mask, int width, int height,
                   int start_y, int start_x, std::vector<uint8_t> &visited,
                   std::vector<Point> &cluster_points);
};
