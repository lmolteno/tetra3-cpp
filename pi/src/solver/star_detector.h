#pragma once

#include <tetra3/types.h>
#include "capture/frame_source.h"
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
    int max_cluster_size = 30;
    int max_stars = 200;
    int bg_tile_size = 64;   // tile size for local background estimation
    // Reject low-contrast clusters: ratio of (peak - 7x7 ring mean) to
    // (ring mean - local background). Real stars are sharply peaked
    // (~6+ on twilight scenes); diffuse structures and cloud edges score <1.
    float min_sharpness = 2.0f;
    bool verbose = false;
};

// Per-call stage breakdown. Filled in by detect() when out is non-null.
struct StarDetectorStats {
    double bg_ms = 0;       // tiled sigma-clipped background estimation
    double mask_ms = 0;     // bilinear-interpolated threshold mask
    double cluster_ms = 0;  // flood fill + centroiding + SNR
};

class StarDetector {
public:
    explicit StarDetector(const StarDetectorConfig &config = {});

    std::vector<DetectedStar> detect(const Frame &frame,
                                     StarDetectorStats *stats = nullptr);

private:
    StarDetectorConfig config;

    struct BackgroundStats {
        double mean, stddev;
    };

    BackgroundStats estimate_background(const uint16_t *data, int width, int height);

    struct TiledBackground {
        std::vector<float> tile_mean;
        std::vector<float> tile_std;
        int tiles_x, tile_size;

        float mean_at(int y, int x) const {
            return tile_mean[(y / tile_size) * tiles_x + x / tile_size];
        }
        float std_at(int y, int x) const {
            return tile_std[(y / tile_size) * tiles_x + x / tile_size];
        }
    };

    TiledBackground compute_local_background(const uint16_t *data, int width, int height);

    struct Point { int y, x; };

    int flood_fill(const std::vector<uint8_t> &mask, int width, int height,
                   int start_y, int start_x, std::vector<uint8_t> &visited,
                   std::vector<Point> &cluster_points);
};
