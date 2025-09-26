#pragma once

#include <esp_camera.h>
#include "star_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

// Star detection results
typedef struct {
    int num_centroids;
    centroid_t* centroids;
    float* sizes; // cluster sizes instead of magnitudes
} star_detection_result_t;

// Simple star detection from grayscale frame buffer
// Uses median subtraction, rescaling, and flood fill algorithm
star_detection_result_t* detect_stars_simple(const camera_fb_t* fb,
                                            float threshold,      // threshold after rescaling (default: 28.0)
                                            int min_area,         // minimum cluster size (default: 3)
                                            int max_stars);       // maximum stars to detect (default: 100)

// Free star detection results
void star_detection_free(star_detection_result_t* result);

#ifdef __cplusplus
}
#endif