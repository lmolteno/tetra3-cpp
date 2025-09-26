#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double y, x; // pixel coordinates (y=row, x=col)
} centroid_t;

typedef struct {
    float ra, dec; // radians
    float x, y, z; // unit vector components
    float magnitude;
} star_entry_t;

typedef struct {
    uint16_t star_indices[4]; // indices into star catalog
} pattern_entry_t;

typedef struct {
    bool solved;
    float ra, dec, roll; // radians
    float fov; // horizontal FOV in radians
    float rmse; // RMS error in arcseconds
    int num_matches;
    float solve_time_ms;
    float distortion_k;
} solve_result_t;

// Opaque handle for the star solver
typedef struct star_solver* star_solver_handle_t;

// Create a new star solver with catalog data
star_solver_handle_t star_solver_create(
    const star_entry_t* stars, size_t num_stars,
    const pattern_entry_t* patterns, size_t num_patterns
);

// Destroy the star solver
void star_solver_destroy(star_solver_handle_t handle);

// Main solving function
solve_result_t star_solver_solve_from_centroids(
    star_solver_handle_t handle,
    const centroid_t* centroids, size_t num_centroids,
    int height, int width,
    double fov_estimate_deg,
    int pattern_checking_stars,
    float match_radius,
    float match_threshold,
    bool use_distortion,
    float distortion_coeff
);

// Get memory usage estimate
size_t star_solver_get_memory_usage(star_solver_handle_t handle);

// Print memory usage
void star_solver_print_memory_usage(star_solver_handle_t handle);

#ifdef __cplusplus
}
#endif