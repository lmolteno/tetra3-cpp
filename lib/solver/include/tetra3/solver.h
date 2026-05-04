#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <functional>
#include <set>
#include <iostream>
#include <iomanip>
#include <optional>
#include <limits>
#include <Eigen/Dense>
#include "types.h"

class SimpleStarSolver {
private:
    std::vector<StarEntry> star_catalog;
    std::vector<PatternEntry> pattern_catalog;

    int pattern_size = 4;
    int pattern_bins = 50;
    float pattern_max_error = 0.005f;
    float max_fov_deg = 30.0f;
    size_t num_patterns_in_catalog = 0;

    // Lens calibration: optical-axis offset from geometric center, plus a
    // radial-only distortion polynomial.
    //   F(r²) = sum_i k[i] * r^(2*(i+1))
    //   undistort: rx_pred = rx_obs * (1 - F(r_obs²))
    // r is normalized by half-width. When lens_k is empty, the solver
    // behaves as a plain pinhole.
    float lens_cx_offset = 0.0f;
    float lens_cy_offset = 0.0f;
    std::vector<float> lens_k;

    static constexpr uint64_t MAGIC_RAND = 2654435761ULL;

    std::vector<uint64_t> calculate_powers(int base, unsigned long count);
    std::vector<uint64_t> key_to_index(const std::vector<std::vector<int>> &keys, int bin_factor, uint64_t max_index);

    // Single-k undistortion (legacy path used by the solver's per-frame fit).
    std::vector<Centroid> _undistort_centroids(
        const std::vector<Centroid> &centroids,
        int height, int width, float k_distortion);

    // Polynomial undistortion using the lens calibration set on the solver.
    std::vector<Centroid> _undistort_centroids_lens(
        const std::vector<Centroid> &centroids,
        int height, int width);

    std::vector<std::array<double, 3>> sort_pattern_by_centroid(const std::vector<std::array<double, 3>> &vectors);

    std::vector<std::array<double, 3>> compute_vectors(
        const std::vector<Centroid> &centroids,
        int height, int width, double fov);

    double vector_distance(const std::array<double, 3> &v1, const std::array<double, 3> &v2);
    double vector_angle(const std::array<double, 3> &v1, const std::array<double, 3> &v2);

    std::vector<double> calculate_edge_ratios(const std::vector<std::array<double, 3>> &pattern_vectors);

    std::vector<std::vector<int>> generate_hash_code_combinations(
        const std::vector<double> &edge_ratios,
        int pattern_bins,
        float tolerance);

    std::array<std::array<double, 3>, 3> find_rotation_matrix(
        const std::vector<std::array<double, 3>> &image_vectors,
        const std::vector<std::array<double, 3>> &catalog_vectors);

    std::vector<int> get_nearby_stars(const std::array<float, 3> &center_vector, float radius);

    std::vector<std::array<int, 2>> _find_centroid_matches(
        const std::vector<Centroid> &image_centroids,
        const std::vector<Centroid> &catalog_centroids,
        float match_radius_pixels);

    std::pair<std::vector<Centroid>, std::vector<bool>> _compute_centroids_from_vectors(
        const std::vector<std::array<float, 3>> &vectors_derot,
        int height, int width, float fov);

public:
    void load_star_catalog(const std::vector<StarEntry> &stars);
    void load_pattern_catalog(const std::vector<PatternEntry> &patterns);

    // Configure a fixed lens calibration: optical-axis offset from the
    // geometric image center (pixels) and a radial polynomial in r²
    // (with r = pixel_distance / half_width). When set, the solver uses
    // this throughout: undistortion of input centroids, projection of
    // catalog stars to image, and matching. Pass an empty k vector and
    // zero offsets to disable.
    void set_lens_calibration(float cx_offset, float cy_offset,
                              const std::vector<float> &k_radial);
    bool has_lens_calibration() const { return !lens_k.empty(); }

    std::vector<uint64_t> get_table_indices_from_hash(uint64_t hash_index, uint64_t table_size);

    SolveResult solve_from_centroids(
        const std::vector<Centroid> &centroids,
        int height, int width,
        double fov_estimate_deg = 20.0f,
        int pattern_checking_stars = 8,
        float match_radius = 0.01,
        float match_threshold = 0.001,
        std::optional<float> distortion_coeff_in = std::nullopt,
        float fov_max_error_deg = 0.0f,
        std::vector<MatchInfo> *matches_out = nullptr,
        bool freeze_distortion = false,
        float extended_match_radius = 0.0f);

    size_t get_memory_usage() const;
    void print_memory_usage() const;
};
