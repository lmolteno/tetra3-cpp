#include "star_solver.h"
#include <cmath>
#include <algorithm>
#include <chrono>
#include <functional>
#include <set>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <array>
#include <optional>
#include <memory>
#include <eigen3/Eigen/Dense>

struct star_solver {
    const star_entry_t* star_catalog_ptr;
    const pattern_entry_t* pattern_catalog_ptr;
    size_t num_stars;
    size_t num_patterns_in_catalog;

    // Database properties
    int pattern_size = 4;
    int pattern_bins = 50;
    float pattern_max_error = 0.005f;
    float max_fov_deg = 30.0f;

    static constexpr uint64_t MAGIC_RAND = 2654435761ULL;

    star_solver(const star_entry_t* stars, size_t num_stars_param,
                const pattern_entry_t* patterns, size_t num_patterns) {
        star_catalog_ptr = stars;
        pattern_catalog_ptr = patterns;
        num_stars = num_stars_param;
        num_patterns_in_catalog = num_patterns;
    }

    // Convert centroids from C arrays to C++ vectors
    std::vector<centroid_t> convert_centroids(const centroid_t* centroids, size_t num_centroids) {
        return std::vector<centroid_t>(centroids, centroids + num_centroids);
    }

    long double combinations(int n, int k) {
        if (k < 0 || k > n) {
            return 0;
        }
        if (k == 0 || k == n) {
            return 1;
        }
        if (k > n / 2) {
            k = n - k;
        }
        long double res = 1.0L;
        for (int i = 1; i <= k; ++i) {
            res = res * (n - i + 1) / i;
        }
        return res;
    }

    long double calculate_binomial_cdf(int k, int n, long double p) {
        if (k < 0) return 0.0L;
        if (k >= n) return 1.0L;
        if (p < 0.0L || p > 1.0L) return 0.0L;

        long double cdf = 0.0L;
        long double one_minus_p = 1.0L - p;

        for (int i = 0; i <= k; ++i) {
            long double term = combinations(n, i);
            term *= std::pow(p, i);
            term *= std::pow(one_minus_p, n - i);
            cdf += term;
        }
        return cdf;
    }

    std::vector<uint64_t> calculate_powers(int base, unsigned long count) {
        std::vector<uint64_t> powers(count);
        if (count > 0) {
            powers[0] = 1;
            for (int i = 1; i < count; ++i) {
                powers[i] = powers[i - 1] * base;
            }
        }
        return powers;
    }

    std::vector<uint64_t> key_to_index(const std::vector<std::vector<int>>& keys, int bin_factor, uint64_t max_index) {
        std::vector<uint64_t> hash_indices;
        if (keys.empty()) {
            return hash_indices;
        }

        size_t key_length = keys[0].size();
        if (key_length == 0) {
            return hash_indices;
        }

        std::vector<uint64_t> powers = calculate_powers(bin_factor, key_length);
        hash_indices.reserve(keys.size());

        for (const auto& key : keys) {
            uint64_t hash_val = 0;
            for (size_t i = 0; i < key_length; ++i) {
                hash_val += static_cast<uint64_t>(key[i]) * powers[i];
            }

            hash_val = (hash_val * MAGIC_RAND) % max_index;
            hash_indices.push_back(hash_val);
        }

        return hash_indices;
    }

    std::vector<centroid_t> _undistort_centroids(const std::vector<centroid_t>& centroids, int height, int width, float k_distortion) {
        std::vector<centroid_t> undistorted_centroids;
        undistorted_centroids.reserve(centroids.size());

        float img_center_y = static_cast<float>(height) / 2.0f;
        float img_center_x = static_cast<float>(width) / 2.0f;
        float half_width_f = static_cast<float>(width) / 2.0f;

        for (const auto& c_distorted : centroids) {
            float dx = c_distorted.x - img_center_x;
            float dy = c_distorted.y - img_center_y;
            float r_distorted_pixels = std::sqrt(dx * dx + dy * dy);

            float r_d = r_distorted_pixels / half_width_f;
            float r_u_scaled = r_d * (1.0f - k_distortion * r_d * r_d);
            float scale_factor = (r_distorted_pixels > 1e-9) ? (r_u_scaled * half_width_f / r_distorted_pixels) : 0.0f;

            centroid_t c_undistorted;
            c_undistorted.x = img_center_x + dx * scale_factor;
            c_undistorted.y = img_center_y + dy * scale_factor;
            undistorted_centroids.push_back(c_undistorted);
        }
        return undistorted_centroids;
    }

    double vector_distance(const std::array<double, 3>& v1, const std::array<double, 3>& v2) {
        double diff_x = v1[0] - v2[0];
        double diff_y = v1[1] - v2[1];
        double diff_z = v1[2] - v2[2];
        return std::sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    }

    double vector_angle(const std::array<double, 3>& v1, const std::array<double, 3>& v2) {
        double dist = vector_distance(v1, v2);
        return 2.0f * std::asin(dist / 2.0f);
    }

    std::vector<std::array<double, 3>> compute_vectors(const std::vector<centroid_t>& centroids, int height, int width, double fov) {
        if (width == 0) {
            return {};
        }

        double scale_factor = std::tan(fov / 2.0) / static_cast<double>(width) * 2.0;
        std::vector<std::array<double, 3>> star_vectors(centroids.size());

        double img_center_y = static_cast<double>(height) / 2.0;
        double img_center_x = static_cast<double>(width) / 2.0;

        for (size_t i = 0; i < centroids.size(); ++i) {
            star_vectors[i][0] = 1.0;
            star_vectors[i][2] = (img_center_y - static_cast<double>(centroids[i].y)) * scale_factor;
            star_vectors[i][1] = (img_center_x - static_cast<double>(centroids[i].x)) * scale_factor;

            double norm = std::sqrt(
                star_vectors[i][0] * star_vectors[i][0] +
                star_vectors[i][1] * star_vectors[i][1] +
                star_vectors[i][2] * star_vectors[i][2]
            );

            if (norm > 1e-12) {
                star_vectors[i][0] /= norm;
                star_vectors[i][1] /= norm;
                star_vectors[i][2] /= norm;
            } else {
                star_vectors[i][0] = 0.0;
                star_vectors[i][1] = 0.0;
                star_vectors[i][2] = 0.0;
            }
        }

        return star_vectors;
    }

    std::vector<double> calculate_edge_ratios(const std::vector<std::array<double, 3>>& pattern_vectors) {
        std::vector<double> edge_angles;

        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                float angle = vector_angle(pattern_vectors[i], pattern_vectors[j]);
                edge_angles.push_back(angle);
            }
        }

        std::sort(edge_angles.begin(), edge_angles.end());
        float largest_edge = edge_angles.back();

        std::vector<double> edge_ratios;
        for (int i = 0; i < 5; i++) {
            edge_ratios.push_back(edge_angles[i] / largest_edge);
        }

        return edge_ratios;
    }

    std::vector<std::vector<int>> generate_hash_code_combinations(const std::vector<double>& edge_ratios, int pattern_bins, float tolerance) {
        if (edge_ratios.size() != 5) {
            return {};
        }

        std::vector<int> hash_code_space_min(5);
        std::vector<int> hash_code_space_max(5);

        for (size_t i = 0; i < edge_ratios.size(); ++i) {
            double ratio_min_val = edge_ratios[i] - tolerance;
            double ratio_max_val = edge_ratios[i] + tolerance;

            hash_code_space_min[i] = std::max(0, static_cast<int>(ratio_min_val * pattern_bins));
            hash_code_space_max[i] = std::min(pattern_bins - 1, static_cast<int>(ratio_max_val * pattern_bins));
        }

        size_t num_dimensions = hash_code_space_min.size();

        std::vector<std::vector<int>> hash_code_ranges(num_dimensions);
        for (size_t i = 0; i < num_dimensions; ++i) {
            for (int val = hash_code_space_min[i]; val <= hash_code_space_max[i]; ++val) {
                hash_code_ranges[i].push_back(val);
            }
        }

        std::vector<std::vector<int>> all_combinations;

        std::function<void(size_t, std::vector<int>&)> generate_product = [&](size_t dimension, std::vector<int>& current_combination) {
            if (dimension == num_dimensions) {
                all_combinations.push_back(current_combination);
                return;
            }

            for (int val : hash_code_ranges[dimension]) {
                current_combination.push_back(val);
                generate_product(dimension + 1, current_combination);
                current_combination.pop_back();
            }
        };

        std::vector<int> temp;
        generate_product(0, temp);

        for (auto& combo : all_combinations) {
            std::sort(combo.begin(), combo.end());
        }

        std::set<std::vector<int>> unique_combinations(all_combinations.begin(), all_combinations.end());

        return std::vector<std::vector<int>>(unique_combinations.begin(), unique_combinations.end());
    }

    std::vector<uint64_t> get_table_indices_from_hash(uint64_t hash_index, uint64_t table_size) {
        std::vector<uint64_t> found_indices;

        for (uint64_t probe = 0; probe < table_size; probe++) {
            uint64_t probed_index = (hash_index + probe * probe) % table_size;

            if (probed_index >= num_patterns_in_catalog) break;

            const auto& test_pattern = pattern_catalog_ptr[probed_index];
            bool is_empty = true;
            for (int j = 0; j < 4; j++) {
                if (test_pattern.star_indices[j] != 0) {
                    is_empty = false;
                    break;
                }
            }

            if (is_empty) {
                break;
            } else {
                found_indices.push_back(probed_index);
            }

            if (probe > 1000) {
                std::cout << "Warning: Excessive probing, breaking at " << probe << std::endl;
                break;
            }
        }

        return found_indices;
    }

    std::vector<std::array<double, 3>> sort_pattern_by_centroid(const std::vector<std::array<double, 3>>& vectors) {
        if (vectors.empty()) {
            return {};
        }

        std::array<double, 3> centroid = {0.0, 0.0, 0.0};
        for (const auto& vec : vectors) {
            centroid[0] += vec[0];
            centroid[1] += vec[1];
            centroid[2] += vec[2];
        }
        centroid[0] /= vectors.size();
        centroid[1] /= vectors.size();
        centroid[2] /= vectors.size();

        std::vector<std::pair<double, std::array<double, 3>>> tagged_vectors;
        tagged_vectors.reserve(vectors.size());
        for (const auto& vec : vectors) {
            double dist = vector_distance(vec, centroid);
            tagged_vectors.push_back({dist, vec});
        }

        std::sort(tagged_vectors.begin(), tagged_vectors.end(),
                  [](const std::pair<double, std::array<double, 3>>& a,
                     const std::pair<double, std::array<double, 3>>& b) {
                      if (a.first != b.first) {
                          return a.first < b.first;
                      }
                      if (a.second[0] != b.second[0]) return a.second[0] < b.second[0];
                      if (a.second[1] != b.second[1]) return a.second[1] < b.second[1];
                      return a.second[2] < b.second[2];
                  });

        std::vector<std::array<double, 3>> sorted_vectors;
        sorted_vectors.reserve(vectors.size());
        for (const auto& tagged_vec : tagged_vectors) {
            sorted_vectors.push_back(tagged_vec.second);
        }

        return sorted_vectors;
    }

    std::array<std::array<double, 3>, 3> find_rotation_matrix(
        const std::vector<std::array<double, 3>>& image_vectors,
        const std::vector<std::array<double, 3>>& catalog_vectors) {
        int num_vectors = image_vectors.size();

        if (num_vectors < 3) {
            return {{{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}}};
        }

        Eigen::MatrixXf A(num_vectors, 3);
        Eigen::MatrixXf B(num_vectors, 3);

        for (int i = 0; i < num_vectors; ++i) {
            for (int j = 0; j < 3; ++j) {
                A(i, j) = image_vectors[i][j];
                B(i, j) = catalog_vectors[i][j];
            }
        }

        Eigen::Matrix3f H = A.transpose() * B;
        Eigen::BDCSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Matrix3f V = svd.matrixV();
        Eigen::Matrix3f R = (V * U.transpose()).transpose();

        if (R.determinant() < 0) {
            V.col(2) *= -1;
            R = V * U.transpose();
        }

        std::array<std::array<double, 3>, 3> rotation_matrix_std;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotation_matrix_std[i][j] = R(i, j);
            }
        }

        return rotation_matrix_std;
    }

    std::vector<int> get_nearby_stars(const std::array<float, 3>& center_vector, float radius) {
        std::vector<int> nearby_indices;
        float cos_radius = std::cos(radius);

        for (size_t i = 0; i < num_stars; i++) {
            const auto& star = star_catalog_ptr[i];
            float dot_product = center_vector[0] * star.x +
                center_vector[1] * star.y +
                center_vector[2] * star.z;

            if (dot_product > cos_radius) {
                nearby_indices.push_back(i);
            }
        }

        return nearby_indices;
    }

    std::vector<std::array<int, 2>> _find_centroid_matches(
        const std::vector<centroid_t>& image_centroids,
        const std::vector<centroid_t>& catalog_centroids,
        float match_radius_pixels) {
        std::vector<std::array<int, 2>> matched_stars;
        std::vector<bool> catalog_matched(catalog_centroids.size(), false);

        float match_radius_sq = match_radius_pixels * match_radius_pixels;

        for (size_t i = 0; i < image_centroids.size(); ++i) {
            float min_dist_sq = std::numeric_limits<float>::max();
            int best_match_idx = -1;

            for (size_t j = 0; j < catalog_centroids.size(); ++j) {
                if (catalog_matched[j])
                    continue;

                float dx = image_centroids[i].x - catalog_centroids[j].x;
                float dy = image_centroids[i].y - catalog_centroids[j].y;
                float dist_sq = dx * dx + dy * dy;

                if (dist_sq < match_radius_sq && dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                    best_match_idx = static_cast<int>(j);
                }
            }

            if (best_match_idx != -1) {
                matched_stars.push_back({static_cast<int>(i), best_match_idx});
                catalog_matched[best_match_idx] = true;
            }
        }
        return matched_stars;
    }

    std::pair<std::vector<centroid_t>, std::vector<bool>> _compute_centroids_from_vectors(
        const std::vector<std::array<float, 3>>& vectors_derot,
        int height, int width, float fov) {
        std::vector<centroid_t> centroids;
        std::vector<bool> kept_mask(vectors_derot.size(), false);

        if (width == 0 || height == 0) {
            return {centroids, kept_mask};
        }

        float tan_half_fov = std::tan(fov / 2.0f);
        float scale_factor_inv = static_cast<float>(width) / (2.0f * tan_half_fov);

        float img_center_y = static_cast<float>(height) / 2.0f;
        float img_center_x = static_cast<float>(width) / 2.0f;

        for (size_t i = 0; i < vectors_derot.size(); ++i) {
            const auto& vec = vectors_derot[i];
            if (vec[0] <= 1e-6) {
                continue;
            }

            float normalized_y = vec[2] / vec[0];
            float normalized_x = vec[1] / vec[0];

            centroid_t c{};
            c.y = img_center_y - normalized_y * scale_factor_inv;
            c.x = img_center_x - normalized_x * scale_factor_inv;

            if (c.x >= 0 && c.x < width && c.y >= 0 && c.y < height) {
                centroids.push_back(c);
                kept_mask[i] = true;
            }
        }
        return {centroids, kept_mask};
    }

    // Complete solving implementation - main solving function
    solve_result_t solve(const centroid_t* centroids, size_t num_centroids,
                        int height, int width, double fov_estimate_deg,
                        int pattern_checking_stars, float match_radius,
                        float match_threshold, bool use_distortion, float distortion_coeff) {

        auto start_time = std::chrono::high_resolution_clock::now();
        solve_result_t result = {false, 0, 0, 0, 0, 0, 0, 0, 0.0f};

        if (num_centroids < 4 || !star_catalog_ptr || !pattern_catalog_ptr || num_stars == 0 || num_patterns_in_catalog == 0) {
            return result;
        }

        std::vector<centroid_t> centroid_vec = convert_centroids(centroids, num_centroids);
        double fov = fov_estimate_deg * M_PI / 180.0f;
        float k_distortion = use_distortion ? distortion_coeff : 0.0f;

        std::vector<centroid_t> image_centroids_undist;
        if (use_distortion) {
            image_centroids_undist = _undistort_centroids(centroid_vec, height, width, k_distortion);
        } else {
            image_centroids_undist = centroid_vec;
        }

        int max_stars_to_check = std::min(static_cast<int>(num_centroids), pattern_checking_stars);

        // Try pattern matching
        for (int i = 0; i < max_stars_to_check - 3; i++) {
            for (int j = i + 1; j < max_stars_to_check - 2; j++) {
                for (int k = j + 1; k < max_stars_to_check - 1; k++) {
                    for (int l = k + 1; l < max_stars_to_check; l++) {
                        std::vector<centroid_t> pattern_centroids = {
                            centroid_vec[i], centroid_vec[j], centroid_vec[k], centroid_vec[l]
                        };
                        auto pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);
                        auto edge_ratios = calculate_edge_ratios(pattern_vectors);
                        auto hash_code_list = generate_hash_code_combinations(edge_ratios, pattern_bins, pattern_max_error);
                        auto hash_indices = key_to_index(hash_code_list, pattern_bins, num_patterns_in_catalog);

                        for (auto hash_index : hash_indices) {
                            std::vector<uint64_t> hash_match_indices = get_table_indices_from_hash(hash_index, num_patterns_in_catalog);

                            for (uint64_t catalog_index : hash_match_indices) {
                                if (catalog_index >= num_patterns_in_catalog) continue;

                                const auto& catalog_pattern = pattern_catalog_ptr[catalog_index];
                                std::vector<std::array<double, 3>> catalog_vectors;
                                bool valid_catalog_pattern = true;

                                for (int idx : catalog_pattern.star_indices) {
                                    if (idx < num_stars) {
                                        const auto& star = star_catalog_ptr[idx];
                                        catalog_vectors.push_back({star.x, star.y, star.z});
                                    } else {
                                        valid_catalog_pattern = false;
                                        break;
                                    }
                                }

                                if (!valid_catalog_pattern || catalog_vectors.size() != 4) {
                                    continue;
                                }

                                auto catalog_edge_ratios = calculate_edge_ratios(catalog_vectors);

                                bool match = true;
                                float max_error = 0.0f;
                                for (size_t m = 0; m < edge_ratios.size() && m < catalog_edge_ratios.size(); m++) {
                                    float error = std::abs(edge_ratios[m] - catalog_edge_ratios[m]);
                                    max_error = std::max(max_error, error);
                                    if (error > pattern_max_error) {
                                        match = false;
                                        break;
                                    }
                                }

                                if (match) {
                                    std::cout << "Pattern match found! Max error: " << std::fixed << std::setprecision(7) << max_error << std::endl;

                                    // Refine FOV estimate using the matched pattern
                                    float catalog_largest_edge = 0;
                                    for (int p = 0; p < 4; p++) {
                                        for (int q = p + 1; q < 4; q++) {
                                            float angle = vector_angle(catalog_vectors[p], catalog_vectors[q]);
                                            if (angle > catalog_largest_edge) {
                                                catalog_largest_edge = angle;
                                            }
                                        }
                                    }

                                    float image_largest_edge = 0;
                                    for (int p = 0; p < 4; p++) {
                                        for (int q = p + 1; q < 4; q++) {
                                            float angle = vector_angle(pattern_vectors[p], pattern_vectors[q]);
                                            if (angle > image_largest_edge) {
                                                image_largest_edge = angle;
                                            }
                                        }
                                    }

                                    if (image_largest_edge > 0.001f) {
                                        // Refine FOV based on pattern scale
                                        fov = catalog_largest_edge / image_largest_edge * fov;
                                    }

                                    std::cout << "Refined FOV: " << fov * 180.0f / M_PI << " degrees" << std::endl;

                                    // Recompute pattern vectors with refined FOV
                                    pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);

                                    // Sort both patterns by consistent criteria
                                    auto sorted_image_vectors = sort_pattern_by_centroid(pattern_vectors);
                                    auto sorted_catalog_vectors = catalog_vectors;

                                    // Calculate initial rotation matrix
                                    auto rotation_matrix_std_array = find_rotation_matrix(sorted_image_vectors, sorted_catalog_vectors);

                                    // Convert to Eigen for easier operations
                                    Eigen::Matrix3f rotation_matrix_eigen;
                                    for (int row = 0; row < 3; ++row) {
                                        for (int col = 0; col < 3; ++col) {
                                            rotation_matrix_eigen(row, col) = rotation_matrix_std_array[row][col];
                                        }
                                    }

                                    // Calculate diagonal FOV for nearby star search
                                    float fov_diagonal_rad = fov * std::sqrt(
                                        static_cast<float>(width) * width + static_cast<float>(height) * height) /
                                        static_cast<float>(width);
                                    float search_radius_rad = fov_diagonal_rad / 2.0f;

                                    // Find all star vectors inside the field of view for matching
                                    std::array<float, 3> image_center_vector_arr = {
                                        rotation_matrix_eigen(0, 0), rotation_matrix_eigen(0, 1), rotation_matrix_eigen(0, 2)
                                    };

                                    std::vector<int> nearby_star_inds = get_nearby_stars(image_center_vector_arr, search_radius_rad);
                                    std::cout << "Found " << nearby_star_inds.size() << " nearby catalog stars" << std::endl;

                                    std::vector<std::array<double, 3>> nearby_star_vectors;
                                    nearby_star_vectors.reserve(nearby_star_inds.size());
                                    for (int idx : nearby_star_inds) {
                                        if (idx < num_stars) {
                                            const auto& star = star_catalog_ptr[idx];
                                            nearby_star_vectors.push_back({star.x, star.y, star.z});
                                        }
                                    }

                                    // Derotate nearby stars to get their centroids
                                    Eigen::MatrixXf nearby_star_vectors_eigen(nearby_star_vectors.size(), 3);
                                    for (size_t r = 0; r < nearby_star_vectors.size(); ++r) {
                                        for (int c = 0; c < 3; ++c) {
                                            nearby_star_vectors_eigen(r, c) = nearby_star_vectors[r][c];
                                        }
                                    }
                                    Eigen::MatrixXf nearby_star_vectors_derot_eigen = (rotation_matrix_eigen * nearby_star_vectors_eigen.transpose()).transpose();

                                    std::vector<std::array<float, 3>> nearby_star_vectors_derot_std(nearby_star_vectors_derot_eigen.rows());
                                    for (int r = 0; r < nearby_star_vectors_derot_eigen.rows(); ++r) {
                                        for (int c = 0; c < 3; ++c) {
                                            nearby_star_vectors_derot_std[r][c] = nearby_star_vectors_derot_eigen(r, c);
                                        }
                                    }

                                    auto centroids_kept_pair = _compute_centroids_from_vectors(nearby_star_vectors_derot_std, height, width, fov);
                                    std::vector<centroid_t> nearby_star_centroids = centroids_kept_pair.first;
                                    const std::vector<bool>& kept_mask = centroids_kept_pair.second;

                                    // Filter vectors and indices using the kept mask
                                    std::vector<std::array<double, 3>> filtered_nearby_star_vectors;
                                    std::vector<int> filtered_nearby_star_inds;
                                    for (size_t idx = 0; idx < kept_mask.size(); ++idx) {
                                        if (kept_mask[idx]) {
                                            filtered_nearby_star_vectors.push_back(nearby_star_vectors[idx]);
                                            filtered_nearby_star_inds.push_back(nearby_star_inds[idx]);
                                        }
                                    }

                                    std::cout << "Kept " << nearby_star_centroids.size() << " stars in image bounds" << std::endl;

                                    // Match centroids to image stars
                                    float match_radius_pixels = static_cast<float>(width) * match_radius;
                                    std::vector<std::array<int, 2>> matched_stars = _find_centroid_matches(
                                        image_centroids_undist, nearby_star_centroids, match_radius_pixels);

                                    std::cout << "Matched " << matched_stars.size() << " stars" << std::endl;

                                    // Calculate probability of false match
                                    int num_extracted_stars = image_centroids_undist.size();
                                    int num_nearby_catalog_stars = nearby_star_centroids.size();
                                    int num_star_matches = matched_stars.size();

                                    long double prob_single_star_mismatch = static_cast<long double>(num_nearby_catalog_stars) *
                                        static_cast<long double>(match_radius) * static_cast<long double>(match_radius);

                                    int k_binom = num_extracted_stars - (num_star_matches - 2);
                                    int n_binom = num_extracted_stars;
                                    long double p_binom = 1.0L - prob_single_star_mismatch;

                                    long double prob_mismatch = calculate_binomial_cdf(k_binom, n_binom, p_binom);

                                    std::cout << "Mismatch probability = " << std::scientific << std::setprecision(2) << prob_mismatch << std::endl;

                                    if (prob_mismatch < match_threshold) {
                                        std::cout << "MATCH ACCEPTED" << std::endl;

                                        // Get vectors for all matches using refined FOV
                                        std::vector<centroid_t> matched_image_centroids;
                                        matched_image_centroids.reserve(num_star_matches);
                                        for (const auto& match_pair : matched_stars) {
                                            matched_image_centroids.push_back(image_centroids_undist[match_pair[0]]);
                                        }

                                        std::vector<std::array<double, 3>> matched_image_vectors = compute_vectors(matched_image_centroids, height, width, fov);

                                        std::vector<std::array<double, 3>> matched_catalog_vectors;
                                        matched_catalog_vectors.reserve(num_star_matches);
                                        for (const auto& match_pair : matched_stars) {
                                            matched_catalog_vectors.push_back(filtered_nearby_star_vectors[match_pair[1]]);
                                        }

                                        // Recompute rotation matrix with all matches for better accuracy
                                        rotation_matrix_std_array = find_rotation_matrix(matched_image_vectors, matched_catalog_vectors);
                                        for (int row = 0; row < 3; ++row) {
                                            for (int col = 0; col < 3; ++col) {
                                                rotation_matrix_eigen(row, col) = rotation_matrix_std_array[row][col];
                                            }
                                        }

                                        // Extract final orientation
                                        float boresight_x = rotation_matrix_eigen(0, 0);
                                        float boresight_y = rotation_matrix_eigen(0, 1);
                                        float boresight_z = rotation_matrix_eigen(0, 2);

                                        float ra = std::atan2(boresight_y, boresight_x);
                                        if (ra < 0) ra += 2.0f * M_PI;

                                        float dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));
                                        float roll = std::atan2(rotation_matrix_eigen(1, 2), rotation_matrix_eigen(2, 2));
                                        if (roll < 0) roll += 2.0f * M_PI;

                                        // Calculate RMSE
                                        std::vector<std::array<double, 3>> final_match_vectors = compute_vectors(matched_image_centroids, height, width, fov);

                                        // Rotate to sky coordinates
                                        Eigen::MatrixXf final_match_vectors_eigen(final_match_vectors.size(), 3);
                                        for (size_t r = 0; r < final_match_vectors.size(); ++r) {
                                            for (int c = 0; c < 3; ++c) {
                                                final_match_vectors_eigen(r, c) = final_match_vectors[r][c];
                                            }
                                        }
                                        Eigen::MatrixXf rotated_final_match_vectors_eigen = final_match_vectors_eigen * rotation_matrix_eigen;

                                        // Calculate residual angles
                                        float sum_angle_sq = 0.0f;
                                        for (size_t res_idx = 0; res_idx < matched_catalog_vectors.size(); ++res_idx) {
                                            std::array<double, 3> rotated_vec = {
                                                rotated_final_match_vectors_eigen(res_idx, 0),
                                                rotated_final_match_vectors_eigen(res_idx, 1),
                                                rotated_final_match_vectors_eigen(res_idx, 2)
                                            };
                                            float dist = vector_distance(rotated_vec, matched_catalog_vectors[res_idx]);
                                            float angle_rad = 2.0f * std::asin(std::min(1.0f, dist / 2.0f));
                                            sum_angle_sq += angle_rad * angle_rad;
                                        }

                                        float rmse_rad = std::sqrt(sum_angle_sq / static_cast<float>(matched_catalog_vectors.size()));
                                        float residual_arcsec = rmse_rad * 180.0f / M_PI * 3600.0f;

                                        auto end_time = std::chrono::high_resolution_clock::now();
                                        std::chrono::duration<double, std::milli> solve_time = end_time - start_time;

                                        result.solved = true;
                                        result.ra = ra;
                                        result.dec = dec;
                                        result.roll = roll;
                                        result.fov = fov;
                                        result.rmse = residual_arcsec;
                                        result.num_matches = num_star_matches;
                                        result.solve_time_ms = solve_time.count();
                                        result.distortion_k = k_distortion;
                                        return result;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> solve_time = end_time - start_time;
        result.solve_time_ms = solve_time.count();
        return result;
    }

    size_t get_memory_usage() const {
        size_t star_memory = num_stars * sizeof(star_entry_t);
        size_t pattern_memory = num_patterns_in_catalog * sizeof(pattern_entry_t);
        return star_memory + pattern_memory;
    }

    void print_memory_usage() const {
        size_t total = get_memory_usage();
        std::cout << "Memory usage:\n";
        std::cout << "  Stars: " << num_stars << " entries, " << (num_stars * sizeof(star_entry_t)) << " bytes\n";
        std::cout << "  Patterns: " << num_patterns_in_catalog << " entries, " << (num_patterns_in_catalog * sizeof(pattern_entry_t)) << " bytes\n";
        std::cout << "  Total: " << total << " bytes (" << total / 1024.0f << " KB)\n";
    }
};

extern "C" {

star_solver_handle_t star_solver_create(const star_entry_t* stars, size_t num_stars, const pattern_entry_t* patterns, size_t num_patterns) {
    star_solver* solver = new(std::nothrow) star_solver(stars, num_stars, patterns, num_patterns);
    return solver;
}

void star_solver_destroy(star_solver_handle_t handle) {
    delete static_cast<star_solver*>(handle);
}

solve_result_t star_solver_solve_from_centroids(star_solver_handle_t handle, const centroid_t* centroids, size_t num_centroids,
                                               int height, int width, double fov_estimate_deg, int pattern_checking_stars,
                                               float match_radius, float match_threshold, bool use_distortion, float distortion_coeff) {
    if (!handle) {
        return {false, 0, 0, 0, 0, 0, 0, 0, 0.0f};
    }

    auto* solver = static_cast<star_solver*>(handle);
    return solver->solve(centroids, num_centroids, height, width, fov_estimate_deg, pattern_checking_stars, match_radius, match_threshold, use_distortion, distortion_coeff);
}

size_t star_solver_get_memory_usage(star_solver_handle_t handle) {
    if (!handle) return 0;
    auto* solver = static_cast<star_solver*>(handle);
    return solver->get_memory_usage();
}

void star_solver_print_memory_usage(star_solver_handle_t handle) {
    if (!handle) return;
    auto* solver = static_cast<star_solver*>(handle);
    solver->print_memory_usage();
}

}