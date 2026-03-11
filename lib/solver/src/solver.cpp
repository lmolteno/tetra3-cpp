#include <tetra3/solver.h>

std::vector<uint64_t> SimpleStarSolver::calculate_powers(int base, unsigned long count) {
    std::vector<uint64_t> powers(count);
    if (count > 0) {
        powers[0] = 1;
        for (unsigned long i = 1; i < count; ++i) {
            powers[i] = powers[i - 1] * base;
        }
    }
    return powers;
}

std::vector<uint64_t> SimpleStarSolver::key_to_index(const std::vector<std::vector<int>> &keys, int bin_factor,
                                                      uint64_t max_index) {
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

    for (const auto &key: keys) {
        uint64_t hash_val = 0;
        for (size_t i = 0; i < key_length; ++i) {
            hash_val += static_cast<uint64_t>(key[i]) * powers[i];
        }
        hash_val = (hash_val * MAGIC_RAND) % max_index;
        hash_indices.push_back(hash_val);
    }

    return hash_indices;
}

std::vector<Centroid> SimpleStarSolver::_undistort_centroids(
    const std::vector<Centroid> &centroids,
    int height, int width, float k_distortion) {
    std::vector<Centroid> undistorted_centroids;
    undistorted_centroids.reserve(centroids.size());

    float img_center_y = static_cast<float>(height) / 2.0f;
    float img_center_x = static_cast<float>(width) / 2.0f;
    float half_width_f = static_cast<float>(width) / 2.0f;

    for (const auto &c_distorted: centroids) {
        float dx = c_distorted.x - img_center_x;
        float dy = c_distorted.y - img_center_y;
        float r_distorted_pixels = std::sqrt(dx * dx + dy * dy);
        float r_d = r_distorted_pixels / half_width_f;
        float r_u_scaled = r_d * (1.0f - k_distortion * r_d * r_d);
        float scale_factor = (r_distorted_pixels > 1e-9) ? (r_u_scaled * half_width_f / r_distorted_pixels) : 0.0f;

        Centroid c_undistorted;
        c_undistorted.x = img_center_x + dx * scale_factor;
        c_undistorted.y = img_center_y + dy * scale_factor;
        undistorted_centroids.push_back(c_undistorted);
    }
    return undistorted_centroids;
}

std::vector<std::array<double, 3>> SimpleStarSolver::sort_pattern_by_centroid(const std::vector<std::array<double, 3>> &vectors) {
    if (vectors.empty()) {
        return {};
    }

    std::array<double, 3> centroid = {0.0, 0.0, 0.0};
    for (const auto &vec: vectors) {
        centroid[0] += vec[0];
        centroid[1] += vec[1];
        centroid[2] += vec[2];
    }
    centroid[0] /= vectors.size();
    centroid[1] /= vectors.size();
    centroid[2] /= vectors.size();

    std::vector<std::pair<double, std::array<double, 3>>> tagged_vectors;
    tagged_vectors.reserve(vectors.size());
    for (const auto &vec: vectors) {
        double dist = vector_distance(vec, centroid);
        tagged_vectors.push_back({dist, vec});
    }

    std::sort(tagged_vectors.begin(), tagged_vectors.end(),
              [](const std::pair<double, std::array<double, 3>> &a,
                 const std::pair<double, std::array<double, 3>> &b) {
                  if (a.first != b.first) {
                      return a.first < b.first;
                  }
                  if (a.second[0] != b.second[0]) return a.second[0] < b.second[0];
                  if (a.second[1] != b.second[1]) return a.second[1] < b.second[1];
                  return a.second[2] < b.second[2];
              });

    std::vector<std::array<double, 3>> sorted_vectors;
    sorted_vectors.reserve(vectors.size());
    for (const auto &tagged_vec: tagged_vectors) {
        sorted_vectors.push_back(tagged_vec.second);
    }

    return sorted_vectors;
}

std::vector<std::array<double, 3>> SimpleStarSolver::compute_vectors(
    const std::vector<Centroid> &centroids,
    int height, int width, double fov) {
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

double SimpleStarSolver::vector_distance(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    double diff_x = v1[0] - v2[0];
    double diff_y = v1[1] - v2[1];
    double diff_z = v1[2] - v2[2];
    return std::sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
}

double SimpleStarSolver::vector_angle(const std::array<double, 3> &v1, const std::array<double, 3> &v2) {
    double dist = vector_distance(v1, v2);
    return 2.0 * std::asin(dist / 2.0);
}

std::vector<double> SimpleStarSolver::calculate_edge_ratios(const std::vector<std::array<double, 3>> &pattern_vectors) {
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

std::vector<std::vector<int>> SimpleStarSolver::generate_hash_code_combinations(
    const std::vector<double> &edge_ratios,
    int pattern_bins,
    float tolerance) {
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

    std::function<void(size_t, std::vector<int> &)> generate_product =
            [&](size_t dimension, std::vector<int> &current_combination) {
        if (dimension == num_dimensions) {
            all_combinations.push_back(current_combination);
            return;
        }

        for (int val: hash_code_ranges[dimension]) {
            current_combination.push_back(val);
            generate_product(dimension + 1, current_combination);
            current_combination.pop_back();
        }
    };

    std::vector<int> temp;
    generate_product(0, temp);

    for (auto &combo: all_combinations) {
        std::sort(combo.begin(), combo.end());
    }

    std::set<std::vector<int>> unique_combinations(all_combinations.begin(), all_combinations.end());

    return std::vector<std::vector<int>>(unique_combinations.begin(), unique_combinations.end());
}

std::array<std::array<double, 3>, 3> SimpleStarSolver::find_rotation_matrix(
    const std::vector<std::array<double, 3>> &image_vectors,
    const std::vector<std::array<double, 3>> &catalog_vectors) {
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

std::vector<int> SimpleStarSolver::get_nearby_stars(const std::array<float, 3> &center_vector, float radius) {
    std::vector<int> nearby_indices;
    float cos_radius = std::cos(radius);

    for (size_t i = 0; i < star_catalog.size(); i++) {
        const auto &star = star_catalog[i];
        float dot_product = center_vector[0] * star.x +
                            center_vector[1] * star.y +
                            center_vector[2] * star.z;

        if (dot_product > cos_radius) {
            nearby_indices.push_back(i);
        }
    }

    return nearby_indices;
}

void SimpleStarSolver::load_star_catalog(const std::vector<StarEntry> &stars) {
    star_catalog = stars;
}

void SimpleStarSolver::load_pattern_catalog(const std::vector<PatternEntry> &patterns) {
    pattern_catalog = patterns;
    num_patterns_in_catalog = patterns.size();
}

std::vector<uint64_t> SimpleStarSolver::get_table_indices_from_hash(uint64_t hash_index, uint64_t table_size) {
    std::vector<uint64_t> found_indices;

    for (uint64_t probe = 0; probe < table_size; probe++) {
        uint64_t probed_index = (hash_index + probe * probe) % table_size;

        const auto &test_pattern = pattern_catalog[probed_index];
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

std::vector<std::array<int, 2>> SimpleStarSolver::_find_centroid_matches(
    const std::vector<Centroid> &image_centroids,
    const std::vector<Centroid> &catalog_centroids,
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

std::pair<std::vector<Centroid>, std::vector<bool>> SimpleStarSolver::_compute_centroids_from_vectors(
    const std::vector<std::array<float, 3>> &vectors_derot,
    int height, int width, float fov) {
    std::vector<Centroid> centroids;
    std::vector<bool> kept_mask(vectors_derot.size(), false);

    if (width == 0 || height == 0) {
        return {centroids, kept_mask};
    }

    float tan_half_fov = std::tan(fov / 2.0f);
    float scale_factor_inv = static_cast<float>(width) / (2.0f * tan_half_fov);

    float img_center_y = static_cast<float>(height) / 2.0f;
    float img_center_x = static_cast<float>(width) / 2.0f;

    for (size_t i = 0; i < vectors_derot.size(); ++i) {
        const auto &vec = vectors_derot[i];
        if (vec[0] <= 1e-6) {
            continue;
        }

        float normalized_y = vec[2] / vec[0];
        float normalized_x = vec[1] / vec[0];

        Centroid c{};
        c.y = img_center_y - normalized_y * scale_factor_inv;
        c.x = img_center_x - normalized_x * scale_factor_inv;

        if (c.x >= 0 && c.x < width && c.y >= 0 && c.y < height) {
            centroids.push_back(c);
            kept_mask[i] = true;
        }
    }
    return {centroids, kept_mask};
}

SolveResult SimpleStarSolver::solve_from_centroids(
    const std::vector<Centroid> &centroids,
    int height, int width,
    double fov_estimate_deg,
    int pattern_checking_stars,
    float match_radius,
    float match_threshold,
    std::optional<float> distortion_coeff_in,
    float fov_max_error_deg) {
    auto start_time = std::chrono::high_resolution_clock::now();
    SolveResult result = {false, 0, 0, 0, 0, 0, 0, 0, 0.0f};

    if (centroids.size() < 4) {
        return result;
    }

    double fov = fov_estimate_deg * M_PI / 180.0;
    double fov_max_error_rad = fov_max_error_deg > 0
        ? fov_max_error_deg * M_PI / 180.0
        : 0.0;
    double k_distortion = distortion_coeff_in.value_or(0.0f);

    std::vector<Centroid> image_centroids_undist;
    if (distortion_coeff_in.has_value()) {
        image_centroids_undist = _undistort_centroids(centroids, height, width, k_distortion);
    } else {
        image_centroids_undist = centroids;
    }

    int num_extracted_stars_raw = centroids.size();
    int max_stars_to_check = std::min(num_extracted_stars_raw, pattern_checking_stars);

    for (int i = 0; i < max_stars_to_check - 3; i++) {
        for (int j = i + 1; j < max_stars_to_check - 2; j++) {
            for (int k = j + 1; k < max_stars_to_check - 1; k++) {
                for (int l = k + 1; l < max_stars_to_check; l++) {
                    std::vector<Centroid> pattern_centroids = {
                        centroids[i], centroids[j], centroids[k], centroids[l]
                    };
                    auto pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);
                    auto edge_ratios = calculate_edge_ratios(pattern_vectors);
                    auto hash_code_list = generate_hash_code_combinations(
                        edge_ratios, pattern_bins, pattern_max_error);

                    auto hash_indices = key_to_index(hash_code_list, pattern_bins, pattern_catalog.size());

                    for (auto hash_index: hash_indices) {
                        std::vector<uint64_t> hash_match_indices = get_table_indices_from_hash(
                            hash_index, pattern_catalog.size());

                        if (hash_match_indices.empty()) {
                            continue;
                        }

                        for (uint64_t catalog_index: hash_match_indices) {
                            const auto &catalog_pattern = pattern_catalog[catalog_index];

                            std::vector<std::array<double, 3>> catalog_vectors;
                            bool valid_catalog_pattern = true;

                            for (int idx: catalog_pattern.star_indices) {
                                if (static_cast<size_t>(idx) < star_catalog.size()) {
                                    const auto &star = star_catalog[idx];
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
                                std::cout << "Pattern match found! Max error: " << std::fixed <<
                                        std::setprecision(7) << max_error << std::endl;

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
                                    double new_fov = catalog_largest_edge / image_largest_edge * fov;
                                    // Reject if rescaled FOV is unreasonably far from estimate
                                    if (fov_max_error_rad > 0 &&
                                        std::abs(new_fov - fov_estimate_deg * M_PI / 180.0) > fov_max_error_rad) {
                                        continue;
                                    }
                                    fov = new_fov;
                                }

                                pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);

                                auto sorted_image_vectors = sort_pattern_by_centroid(pattern_vectors);
                                auto sorted_catalog_vectors = catalog_vectors;

                                auto rotation_matrix_std_array = find_rotation_matrix(
                                    sorted_image_vectors, sorted_catalog_vectors);

                                Eigen::Matrix3f rotation_matrix_eigen;
                                for (int row = 0; row < 3; ++row) {
                                    for (int col = 0; col < 3; ++col) {
                                        rotation_matrix_eigen(row, col) = rotation_matrix_std_array[row][col];
                                    }
                                }

                                float fov_diagonal_rad = fov * std::sqrt(
                                                             static_cast<float>(width) * width + static_cast<float>(
                                                                 height) * height) /
                                                         static_cast<float>(width);
                                float search_radius_rad = fov_diagonal_rad / 2.0f;

                                std::array<float, 3> image_center_vector_arr = {
                                    rotation_matrix_eigen(0, 0), rotation_matrix_eigen(0, 1),
                                    rotation_matrix_eigen(0, 2)
                                };

                                std::vector<int> nearby_star_inds = get_nearby_stars(
                                    image_center_vector_arr, search_radius_rad);

                                std::vector<std::array<double, 3>> nearby_star_vectors;
                                nearby_star_vectors.reserve(nearby_star_inds.size());
                                for (int idx: nearby_star_inds) {
                                    if (static_cast<size_t>(idx) < star_catalog.size()) {
                                        const auto &star = star_catalog[idx];
                                        nearby_star_vectors.push_back({star.x, star.y, star.z});
                                    }
                                }

                                Eigen::MatrixXf nearby_star_vectors_eigen(nearby_star_vectors.size(), 3);
                                for (size_t r = 0; r < nearby_star_vectors.size(); ++r) {
                                    for (int c = 0; c < 3; ++c) {
                                        nearby_star_vectors_eigen(r, c) = nearby_star_vectors[r][c];
                                    }
                                }
                                Eigen::MatrixXf nearby_star_vectors_derot_eigen = (
                                    rotation_matrix_eigen * nearby_star_vectors_eigen.transpose()).transpose();

                                std::vector<std::array<float, 3>> nearby_star_vectors_derot_std(
                                    nearby_star_vectors_derot_eigen.rows());
                                for (int r = 0; r < nearby_star_vectors_derot_eigen.rows(); ++r) {
                                    for (int c = 0; c < 3; ++c) {
                                        nearby_star_vectors_derot_std[r][c] = nearby_star_vectors_derot_eigen(r, c);
                                    }
                                }

                                auto centroids_kept_pair = _compute_centroids_from_vectors(
                                    nearby_star_vectors_derot_std, height, width, fov);
                                std::vector<Centroid> nearby_star_centroids = centroids_kept_pair.first;
                                const std::vector<bool> &kept_mask = centroids_kept_pair.second;

                                std::vector<std::array<double, 3>> filtered_nearby_star_vectors;
                                std::vector<int> filtered_nearby_star_inds;
                                for (size_t idx = 0; idx < kept_mask.size(); ++idx) {
                                    if (kept_mask[idx]) {
                                        filtered_nearby_star_vectors.push_back(nearby_star_vectors[idx]);
                                        filtered_nearby_star_inds.push_back(nearby_star_inds[idx]);
                                    }
                                }
                                nearby_star_vectors = filtered_nearby_star_vectors;
                                nearby_star_inds = filtered_nearby_star_inds;

                                int num_image_centroids = image_centroids_undist.size();
                                if (static_cast<int>(nearby_star_centroids.size()) > num_image_centroids) {
                                    nearby_star_centroids.resize(num_image_centroids);
                                    nearby_star_vectors.resize(num_image_centroids);
                                    nearby_star_inds.resize(num_image_centroids);
                                }

                                float match_radius_pixels = static_cast<float>(width) * match_radius;
                                std::vector<std::array<int, 2>> matched_stars = _find_centroid_matches(
                                    image_centroids_undist, nearby_star_centroids, match_radius_pixels);

                                float boresight_x = rotation_matrix_eigen(0, 0);
                                float boresight_y = rotation_matrix_eigen(1, 0);
                                float boresight_z = rotation_matrix_eigen(2, 0);

                                float ra = std::atan2(boresight_y, boresight_x);
                                if (ra < 0) ra += 2.0f * M_PI;

                                float dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));

                                float roll = std::atan2(rotation_matrix_eigen(1, 2),
                                                        rotation_matrix_eigen(2, 2));
                                if (roll < 0) roll += 2.0f * M_PI;

                                int num_extracted_stars = image_centroids_undist.size();
                                int num_nearby_catalog_stars = nearby_star_centroids.size();
                                int num_star_matches = matched_stars.size();

                                std::cout << "Number of nearby stars: " << num_nearby_catalog_stars
                                        << ", total matched: " << num_star_matches << std::endl;

                                long double prob_single_star_mismatch = static_cast<long double>(
                                                                            num_nearby_catalog_stars) * static_cast<
                                                                            long double>(match_radius) *
                                                                        static_cast<long double>(match_radius);

                                int k_binom = num_extracted_stars - (num_star_matches - 2);
                                int n_binom = num_extracted_stars;
                                long double p_binom = 1.0L - prob_single_star_mismatch;

                                long double prob_mismatch = calculate_binomial_cdf(k_binom, n_binom, p_binom);

                                std::cout << "Mismatch probability = " << std::scientific << std::setprecision(2) <<
                                        prob_mismatch
                                        << ", at FOV = " << std::fixed << std::setprecision(5) << fov * 180.0 /
                                        M_PI <<
                                        "deg" << std::endl;

                                if (prob_mismatch < match_threshold) {
                                    std::cout << "MATCH ACCEPTED" << std::endl;

                                    std::vector<Centroid> matched_image_centroids;
                                    matched_image_centroids.reserve(num_star_matches);
                                    for (const auto &match_pair: matched_stars) {
                                        matched_image_centroids.push_back(image_centroids_undist[match_pair[0]]);
                                    }

                                    std::vector<std::array<double, 3>> matched_image_vectors = compute_vectors(
                                        matched_image_centroids, height, width, fov);

                                    std::vector<std::array<double, 3>> matched_catalog_vectors;
                                    matched_catalog_vectors.reserve(num_star_matches);
                                    for (const auto &match_pair: matched_stars) {
                                        matched_catalog_vectors.push_back(nearby_star_vectors[match_pair[1]]);
                                    }

                                    rotation_matrix_std_array = find_rotation_matrix(
                                        matched_image_vectors, matched_catalog_vectors);
                                    for (int row = 0; row < 3; ++row) {
                                        for (int col = 0; col < 3; ++col) {
                                            rotation_matrix_eigen(row, col) = rotation_matrix_std_array[row][col];
                                        }
                                    }

                                    boresight_x = rotation_matrix_eigen(0, 0);
                                    boresight_y = rotation_matrix_eigen(0, 1);
                                    boresight_z = rotation_matrix_eigen(0, 2);

                                    ra = std::atan2(boresight_y, boresight_x);
                                    if (ra < 0) ra += 2.0f * M_PI;

                                    dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));

                                    roll = std::atan2(rotation_matrix_eigen(1, 2),
                                                      rotation_matrix_eigen(2, 2));
                                    if (roll < 0) roll += 2.0f * M_PI;

                                    if (!distortion_coeff_in.has_value()) {
                                        std::vector<float> angles_camera;
                                        std::vector<float> angles_catalogue;

                                        for (size_t m1 = 0; m1 < matched_image_vectors.size(); ++m1) {
                                            for (size_t m2 = m1 + 1; m2 < matched_image_vectors.size(); ++m2) {
                                                angles_camera.push_back(vector_angle(
                                                    matched_image_vectors[m1], matched_image_vectors[m2]));
                                                angles_catalogue.push_back(vector_angle(
                                                    matched_catalog_vectors[m1], matched_catalog_vectors[m2]));
                                            }
                                        }

                                        if (!angles_camera.empty()) {
                                            float mean_ratio = 0.0f;
                                            int count_ratios = 0;
                                            for (size_t ratio_idx = 0; ratio_idx < angles_camera.size(); ++ratio_idx) {
                                                if (angles_camera[ratio_idx] > 1e-9) {
                                                    mean_ratio += angles_catalogue[ratio_idx] / angles_camera[ratio_idx];
                                                    count_ratios++;
                                                }
                                            }
                                            if (count_ratios > 0) {
                                                fov *= (mean_ratio / count_ratios);
                                            }
                                        }
                                        k_distortion = 0.0f;
                                        image_centroids_undist = matched_image_centroids;
                                    } else {
                                        Eigen::MatrixXf matched_catalog_vectors_eigen(
                                            matched_catalog_vectors.size(), 3);
                                        for (size_t r = 0; r < matched_catalog_vectors.size(); ++r) {
                                            for (int c = 0; c < 3; ++c) {
                                                matched_catalog_vectors_eigen(r, c) = matched_catalog_vectors[r][c];
                                            }
                                        }
                                        Eigen::MatrixXf matched_catalog_vectors_derot_eigen =
                                                matched_catalog_vectors_eigen * rotation_matrix_eigen.transpose();

                                        Eigen::VectorXf tangent_matched_catalog_vectors(
                                            matched_catalog_vectors_derot_eigen.rows());
                                        for (int r = 0; r < matched_catalog_vectors_derot_eigen.rows(); ++r) {
                                            float norm_yz = std::sqrt(
                                                matched_catalog_vectors_derot_eigen(r, 1) *
                                                matched_catalog_vectors_derot_eigen(r, 1) +
                                                matched_catalog_vectors_derot_eigen(r, 2) *
                                                matched_catalog_vectors_derot_eigen(r, 2));
                                            if (matched_catalog_vectors_derot_eigen(r, 0) > 1e-9) {
                                                tangent_matched_catalog_vectors(r) = norm_yz /
                                                    matched_catalog_vectors_derot_eigen(r, 0);
                                            } else {
                                                tangent_matched_catalog_vectors(r) = 0.0f;
                                            }
                                        }

                                        Eigen::VectorXf radius_matched_image_centroids(
                                            matched_image_centroids.size());
                                        float img_center_y_f = static_cast<float>(height) / 2.0f;
                                        float img_center_x_f = static_cast<float>(width) / 2.0f;
                                        float half_width_f = static_cast<float>(width) / 2.0f;
                                        for (size_t r = 0; r < matched_image_centroids.size(); ++r) {
                                            float dx_pix = matched_image_centroids[r].x - img_center_x_f;
                                            float dy_pix = matched_image_centroids[r].y - img_center_y_f;
                                            radius_matched_image_centroids(r) = std::sqrt(
                                                dx_pix * dx_pix + dy_pix * dy_pix) / half_width_f;
                                        }

                                        Eigen::MatrixXf A_ls(num_star_matches, 2);
                                        A_ls.col(0) = tangent_matched_catalog_vectors;
                                        A_ls.col(1) = radius_matched_image_centroids.array().pow(3).matrix();

                                        Eigen::VectorXf b_ls = radius_matched_image_centroids;

                                        Eigen::VectorXf solution = A_ls.colPivHouseholderQr().solve(b_ls);
                                        float f_focal = solution(0);
                                        k_distortion = solution(1);

                                        if (std::abs(1.0f - k_distortion) > 1e-9) {
                                            f_focal = f_focal / (1.0f - k_distortion);
                                        }

                                        fov = 2.0 * std::atan(1.0 / f_focal);

                                        image_centroids_undist = _undistort_centroids(
                                            matched_image_centroids, height, width, k_distortion);
                                    }

                                    std::vector<std::array<double, 3>> final_match_vectors = compute_vectors(
                                        image_centroids_undist, height, width, fov);

                                    Eigen::MatrixXf final_match_vectors_eigen(final_match_vectors.size(), 3);
                                    for (size_t r = 0; r < final_match_vectors.size(); ++r) {
                                        for (int c = 0; c < 3; ++c) {
                                            final_match_vectors_eigen(r, c) = final_match_vectors[r][c];
                                        }
                                    }
                                    Eigen::MatrixXf rotated_final_match_vectors_eigen = final_match_vectors_eigen *
                                        rotation_matrix_eigen;

                                    std::vector<std::array<double, 3>> rotated_final_match_vectors_std(
                                        rotated_final_match_vectors_eigen.rows());
                                    for (int r = 0; r < rotated_final_match_vectors_eigen.rows(); ++r) {
                                        for (int c = 0; c < 3; ++c) {
                                            rotated_final_match_vectors_std[r][c] =
                                                    rotated_final_match_vectors_eigen(r, c);
                                        }
                                    }

                                    std::vector<float> angles_residual_rad;
                                    float sum_angle_sq = 0.0f;
                                    for (size_t res_idx = 0; res_idx < rotated_final_match_vectors_std.size(); ++res_idx) {
                                        float dist = vector_distance(
                                            rotated_final_match_vectors_std[res_idx],
                                            matched_catalog_vectors[res_idx]);
                                        float angle_rad = 2.0f * std::asin(std::min(1.0f, dist / 2.0f));
                                        angles_residual_rad.push_back(angle_rad);
                                        sum_angle_sq += angle_rad * angle_rad;
                                    }

                                    float rmse_rad = 0.0f;
                                    if (!angles_residual_rad.empty()) {
                                        rmse_rad = std::sqrt(
                                            sum_angle_sq / static_cast<float>(angles_residual_rad.size()));
                                    }
                                    float residual_arcsec = rmse_rad * 180.0f / M_PI * 3600.0f;

                                    auto end_time = std::chrono::high_resolution_clock::now();
                                    std::chrono::duration<double, std::milli> solve_time = end_time - start_time;

                                    std::cout << "Solution found!" << std::endl;
                                    std::cout << "RA: " << ra * 180.0f / M_PI << " degrees" << std::endl;
                                    std::cout << "Dec: " << dec * 180.0f / M_PI << " degrees" << std::endl;
                                    std::cout << "Roll: " << roll * 180.0f / M_PI << " degrees" << std::endl;
                                    std::cout << "FOV: " << fov * 180.0 / M_PI << " degrees" << std::endl;
                                    std::cout << "Matches: " << num_star_matches << std::endl;
                                    std::cout << "RMSE: " << std::fixed << std::setprecision(3) <<
                                            residual_arcsec << " arcsec" << std::endl;
                                    std::cout << "Solve time: " << std::fixed << std::setprecision(3) <<
                                            solve_time.count() << " ms" << std::endl;

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

size_t SimpleStarSolver::get_memory_usage() const {
    size_t star_memory = star_catalog.size() * sizeof(StarEntry);
    size_t pattern_memory = pattern_catalog.size() * sizeof(PatternEntry);
    return star_memory + pattern_memory;
}

void SimpleStarSolver::print_memory_usage() const {
    size_t total = get_memory_usage();
    std::cout << "Memory usage:\n";
    std::cout << "  Stars: " << star_catalog.size() << " entries, "
            << (star_catalog.size() * sizeof(StarEntry)) << " bytes\n";
    std::cout << "  Patterns: " << pattern_catalog.size() << " entries, "
            << (pattern_catalog.size() * sizeof(PatternEntry)) << " bytes\n";
    std::cout << "  Total: " << total << " bytes (" << total / 1024.0f << " KB)\n";
}
