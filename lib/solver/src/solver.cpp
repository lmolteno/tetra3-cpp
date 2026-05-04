#include <tetra3/solver.h>
#include <cstdlib>

namespace {
// Set TETRA3_DEBUG=1 in the environment to print a per-call funnel summary
// (combos tried, hash hits, ratio matches, binomial passes).
inline bool debug_enabled() {
    static const bool on = []{
        const char *v = std::getenv("TETRA3_DEBUG");
        return v && v[0] && v[0] != '0';
    }();
    return on;
}
}

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

std::vector<Centroid> SimpleStarSolver::_undistort_centroids_lens(
    const std::vector<Centroid> &centroids,
    int height, int width) {
    // Polynomial radial undistortion centered at the calibrated optical
    // axis. After this, every centroid is shifted to the coordinate system
    // where the geometric image center IS the optical axis, so the rest of
    // the solver (which uses width/2, height/2 as the projection origin)
    // works correctly without further changes.
    std::vector<Centroid> out;
    out.reserve(centroids.size());

    float img_center_x = static_cast<float>(width) / 2.0f;
    float img_center_y = static_cast<float>(height) / 2.0f;
    float cx_t = img_center_x + lens_cx_offset;
    float cy_t = img_center_y + lens_cy_offset;
    float half_w = img_center_x;
    float half_w2 = half_w * half_w;

    for (const auto &c : centroids) {
        float rx = static_cast<float>(c.x) - cx_t;
        float ry = static_cast<float>(c.y) - cy_t;
        float r2 = (rx * rx + ry * ry) / half_w2;
        // Clamp r²>1 to its boundary value to keep the polynomial well-
        // behaved past the calibrated region (corners reach r≈1.15).
        float r2_eval = r2 > 1.0f ? 1.0f : r2;
        float F = 0.0f;
        float p = r2_eval;
        for (float k : lens_k) {
            F += k * p;
            p *= r2_eval;
        }
        Centroid u;
        // Undistort to the optical axis frame, then translate to
        // geometric center. Net: x_out = (W/2) + rx*(1-F).
        u.x = img_center_x + rx * (1.0f - F);
        u.y = img_center_y + ry * (1.0f - F);
        out.push_back(u);
    }
    return out;
}

void SimpleStarSolver::set_lens_calibration(
    float cx_offset, float cy_offset,
    const std::vector<float> &k_radial) {
    lens_cx_offset = cx_offset;
    lens_cy_offset = cy_offset;
    lens_k = k_radial;
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
    float fov_max_error_deg,
    std::vector<MatchInfo> *matches_out,
    bool freeze_distortion,
    float extended_match_radius) {
    auto start_time = std::chrono::high_resolution_clock::now();
    SolveResult result = {false, 0, 0, 0, 0, 0, 0, 0, 0.0f};

    if (centroids.size() < 4) {
        return result;
    }

    double fov = fov_estimate_deg * M_PI / 180.0;
    const double fov_initial = fov;
    double fov_max_error_rad = fov_max_error_deg > 0
        ? fov_max_error_deg * M_PI / 180.0
        : 0.0;
    double k_distortion = distortion_coeff_in.value_or(0.0f);

    // Lens calibration takes precedence: pre-correct all input centroids
    // using the polynomial + offset center, putting the optical axis at the
    // geometric center. The rest of the solver then sees an "ideal pinhole"
    // and does NOT need to fit distortion per frame.
    bool lens_cal = has_lens_calibration();
    std::vector<Centroid> image_centroids_undist;
    if (lens_cal) {
        image_centroids_undist = _undistort_centroids_lens(centroids, height, width);
    } else if (distortion_coeff_in.has_value()) {
        image_centroids_undist = _undistort_centroids(centroids, height, width, k_distortion);
    } else {
        image_centroids_undist = centroids;
    }
    // Force freeze when lens calibration is set: the polynomial already
    // captures the lens, no per-frame k re-fit is wanted.
    if (lens_cal) {
        freeze_distortion = true;
        k_distortion = 0.0;
    }

    int num_extracted_stars_raw = centroids.size();
    int max_stars_to_check = std::min(num_extracted_stars_raw, pattern_checking_stars);

    // Debug funnel counters.
    long long dbg_combos = 0, dbg_hash_probes = 0, dbg_pattern_hits = 0;
    long long dbg_ratio_match = 0, dbg_fov_reject = 0;
    long long dbg_binomial_eval = 0, dbg_binomial_pass = 0;
    long double dbg_best_prob = 1.0L;
    int dbg_best_matches = 0;
    if (debug_enabled()) {
        std::cerr << "[tetra3] solve start: centroids=" << num_extracted_stars_raw
                  << " checking=" << max_stars_to_check
                  << " fov_est=" << fov_estimate_deg << "deg"
                  << " fov_max_err=" << fov_max_error_deg << "deg"
                  << " match_radius=" << match_radius
                  << " match_threshold=" << match_threshold
                  << " image=" << width << "x" << height
                  << "\n";
    }

    // Pre-compute star vectors for pattern-checking stars (avoids redundant work per combo)
    // Use undistorted centroids so distortion correction actually affects pattern matching
    std::vector<std::array<double, 3>> precomputed_vectors;
    {
        std::vector<Centroid> check_centroids(image_centroids_undist.begin(),
                                               image_centroids_undist.begin() + max_stars_to_check);
        precomputed_vectors = compute_vectors(check_centroids, height, width, fov);
    }

    // Pre-compute hash powers: pattern_bins^0 .. pattern_bins^4
    const uint64_t hp0 = 1;
    const uint64_t hp1 = pattern_bins;
    const uint64_t hp2 = hp1 * pattern_bins;
    const uint64_t hp3 = hp2 * pattern_bins;
    const uint64_t hp4 = hp3 * pattern_bins;
    const uint64_t table_size = pattern_catalog.size();

    for (int i = 0; i < max_stars_to_check - 3; i++) {
        for (int j = i + 1; j < max_stars_to_check - 2; j++) {
            for (int k = j + 1; k < max_stars_to_check - 1; k++) {
                for (int l = k + 1; l < max_stars_to_check; l++) {
                    dbg_combos++;
                    fov = fov_initial; // Reset for each pattern attempt
                    // Compute 6 edge angles from precomputed vectors
                    const auto &vi = precomputed_vectors[i];
                    const auto &vj = precomputed_vectors[j];
                    const auto &vk = precomputed_vectors[k];
                    const auto &vl = precomputed_vectors[l];

                    std::array<double, 6> edges = {
                        vector_angle(vi, vj), vector_angle(vi, vk), vector_angle(vi, vl),
                        vector_angle(vj, vk), vector_angle(vj, vl), vector_angle(vk, vl)
                    };
                    std::sort(edges.begin(), edges.end());
                    double largest_edge = edges[5];
                    if (largest_edge < 1e-9) continue;

                    std::array<double, 5> edge_ratios;
                    for (int d = 0; d < 5; d++) edge_ratios[d] = edges[d] / largest_edge;

                    // Compute hash bin ranges and generate indices inline (no heap alloc)
                    std::array<int, 5> bmin, bmax;
                    for (int d = 0; d < 5; d++) {
                        bmin[d] = std::max(0, static_cast<int>((edge_ratios[d] - pattern_max_error) * pattern_bins));
                        bmax[d] = std::min(pattern_bins - 1, static_cast<int>((edge_ratios[d] + pattern_max_error) * pattern_bins));
                    }

                    std::array<uint64_t, 64> hash_buf;
                    int n_hashes = 0;
                    for (int d0 = bmin[0]; d0 <= bmax[0]; d0++)
                    for (int d1 = bmin[1]; d1 <= bmax[1]; d1++)
                    for (int d2 = bmin[2]; d2 <= bmax[2]; d2++)
                    for (int d3 = bmin[3]; d3 <= bmax[3]; d3++)
                    for (int d4 = bmin[4]; d4 <= bmax[4]; d4++) {
                        std::array<int, 5> key = {d0, d1, d2, d3, d4};
                        std::sort(key.begin(), key.end());
                        uint64_t h = key[0]*hp0 + key[1]*hp1 + key[2]*hp2 + key[3]*hp3 + key[4]*hp4;
                        h = (h * MAGIC_RAND) % table_size;
                        if (n_hashes < 64) hash_buf[n_hashes++] = h;
                    }
                    std::sort(hash_buf.begin(), hash_buf.begin() + n_hashes);
                    n_hashes = static_cast<int>(std::unique(hash_buf.begin(), hash_buf.begin() + n_hashes) - hash_buf.begin());

                    for (int hi = 0; hi < n_hashes; hi++) {
                        // Inline hash probing (avoids vector allocation)
                        uint64_t hash_index = hash_buf[hi];
                        for (uint64_t probe = 0; probe < table_size && probe <= 1000; probe++) {
                            dbg_hash_probes++;
                            uint64_t probed_index = (hash_index + probe * probe) % table_size;
                            const auto &catalog_pattern = pattern_catalog[probed_index];

                            bool is_empty = true;
                            for (int x = 0; x < 4; x++) {
                                if (catalog_pattern.star_indices[x] != 0) { is_empty = false; break; }
                            }
                            if (is_empty) break;
                            dbg_pattern_hits++;

                            // Build catalog vectors
                            std::array<std::array<double, 3>, 4> cat_v;
                            bool valid = true;
                            for (int x = 0; x < 4; x++) {
                                auto idx = catalog_pattern.star_indices[x];
                                if (static_cast<size_t>(idx) < star_catalog.size()) {
                                    const auto &star = star_catalog[idx];
                                    cat_v[x] = {star.x, star.y, star.z};
                                } else { valid = false; break; }
                            }
                            if (!valid) continue;

                            // Verify edge ratios against catalog
                            std::vector<std::array<double, 3>> catalog_vectors(cat_v.begin(), cat_v.end());
                            auto catalog_edge_ratios = calculate_edge_ratios(catalog_vectors);

                            bool match = true;
                            for (size_t m = 0; m < 5 && m < catalog_edge_ratios.size(); m++) {
                                if (std::abs(edge_ratios[m] - catalog_edge_ratios[m]) > pattern_max_error) {
                                    match = false;
                                    break;
                                }
                            }

                            if (match) {
                                dbg_ratio_match++;
                                float catalog_largest_edge_val = 0;
                                for (int p = 0; p < 4; p++)
                                    for (int q = p + 1; q < 4; q++) {
                                        float angle = vector_angle(catalog_vectors[p], catalog_vectors[q]);
                                        if (angle > catalog_largest_edge_val) catalog_largest_edge_val = angle;
                                    }

                                if (largest_edge > 0.001) {
                                    double new_fov = catalog_largest_edge_val / largest_edge * fov_initial;
                                    if (fov_max_error_rad > 0 &&
                                        std::abs(new_fov - fov_estimate_deg * M_PI / 180.0) > fov_max_error_rad) {
                                        dbg_fov_reject++;
                                        continue;
                                    }
                                    fov = new_fov;
                                }

                                std::vector<Centroid> pattern_centroids = {
                                    image_centroids_undist[i], image_centroids_undist[j],
                                    image_centroids_undist[k], image_centroids_undist[l]
                                };
                                auto pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);

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

                                // === Iterative re-projection refinement ===
                                // Refit rotation from current matches, re-project all nearby
                                // catalog stars, re-match. Repeat until match count stabilizes.
                                // A real solve grows monotonically here (4 -> 7 -> 10); a
                                // spurious one stays flat or shrinks.
                                int prev_matches = -1;
                                const int max_refine_iter = 5;
                                for (int rit = 0; rit < max_refine_iter; rit++) {
                                    int curr = static_cast<int>(matched_stars.size());
                                    if (curr == prev_matches) break;
                                    if (curr < 4) break;
                                    prev_matches = curr;

                                    // Refit R from latest matches
                                    std::vector<Centroid> img_c;
                                    std::vector<std::array<double, 3>> cat_v;
                                    img_c.reserve(curr); cat_v.reserve(curr);
                                    for (const auto &mp : matched_stars) {
                                        img_c.push_back(image_centroids_undist[mp[0]]);
                                        cat_v.push_back(nearby_star_vectors[mp[1]]);
                                    }
                                    auto img_v = compute_vectors(img_c, height, width, fov);

                                    // Refit FOV from catalog/image angle ratios across
                                    // all matched pairs. Averages out the noise the
                                    // initial single-pattern FOV estimate carried.
                                    if (img_v.size() >= 2) {
                                        double ratio_sum = 0.0;
                                        int ratio_n = 0;
                                        for (size_t a = 0; a < img_v.size(); a++) {
                                            for (size_t b = a + 1; b < img_v.size(); b++) {
                                                double img_ang = vector_angle(img_v[a], img_v[b]);
                                                double cat_ang = vector_angle(cat_v[a], cat_v[b]);
                                                if (img_ang > 1e-9) {
                                                    ratio_sum += cat_ang / img_ang;
                                                    ratio_n++;
                                                }
                                            }
                                        }
                                        if (ratio_n > 0) {
                                            double new_fov = fov * (ratio_sum / ratio_n);
                                            // Recompute image vectors at refined FOV so R
                                            // is fit at consistent scale.
                                            if (std::abs(new_fov - fov) > 1e-6) {
                                                fov = new_fov;
                                                img_v = compute_vectors(img_c, height, width, fov);
                                            }
                                        }
                                    }

                                    auto rot = find_rotation_matrix(img_v, cat_v);
                                    for (int rr = 0; rr < 3; rr++)
                                        for (int cc = 0; cc < 3; cc++)
                                            rotation_matrix_eigen(rr, cc) = rot[rr][cc];

                                    // Search radius depends on FOV — refresh it
                                    float fov_diag_rad = fov * std::sqrt(
                                        static_cast<float>(width) * width +
                                        static_cast<float>(height) * height) /
                                        static_cast<float>(width);
                                    search_radius_rad = fov_diag_rad / 2.0f;

                                    // Re-find nearby catalog stars with the new boresight
                                    std::array<float, 3> bs = {
                                        rotation_matrix_eigen(0, 0),
                                        rotation_matrix_eigen(0, 1),
                                        rotation_matrix_eigen(0, 2)};
                                    auto new_inds = get_nearby_stars(bs, search_radius_rad);
                                    std::vector<std::array<double, 3>> new_vecs;
                                    new_vecs.reserve(new_inds.size());
                                    for (int idx : new_inds) {
                                        if (static_cast<size_t>(idx) < star_catalog.size()) {
                                            const auto &star = star_catalog[idx];
                                            new_vecs.push_back({star.x, star.y, star.z});
                                        }
                                    }
                                    if (new_vecs.empty()) break;

                                    Eigen::MatrixXf nve(new_vecs.size(), 3);
                                    for (size_t r = 0; r < new_vecs.size(); r++)
                                        for (int c = 0; c < 3; c++) nve(r, c) = new_vecs[r][c];
                                    Eigen::MatrixXf nde = (rotation_matrix_eigen * nve.transpose()).transpose();
                                    std::vector<std::array<float, 3>> nde_std(nde.rows());
                                    for (int r = 0; r < nde.rows(); r++)
                                        for (int c = 0; c < 3; c++) nde_std[r][c] = nde(r, c);

                                    auto pair2 = _compute_centroids_from_vectors(nde_std, height, width, fov);
                                    auto &new_centroids = pair2.first;
                                    const auto &new_kept = pair2.second;

                                    std::vector<std::array<double, 3>> kept_vecs;
                                    std::vector<int> kept_inds;
                                    for (size_t idx = 0; idx < new_kept.size(); idx++) {
                                        if (new_kept[idx]) {
                                            kept_vecs.push_back(new_vecs[idx]);
                                            kept_inds.push_back(new_inds[idx]);
                                        }
                                    }

                                    int nimg = static_cast<int>(image_centroids_undist.size());
                                    if (static_cast<int>(new_centroids.size()) > nimg) {
                                        new_centroids.resize(nimg);
                                        kept_vecs.resize(nimg);
                                        kept_inds.resize(nimg);
                                    }

                                    nearby_star_vectors = std::move(kept_vecs);
                                    nearby_star_inds    = std::move(kept_inds);
                                    nearby_star_centroids = std::move(new_centroids);
                                    matched_stars = _find_centroid_matches(
                                        image_centroids_undist, nearby_star_centroids, match_radius_pixels);
                                }

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

                                long double prob_single_star_mismatch = static_cast<long double>(
                                                                            num_nearby_catalog_stars) * static_cast<
                                                                            long double>(match_radius) *
                                                                        static_cast<long double>(match_radius);

                                int k_binom = num_extracted_stars - (num_star_matches - 2);
                                int n_binom = num_extracted_stars;
                                long double p_binom = 1.0L - prob_single_star_mismatch;

                                long double prob_mismatch = calculate_binomial_cdf(k_binom, n_binom, p_binom);
                                dbg_binomial_eval++;

                                // Require enough matches to trust the solve. Spurious solves
                                // typically plateau at ~4 matches even after refinement;
                                // real solves grow well past this floor.
                                int min_matches = std::max(
                                    6, static_cast<int>(0.3 * num_extracted_stars));

                                if (debug_enabled() &&
                                    (num_star_matches > dbg_best_matches ||
                                     (num_star_matches == dbg_best_matches && prob_mismatch < dbg_best_prob))) {
                                    dbg_best_matches = num_star_matches;
                                    dbg_best_prob = prob_mismatch;
                                    std::cerr << "[tetra3]   combo (" << i << "," << j << "," << k << "," << l
                                              << ") matches=" << num_star_matches
                                              << "/" << num_extracted_stars
                                              << " (min=" << min_matches << ")"
                                              << " nearby=" << num_nearby_catalog_stars
                                              << " fov=" << fov * 180.0 / M_PI << "deg"
                                              << " prob_mismatch=" << static_cast<double>(prob_mismatch)
                                              << "\n";
                                }

                                if (prob_mismatch < match_threshold &&
                                    num_star_matches >= min_matches) {
                                    dbg_binomial_pass++;

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

                                    if (!distortion_coeff_in.has_value() || freeze_distortion) {
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
                                        if (!freeze_distortion) {
                                            k_distortion = 0.0f;
                                        }
                                        // matched_image_centroids is already
                                        // undistorted by k_in if freeze is on
                                        // (centroids were undistorted at the
                                        // top of solve_from_centroids).
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

                                    // Outlier rejection: remove matches with residual > 3x median
                                    if (angles_residual_rad.size() > 6) {
                                        std::vector<float> sorted_res = angles_residual_rad;
                                        std::sort(sorted_res.begin(), sorted_res.end());
                                        float median_res = sorted_res[sorted_res.size() / 2];
                                        float outlier_thresh = 3.0f * std::max(median_res, 1e-6f);

                                        std::vector<Centroid> inlier_img_c;
                                        std::vector<std::array<double, 3>> inlier_cat_v;
                                        for (size_t idx = 0; idx < angles_residual_rad.size(); idx++) {
                                            if (angles_residual_rad[idx] <= outlier_thresh) {
                                                inlier_img_c.push_back(image_centroids_undist[idx]);
                                                inlier_cat_v.push_back(matched_catalog_vectors[idx]);
                                            }
                                        }

                                        if (inlier_img_c.size() >= 4 && inlier_img_c.size() < angles_residual_rad.size()) {
                                            int n_rejected = angles_residual_rad.size() - inlier_img_c.size();
                                            auto inlier_img_v = compute_vectors(inlier_img_c, height, width, fov);
                                            auto rot = find_rotation_matrix(inlier_img_v, inlier_cat_v);
                                            for (int rr = 0; rr < 3; rr++)
                                                for (int cc = 0; cc < 3; cc++)
                                                    rotation_matrix_eigen(rr, cc) = rot[rr][cc];

                                            boresight_x = rotation_matrix_eigen(0, 0);
                                            boresight_y = rotation_matrix_eigen(0, 1);
                                            boresight_z = rotation_matrix_eigen(0, 2);
                                            ra = std::atan2(boresight_y, boresight_x);
                                            if (ra < 0) ra += 2.0f * M_PI;
                                            dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));
                                            roll = std::atan2(rotation_matrix_eigen(1, 2), rotation_matrix_eigen(2, 2));
                                            if (roll < 0) roll += 2.0f * M_PI;

                                            // Recompute residuals from inliers
                                            Eigen::MatrixXf inlier_eigen(inlier_img_v.size(), 3);
                                            for (size_t r = 0; r < inlier_img_v.size(); r++)
                                                for (int c = 0; c < 3; c++)
                                                    inlier_eigen(r, c) = inlier_img_v[r][c];
                                            Eigen::MatrixXf rotated_inlier = inlier_eigen * rotation_matrix_eigen;

                                            sum_angle_sq = 0.0f;
                                            angles_residual_rad.clear();
                                            for (int r = 0; r < rotated_inlier.rows(); r++) {
                                                std::array<double, 3> rv = {rotated_inlier(r, 0), rotated_inlier(r, 1), rotated_inlier(r, 2)};
                                                float d = vector_distance(rv, inlier_cat_v[r]);
                                                float a = 2.0f * std::asin(std::min(1.0f, d / 2.0f));
                                                angles_residual_rad.push_back(a);
                                                sum_angle_sq += a * a;
                                            }

                                            num_star_matches = static_cast<int>(inlier_img_c.size());
                                            rmse_rad = std::sqrt(sum_angle_sq / angles_residual_rad.size());
                                            residual_arcsec = rmse_rad * 180.0f / M_PI * 3600.0f;
                                        }
                                    }

                                    if (matches_out && extended_match_radius > 0.0f) {
                                        // Rich match pass for distortion analysis.
                                        // Re-project ALL nearby catalog stars (no
                                        // num_image_centroids cap) using the final
                                        // R + fov + k_distortion, and match them
                                        // against ALL input centroids (undistorted)
                                        // with a generous radius.
                                        matches_out->clear();

                                        std::vector<Centroid> all_undist;
                                        if (lens_cal) {
                                            all_undist = _undistort_centroids_lens(
                                                centroids, height, width);
                                        } else if (k_distortion != 0.0f) {
                                            all_undist = _undistort_centroids(
                                                centroids, height, width, k_distortion);
                                        } else {
                                            all_undist = centroids;
                                        }

                                        std::array<float, 3> bs = {
                                            rotation_matrix_eigen(0, 0),
                                            rotation_matrix_eigen(0, 1),
                                            rotation_matrix_eigen(0, 2)};
                                        float fov_diag = static_cast<float>(fov) * std::sqrt(
                                            static_cast<float>(width) * width +
                                            static_cast<float>(height) * height) /
                                            static_cast<float>(width);
                                        auto rich_inds = get_nearby_stars(bs, fov_diag / 2.0f);

                                        std::vector<Centroid> rich_proj;
                                        std::vector<int> rich_proj_idx;
                                        rich_proj.reserve(rich_inds.size());
                                        rich_proj_idx.reserve(rich_inds.size());

                                        Eigen::MatrixXf rich_cat(rich_inds.size(), 3);
                                        for (size_t i = 0; i < rich_inds.size(); ++i) {
                                            const auto &s = star_catalog[rich_inds[i]];
                                            rich_cat(i, 0) = s.x;
                                            rich_cat(i, 1) = s.y;
                                            rich_cat(i, 2) = s.z;
                                        }
                                        Eigen::MatrixXf rich_cam =
                                            (rotation_matrix_eigen * rich_cat.transpose()).transpose();

                                        float scale_inv = static_cast<float>(width) /
                                            (2.0f * std::tan(static_cast<float>(fov) / 2.0f));
                                        float img_cy = static_cast<float>(height) / 2.0f;
                                        float img_cx = static_cast<float>(width) / 2.0f;

                                        for (size_t i = 0; i < rich_inds.size(); ++i) {
                                            float vx = rich_cam(i, 0);
                                            float vy = rich_cam(i, 1);
                                            float vz = rich_cam(i, 2);
                                            if (vx <= 1e-6f) continue;
                                            Centroid c{};
                                            c.x = img_cx - (vy / vx) * scale_inv;
                                            c.y = img_cy - (vz / vx) * scale_inv;
                                            if (c.x < 0 || c.x >= width || c.y < 0 || c.y >= height)
                                                continue;
                                            rich_proj.push_back(c);
                                            rich_proj_idx.push_back(rich_inds[i]);
                                        }

                                        float ext_pix = static_cast<float>(width) * extended_match_radius;
                                        auto ext_pairs = _find_centroid_matches(
                                            all_undist, rich_proj, ext_pix);

                                        matches_out->reserve(ext_pairs.size());
                                        for (const auto &mp : ext_pairs) {
                                            MatchInfo m{};
                                            int img_i = mp[0];
                                            int proj_i = mp[1];
                                            int cat_global = rich_proj_idx[proj_i];
                                            m.catalog_index = cat_global;
                                            const auto &s = star_catalog[cat_global];
                                            m.catalog_ra = s.ra;
                                            m.catalog_dec = s.dec;
                                            m.magnitude = s.magnitude;
                                            m.img_x = all_undist[img_i].x;
                                            m.img_y = all_undist[img_i].y;
                                            m.pred_x = rich_proj[proj_i].x;
                                            m.pred_y = rich_proj[proj_i].y;
                                            double dx = m.img_x - m.pred_x;
                                            double dy = m.img_y - m.pred_y;
                                            m.residual_pix = std::sqrt(dx * dx + dy * dy);
                                            matches_out->push_back(m);
                                        }
                                    } else if (matches_out) {
                                        // Standard: emit pairs from matched_stars only.
                                        matches_out->clear();
                                        matches_out->reserve(matched_stars.size());

                                        Eigen::MatrixXf cat_eigen(matched_stars.size(), 3);
                                        for (size_t mi = 0; mi < matched_stars.size(); ++mi) {
                                            const auto &nv = nearby_star_vectors[matched_stars[mi][1]];
                                            cat_eigen(mi, 0) = nv[0];
                                            cat_eigen(mi, 1) = nv[1];
                                            cat_eigen(mi, 2) = nv[2];
                                        }
                                        Eigen::MatrixXf cam_eigen =
                                            (rotation_matrix_eigen * cat_eigen.transpose()).transpose();

                                        float scale_inv = static_cast<float>(width) /
                                            (2.0f * std::tan(static_cast<float>(fov) / 2.0f));
                                        float img_cy = static_cast<float>(height) / 2.0f;
                                        float img_cx = static_cast<float>(width) / 2.0f;

                                        for (size_t mi = 0; mi < matched_stars.size(); ++mi) {
                                            MatchInfo m{};
                                            int nearby_idx = matched_stars[mi][1];
                                            int cat_idx = nearby_star_inds[nearby_idx];
                                            m.catalog_index = cat_idx;
                                            if (static_cast<size_t>(cat_idx) < star_catalog.size()) {
                                                const auto &s = star_catalog[cat_idx];
                                                m.catalog_ra = s.ra;
                                                m.catalog_dec = s.dec;
                                                m.magnitude = s.magnitude;
                                            }
                                            m.img_x = image_centroids_undist[mi].x;
                                            m.img_y = image_centroids_undist[mi].y;
                                            float vx = cam_eigen(mi, 0);
                                            float vy = cam_eigen(mi, 1);
                                            float vz = cam_eigen(mi, 2);
                                            if (vx > 1e-9f) {
                                                m.pred_x = img_cx - (vy / vx) * scale_inv;
                                                m.pred_y = img_cy - (vz / vx) * scale_inv;
                                            } else {
                                                m.pred_x = std::numeric_limits<double>::quiet_NaN();
                                                m.pred_y = std::numeric_limits<double>::quiet_NaN();
                                            }
                                            double dx = m.img_x - m.pred_x;
                                            double dy = m.img_y - m.pred_y;
                                            m.residual_pix = std::sqrt(dx * dx + dy * dy);
                                            matches_out->push_back(m);
                                        }
                                    }

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
    if (debug_enabled()) {
        std::cerr << "[tetra3] solve FAILED after " << solve_time.count() << "ms"
                  << "\n  combos_tried     = " << dbg_combos
                  << "\n  hash_probes      = " << dbg_hash_probes
                  << "\n  pattern_hits     = " << dbg_pattern_hits
                  << "\n  edge_ratio_match = " << dbg_ratio_match
                  << "\n  fov_rejected     = " << dbg_fov_reject
                  << "\n  binomial_eval    = " << dbg_binomial_eval
                  << "\n  binomial_pass    = " << dbg_binomial_pass
                  << "\n  best_match_count = " << dbg_best_matches
                  << "\n  best_prob        = " << static_cast<double>(dbg_best_prob)
                  << "\n";
    }
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
