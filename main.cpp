#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <string>
#include <zlib.h>
#include <iostream>
#include <set>
#include <Eigen/Dense>

// #define PLOT

#ifdef PLOT
#include <matplot/matplot.h>
#endif

struct Centroid
{
    double y, x; // pixel coordinates (y=row, x=col)
};

struct StarEntry
{
    float ra, dec; // radians
    float x, y, z; // unit vector components
    float magnitude;
};

struct PatternEntry
{
    uint16_t star_indices[4]; // indices into star catalog
};

struct SolveResult
{
    bool solved;
    float ra, dec, roll; // radians
    float fov; // horizontal FOV in radians
    float rmse; // RMS error in arcseconds
    int num_matches;
    float solve_time_ms;
    float distortion_k;
};

long double combinations(int n, int k)
{
    if (k < 0 || k > n)
    {
        return 0;
    }
    if (k == 0 || k == n)
    {
        return 1;
    }
    if (k > n / 2)
    {
        k = n - k; // C(n, k) = C(n, n-k)
    }
    long double res = 1.0L;
    for (int i = 1; i <= k; ++i)
    {
        // This order helps avoid large intermediate numbers (e.g., n!/k!)
        res = res * (n - i + 1) / i;
    }
    return res;
}

#ifdef PLOT
void plot_2d(std::vector<Centroid> data, std::string marker_color)
{
    using namespace matplot;
    std::vector<double> x_coords;
    std::vector<double> y_coords;

    for (const auto& point : data)
    {
        x_coords.push_back(point.x);
        y_coords.push_back(point.y);
    }

    auto s = scatter(x_coords, y_coords);
    s->marker_color(marker_color);
}
#endif

// Binomial CDF: P(X <= k) = sum_{i=0 to k} C(n, i) * p^i * (1-p)^(n-i)
long double calculate_binomial_cdf(int k, int n, long double p)
{
    if (k < 0)
        return 0.0L;
    if (k >= n)
        return 1.0L;
    if (p < 0.0L || p > 1.0L)
        return 0.0L; // Invalid probability

    long double cdf = 0.0L;
    long double one_minus_p = 1.0L - p;

    for (int i = 0; i <= k; ++i)
    {
        long double term = combinations(n, i);
        term *= std::pow(p, i);
        term *= std::pow(one_minus_p, n - i);
        cdf += term;
    }
    return cdf;
}


class SimpleStarSolver
{
private:
    std::vector<StarEntry> star_catalog;
    std::vector<PatternEntry> pattern_catalog;

    // Database properties
    int pattern_size = 4;
    int pattern_bins = 50; // ~1/4/0.005 from tetra3 defaults
    float pattern_max_error = 0.005f;
    float max_fov_deg = 30.0f;
    size_t num_patterns_in_catalog = 0; // New: To store total num patterns for mismatch prob

    static constexpr uint64_t MAGIC_RAND = 2654435761ULL;


    // Helper function to calculate powers, similar to numpy's bin_factor**np.arange(...)
    // This avoids repeated multiplication inside the loop and can be more numerically stable for large powers.
    std::vector<uint64_t> calculate_powers(int base, unsigned long count)
    {
        std::vector<uint64_t> powers(count);
        if (count > 0)
        {
            powers[0] = 1;
            for (int i = 1; i < count; ++i)
            {
                powers[i] = powers[i - 1] * base;
                // Optional: Add overflow check here if bin_factor can be very large
                // if (powers[i] / base != powers[i-1]) { /* handle overflow */ }
            }
        }
        return powers;
    }

    std::vector<uint64_t> key_to_index(const std::vector<std::vector<int>>& keys, int bin_factor,
                                       uint64_t max_index)
    {
        std::vector<uint64_t> hash_indices;
        if (keys.empty())
        {
            return hash_indices;
        }

        // Determine the length of the keys (assuming all inner vectors have the same size)
        size_t key_length = keys[0].size();
        if (key_length == 0)
        {
            return hash_indices;
        }

        std::vector<uint64_t> powers = calculate_powers(bin_factor, key_length);
        hash_indices.reserve(keys.size());

        for (const auto& key : keys)
        {
            uint64_t hash_val = 0;
            for (size_t i = 0; i < key_length; ++i)
            {
                // Cast both operands to uint64_t before multiplication
                hash_val += static_cast<uint64_t>(key[i]) * powers[i];
            }

            // Apply hashing and mod
            hash_val = (hash_val * MAGIC_RAND) % max_index;
            hash_indices.push_back(hash_val);
        }

        return hash_indices;
    }

    std::vector<Centroid> _undistort_centroids(
        const std::vector<Centroid>& centroids,
        int height, int width, float k_distortion)
    {
        std::vector<Centroid> undistorted_centroids;
        undistorted_centroids.reserve(centroids.size());

        float img_center_y = static_cast<float>(height) / 2.0f;
        float img_center_x = static_cast<float>(width) / 2.0f;
        float half_width_f = static_cast<float>(width) / 2.0f;

        for (const auto& c_distorted : centroids)
        {
            // Calculate distorted pixel distance from image center
            float dx = c_distorted.x - img_center_x;
            float dy = c_distorted.y - img_center_y;
            float r_distorted_pixels = std::sqrt(dx * dx + dy * dy);

            // Convert to scaled radius 'r_d' (relative to width/2) as in Python
            float r_d = r_distorted_pixels / half_width_f;

            // Apply undistortion formula: r_u = r_d * (1 - k * r_d^2)
            float r_u_scaled = r_d * (1.0f - k_distortion * r_d * r_d);

            // Calculate scaling factor for the delta X, Y components
            // This factor transforms the distorted pixel displacement (dx, dy)
            // to the undistorted pixel displacement.
            float scale_factor = (r_distorted_pixels > 1e-9) ? (r_u_scaled * half_width_f / r_distorted_pixels) : 0.0f;

            Centroid c_undistorted;
            c_undistorted.x = img_center_x + dx * scale_factor;
            c_undistorted.y = img_center_y + dy * scale_factor;
            undistorted_centroids.push_back(c_undistorted);
        }
        return undistorted_centroids;
    }

    std::vector<std::array<double, 3>> sort_pattern_by_centroid(const std::vector<std::array<double, 3>>& vectors)
    {
        if (vectors.empty())
        {
            return {};
        }

        // 1. Calculate the centroid
        std::array<double, 3> centroid = {0.0, 0.0, 0.0};
        for (const auto& vec : vectors)
        {
            centroid[0] += vec[0];
            centroid[1] += vec[1];
            centroid[2] += vec[2];
        }
        centroid[0] /= vectors.size();
        centroid[1] /= vectors.size();
        centroid[2] /= vectors.size();

        // 2. Create a vector of pairs: {distance_from_centroid, original_vector}
        std::vector<std::pair<double, std::array<double, 3>>> tagged_vectors;
        tagged_vectors.reserve(vectors.size());
        for (const auto& vec : vectors)
        {
            double dist = vector_distance(vec, centroid);
            tagged_vectors.push_back({dist, vec});
        }

        // 3. Sort based on distances, with lexicographical tie-breaking
        std::sort(tagged_vectors.begin(), tagged_vectors.end(),
                  [](const std::pair<double, std::array<double, 3>>& a,
                     const std::pair<double, std::array<double, 3>>& b)
                  {
                      if (a.first != b.first)
                      {
                          return a.first < b.first; // Sort by distance
                      }
                      // If distances are equal, sort lexicographically
                      if (a.second[0] != b.second[0]) return a.second[0] < b.second[0];
                      if (a.second[1] != b.second[1]) return a.second[1] < b.second[1];
                      return a.second[2] < b.second[2];
                  });

        // Extract the sorted vectors
        std::vector<std::array<double, 3>> sorted_vectors;
        sorted_vectors.reserve(vectors.size());
        for (const auto& tagged_vec : tagged_vectors)
        {
            sorted_vectors.push_back(tagged_vec.second);
        }

        return sorted_vectors;
    }


    /**
     * @brief Get unit vectors from star centroids (pinhole camera model).
     * @param centroids A vector of Centroid structs, each containing (y, x) pixel coordinates.
     * @param height The height of the image in pixels.
     * @param width The width of the image in pixels.
     * @param fov The field-of-view in the x dimension (in radians).
     * @return A vector of std::array<float, 3>, where each array is a 3D unit vector [i, j, k].
     */
    std::vector<std::array<double, 3>> compute_vectors(
        const std::vector<Centroid>& centroids,
        int height, int width, double fov) // fov should also be double for consistency
    {
        // Ensure width is not zero to avoid division by zero
        if (width == 0)
        {
            return {}; // Return empty vector if width is invalid
        }

        // Calculate the scaling factor based on FOV and image width
        // Use double for all calculations to match Python's effective precision
        double scale_factor = std::tan(fov / 2.0) / static_cast<double>(width) * 2.0;

        // Initialize the vector to store the 3D star vectors
        std::vector<std::array<double, 3>> star_vectors(centroids.size());

        // Calculate the pixel center of the image
        double img_center_y = static_cast<double>(height) / 2.0;
        double img_center_x = static_cast<double>(width) / 2.0;

        // Iterate through each centroid to compute its corresponding 3D vector
        for (size_t i = 0; i < centroids.size(); ++i)
        {
            // The first component (index 0) is always 1.0 for the pinhole camera model
            star_vectors[i][0] = 1.0;

            // Calculate the difference from the image center and apply the scale factor
            // Cast centroid values to double for calculations
            star_vectors[i][2] = (img_center_y - static_cast<double>(centroids[i].y)) * scale_factor;
            star_vectors[i][1] = (img_center_x - static_cast<double>(centroids[i].x)) * scale_factor;

            // Calculate the Euclidean norm (magnitude) of the current vector
            double norm = std::sqrt(
                star_vectors[i][0] * star_vectors[i][0] +
                star_vectors[i][1] * star_vectors[i][1] +
                star_vectors[i][2] * star_vectors[i][2]
            );

            // Normalize the vector to make it a unit vector
            // Use a smaller epsilon for double precision comparisons
            if (norm > 1e-12) // A common small epsilon for double precision
            {
                star_vectors[i][0] /= norm;
                star_vectors[i][1] /= norm;
                star_vectors[i][2] /= norm;
            }
            else
            {
                // If norm is zero or very small, set vector to zero (or handle as an error)
                star_vectors[i][0] = 0.0;
                star_vectors[i][1] = 0.0;
                star_vectors[i][2] = 0.0;
            }
        }

        return star_vectors;
    }

    // Calculate distance between two unit vectors
    double vector_distance(const std::array<double, 3>& v1, const std::array<double, 3>& v2)
    {
        double diff_x = v1[0] - v2[0];
        double diff_y = v1[1] - v2[1];
        double diff_z = v1[2] - v2[2];
        return std::sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    }

    // Calculate angle between two unit vectors (more accurate than distance for small angles)
    double vector_angle(const std::array<double, 3>& v1, const std::array<double, 3>& v2)
    {
        double dist = vector_distance(v1, v2);
        return 2.0f * std::asin(dist / 2.0f);
    }

    // Calculate edge ratios for a pattern of 4 stars
    std::vector<double> calculate_edge_ratios(const std::vector<std::array<double, 3>>& pattern_vectors)
    {
        std::vector<double> edge_angles;

        // Calculate all pairwise distances (6 edges for 4 stars)
        for (int i = 0; i < 4; i++)
        {
            for (int j = i + 1; j < 4; j++)
            {
                float angle = vector_angle(pattern_vectors[i], pattern_vectors[j]);
                edge_angles.push_back(angle);
            }
        }

        // Sort and normalize by largest edge
        std::sort(edge_angles.begin(), edge_angles.end());
        float largest_edge = edge_angles.back();

        std::vector<double> edge_ratios;
        for (int i = 0; i < 5; i++)
        {
            // First 5 ratios (excluding largest/largest = 1)
            edge_ratios.push_back(edge_angles[i] / largest_edge);
        }

        return edge_ratios;
    }


    /**
     * @brief Generates a list of all possible hash codes (bin combinations)
     * within a given tolerance range for a set of edge ratios.
     * @param edge_ratios The 5 sorted edge ratios for the current image pattern.
     * @param pattern_bins The number of bins used for hashing (p_bins in Python).
     * @param tolerance A float representing the allowed deviation for each ratio (e.g., 0.01 for 1% tolerance).
     * @return A vector of unique, sorted arrays of hash bins, ready to be converted to hash indices.
     */
    std::vector<std::vector<int>> generate_hash_code_combinations(
        const std::vector<double>& edge_ratios,
        int pattern_bins,
        float tolerance)
    {
        if (edge_ratios.size() != 5)
        {
            // This function is specifically for 5 edge ratios, as implied by the Python.
            // Adjust if your pattern size changes.
            return {};
        }

        // Calculate min/max bounds for each edge ratio
        std::vector<int> hash_code_space_min(5);
        std::vector<int> hash_code_space_max(5);

        for (size_t i = 0; i < edge_ratios.size(); ++i)
        {
            // image_pattern_edge_ratio_min = np.maximum(0, image_pattern_edge_ratio_min*p_bins).astype(int)
            // image_pattern_edge_ratio_max = np.minimum(p_bins, image_pattern_edge_ratio_max*p_bins).astype(int)

            // Apply tolerance to get min/max ratio values
            double ratio_min_val = edge_ratios[i] - tolerance;
            double ratio_max_val = edge_ratios[i] + tolerance;

            hash_code_space_min[i] = std::max(0, static_cast<int>(ratio_min_val * pattern_bins));
            // Note: Python's `p_bins` is often the exclusive upper bound for bins (0 to p_bins-1)
            // np.minimum(p_bins, ...) implies max bin value is p_bins (which is out of range 0 to p_bins-1)
            // Let's assume p_bins corresponds to the count, so max bin index is p_bins - 1.
            // Adjust `pattern_bins - 1` if `pattern_bins` is intended as the exclusive upper limit.
            hash_code_space_max[i] = std::min(pattern_bins - 1, static_cast<int>(ratio_max_val * pattern_bins));
        }

        size_t num_dimensions = hash_code_space_min.size();

        // Create ranges for each dimension: range(low, high + 1)
        std::vector<std::vector<int>> hash_code_ranges(num_dimensions);
        for (size_t i = 0; i < num_dimensions; ++i)
        {
            for (int val = hash_code_space_min[i]; val <= hash_code_space_max[i]; ++val)
            {
                hash_code_ranges[i].push_back(val);
            }
        }

        // Generate Cartesian product (itertools.product equivalent)
        std::vector<std::vector<int>> all_combinations;

        // Recursive helper for Cartesian product
        std::function<void(size_t, std::vector<int>&)> generate_product =
            [&](size_t dimension, std::vector<int>& current_combination)
        {
            if (dimension == num_dimensions)
            {
                all_combinations.push_back(current_combination);
                return;
            }

            for (int val : hash_code_ranges[dimension])
            {
                current_combination.push_back(val);
                generate_product(dimension + 1, current_combination);
                current_combination.pop_back();
            }
        };

        std::vector<int> temp;
        generate_product(0, temp);

        // Sort each combination (axis=1 in numpy)
        for (auto& combo : all_combinations)
        {
            std::sort(combo.begin(), combo.end());
        }

        // Remove duplicates (np.unique equivalent)
        std::set<std::vector<int>> unique_combinations(all_combinations.begin(), all_combinations.end());

        return std::vector<std::vector<int>>(unique_combinations.begin(), unique_combinations.end());
    }

    std::array<std::array<double, 3>, 3> find_rotation_matrix(
        const std::vector<std::array<double, 3>>& image_vectors,
        const std::vector<std::array<double, 3>>& catalog_vectors)
    {
        // Convert to Eigen::MatrixXf
        // Each row of the Eigen matrix will be a vector.
        // image_vectors: Nx3 matrix, catalog_vectors: Nx3 matrix
        int num_vectors = image_vectors.size();

        if (num_vectors < 3)
        {
            // Need at least 3 points/vectors for a robust solution
            return {{{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}}};
        }

        Eigen::MatrixXf A(num_vectors, 3); // Image vectors (source)
        Eigen::MatrixXf B(num_vectors, 3); // Catalog vectors (target)

        for (int i = 0; i < num_vectors; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                A(i, j) = image_vectors[i][j];
                B(i, j) = catalog_vectors[i][j];
            }
        }

        // // Centroiding (optional for unit vectors from origin, but good practice for general case)
        // // For unit vectors from origin, centroids will be near zero, so centering won't change H.
        // Eigen::RowVector3f centroid_A = A.colwise().mean();
        // Eigen::RowVector3f centroid_B = B.colwise().mean();
        //
        // Eigen::MatrixXf A_centered = A.rowwise() - centroid_A;
        // Eigen::MatrixXf B_centered = B.rowwise() - centroid_B;

        // Calculate the covariance matrix H = A_centered^T * B_centered
        // A is Nx3, B is Nx3. A_centered.transpose() is 3xN. B_centered is Nx3.
        // H will be 3x3.
        Eigen::Matrix3f H = A.transpose() * B;

        // Perform Singular Value Decomposition (SVD) on H
        Eigen::BDCSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Matrix3f V = svd.matrixV();

        // Calculate the rotation matrix R = V * U^T
        Eigen::Matrix3f R = (V * U.transpose()).transpose();

        // Handle potential reflection (det(R) == -1)
        // This happens if one set of vectors is a reflection of the other.
        // If the determinant is -1, flip the sign of the last column of V.
        if (R.determinant() < 0)
        {
            V.col(2) *= -1; // Flip the sign of the last column (corresponds to smallest singular value)
            R = V * U.transpose(); // Recompute R
        }

        // Convert Eigen::Matrix3f back
        std::array<std::array<double, 3>, 3> rotation_matrix_std;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                rotation_matrix_std[i][j] = R(i, j);
            }
        }

        return rotation_matrix_std;
    }

    // Find stars near a given direction
    std::vector<int> get_nearby_stars(const std::array<float, 3>& center_vector, float radius)
    {
        std::vector<int> nearby_indices;
        float cos_radius = std::cos(radius);

        for (int i = 0; i < star_catalog.size(); i++)
        {
            const auto& star = star_catalog[i];
            float dot_product = center_vector[0] * star.x +
                center_vector[1] * star.y +
                center_vector[2] * star.z;

            if (dot_product > cos_radius)
            {
                nearby_indices.push_back(i);
            }
        }

        return nearby_indices;
    }

public:
    // Load star catalog (simplified - in practice you'd load from file)
    void load_star_catalog(const std::vector<StarEntry>& stars)
    {
        star_catalog = stars;
    }

    // Load pattern catalog (simplified - in practice you'd load from file)
    void load_pattern_catalog(const std::vector<PatternEntry>& patterns)
    {
        pattern_catalog = patterns;
        num_patterns_in_catalog = patterns.size(); // Set this here
    }

    // Updated probing method that matches tetra3's behavior exactly
    std::vector<uint64_t> get_table_indices_from_hash(uint64_t hash_index, uint64_t table_size)
    {
        std::vector<uint64_t> found_indices;

        for (uint64_t probe = 0; probe < table_size; probe++) // Probe until we find empty or searched all
        {
            uint64_t probed_index = (hash_index + probe * probe) % table_size;

            // Check if this slot is empty (all zeros)
            const auto& test_pattern = pattern_catalog[probed_index];
            bool is_empty = true;
            for (int j = 0; j < 4; j++)
            {
                if (test_pattern.star_indices[j] != 0)
                {
                    is_empty = false;
                    break;
                }
            }

            if (is_empty)
            {
                // Found empty slot - end of probe sequence for this hash
                break;
            }
            else
            {
                // Found a pattern - add to results
                found_indices.push_back(probed_index);
            }

            // Safety check to avoid infinite loops (though shouldn't be needed)
            if (probe > 1000)
            {
                std::cout << "Warning: Excessive probing, breaking at " << probe << std::endl;
                break;
            }
        }

        return found_indices;
    }


    std::vector<std::array<int, 2>> _find_centroid_matches(
        const std::vector<Centroid>& image_centroids,
        const std::vector<Centroid>& catalog_centroids,
        float match_radius_pixels)
    {
        std::vector<std::array<int, 2>> matched_stars;
        std::vector<bool> catalog_matched(catalog_centroids.size(), false); // To ensure 1:1 matching

        float match_radius_sq = match_radius_pixels * match_radius_pixels;

        for (size_t i = 0; i < image_centroids.size(); ++i)
        {
            float min_dist_sq = std::numeric_limits<float>::max(); // Initialize with a very large value
            int best_match_idx = -1;

            for (size_t j = 0; j < catalog_centroids.size(); ++j)
            {
                if (catalog_matched[j])
                    continue; // Skip if this catalog star is already matched

                float dx = image_centroids[i].x - catalog_centroids[j].x;
                float dy = image_centroids[i].y - catalog_centroids[j].y;
                float dist_sq = dx * dx + dy * dy;

                if (dist_sq < match_radius_sq && dist_sq < min_dist_sq)
                {
                    min_dist_sq = dist_sq;
                    best_match_idx = static_cast<int>(j);
                }
            }

            if (best_match_idx != -1)
            {
                matched_stars.push_back({static_cast<int>(i), best_match_idx});
                catalog_matched[best_match_idx] = true; // Mark as matched
            }
        }
        return matched_stars;
    }

    std::pair<std::vector<Centroid>, std::vector<bool>> _compute_centroids_from_vectors(
        const std::vector<std::array<float, 3>>& vectors_derot,
        int height, int width, float fov)
    {
        std::vector<Centroid> centroids;
        std::vector<bool> kept_mask(vectors_derot.size(), false);

        if (width == 0 || height == 0)
        {
            return {centroids, kept_mask};
        }

        float tan_half_fov = std::tan(fov / 2.0f);
        // Inverse of the scale factor used in compute_vectors
        float scale_factor_inv = static_cast<float>(width) / (2.0f * tan_half_fov);

        float img_center_y = static_cast<float>(height) / 2.0f;
        float img_center_x = static_cast<float>(width) / 2.0f;

        for (size_t i = 0; i < vectors_derot.size(); ++i)
        {
            const auto& vec = vectors_derot[i];
            // Ensure the x-component (aligned with boresight) is positive and not near zero.
            // This means the vector is pointing towards the camera.
            if (vec[0] <= 1e-6)
            {
                // Use a small epsilon to handle floating point comparisons
                continue; // Skip this star, it's outside valid camera view (e.g., pointing away or sideways)
            }

            // Calculate normalized projected coordinates (similar to vec[1]/vec[0] and vec[2]/vec[0])
            // These are effectively tangent(angle from boresight) in X and Y directions, scaled.
            float normalized_y = vec[2] / vec[0]; // Corresponds to (img_center_y - centroid.y) * scale_factor
            float normalized_x = vec[1] / vec[0]; // Corresponds to (img_center_x - centroid.x) * scale_factor

            Centroid c{};
            c.y = img_center_y - normalized_y * scale_factor_inv;
            c.x = img_center_x - normalized_x * scale_factor_inv;

            // Check if centroid is within image bounds
            if (c.x >= 0 && c.x < width && c.y >= 0 && c.y < height)
            {
                centroids.push_back(c);
                kept_mask[i] = true;
            }
        }
        return {centroids, kept_mask};
    }


    // Main solving function
    SolveResult solve_from_centroids(
        const std::vector<Centroid>& centroids,
        int height, int width,
        double fov_estimate_deg = 20.0f,
        int pattern_checking_stars = 8,
        float match_radius = 0.01, // Relative match radius (e.g., 0.01 for 1% of width)
        float match_threshold = 0.001, // Probability threshold for match acceptance
        std::optional<float> distortion_coeff_in = std::nullopt // Optional input distortion
    )
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        SolveResult result = {false, 0, 0, 0, 0, 0, 0, 0, 0.0f};

        if (centroids.size() < 4)
        {
            return result;
        }

        double fov = fov_estimate_deg * M_PI / 180.0f; // Current FOV estimate in radians
        double k_distortion = distortion_coeff_in.value_or(0.0f); // Use input or default to 0

        // Undistort image centroids using the initial/estimated distortion
        std::vector<Centroid> image_centroids_undist;
        if (distortion_coeff_in.has_value())
        {
            image_centroids_undist = _undistort_centroids(centroids, height, width, k_distortion);
        }
        else
        {
            image_centroids_undist = centroids; // If no distortion provided, use raw centroids
        }

        int num_extracted_stars_raw = centroids.size(); // Number of extracted stars
        int max_stars_to_check = std::min(num_extracted_stars_raw, pattern_checking_stars);

        for (int i = 0; i < max_stars_to_check - 3; i++)
        {
            for (int j = i + 1; j < max_stars_to_check - 2; j++)
            {
                for (int k = j + 1; k < max_stars_to_check - 1; k++)
                {
                    for (int l = k + 1; l < max_stars_to_check; l++)
                    {
                        std::vector<Centroid> pattern_centroids = {
                            centroids[i], centroids[j], centroids[k], centroids[l]
                        };
                        auto pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);
                        auto edge_ratios = calculate_edge_ratios(pattern_vectors);
                        auto hash_code_list = generate_hash_code_combinations(
                            edge_ratios, pattern_bins, pattern_max_error);

                        // std::cout << "Centroids are" << std::endl;
                        // for (auto& arr : pattern_centroids)
                        // {
                        //     std::cout << "[" << arr.y << ", " << arr.x << "]" << std::endl;
                        // }
                        //
                        // std::cout << "Fov is: " << fov << std::endl;
                        //
                        // std::cout << "hash code lists are [";
                        // for (auto& arr : hash_code_list)
                        // {
                        //     std::cout << "[";
                        //     for (auto& val : arr)
                        //     {
                        //         std::cout << val << ", ";
                        //     }
                        //     std::cout << "]" << std::endl;
                        // }
                        // std::cout << "]" << std::endl;


                        // Convert to hash indices
                        auto hash_indices = key_to_index(hash_code_list, pattern_bins, pattern_catalog.size());

                        for (auto hash_index : hash_indices)
                        {
                            std::vector<uint64_t> hash_match_indices = get_table_indices_from_hash(
                                hash_index, pattern_catalog.size());

                            if (hash_match_indices.empty())
                            {
                                continue;
                            }

                            // Check each pattern found at this hash location
                            for (uint64_t catalog_index : hash_match_indices)
                            {
                                // std::cout << "Checking hash index: " << catalog_index << std::endl;
                                if (catalog_index == 2891333)
                                {
                                    std::cout << "hit breakpoint" << std::endl;
                                }
                                const auto& catalog_pattern = pattern_catalog[catalog_index];

                                // Get catalog star vectors
                                std::vector<std::array<double, 3>> catalog_vectors;
                                bool valid_catalog_pattern = true;

                                for (int idx : catalog_pattern.star_indices)
                                {
                                    if (idx < star_catalog.size())
                                    {
                                        const auto& star = star_catalog[idx];
                                        catalog_vectors.push_back({star.x, star.y, star.z});
                                    }
                                    else
                                    {
                                        valid_catalog_pattern = false;
                                        break;
                                    }
                                }

                                if (!valid_catalog_pattern || catalog_vectors.size() != 4)
                                {
                                    continue;
                                }

                                // Calculate catalog edge ratios
                                auto catalog_edge_ratios = calculate_edge_ratios(catalog_vectors);

                                // Check if ratios match within tolerance
                                bool match = true;
                                float max_error = 0.0f;
                                for (size_t m = 0; m < edge_ratios.size() && m < catalog_edge_ratios.size(); m++)
                                {
                                    float error = std::abs(edge_ratios[m] - catalog_edge_ratios[m]);
                                    max_error = std::max(max_error, error);
                                    if (error > pattern_max_error)
                                    {
                                        match = false;
                                        break;
                                    }
                                }

                                if (match)
                                {
                                    std::cout << "Pattern match found! Max error: " << std::fixed <<
                                        std::setprecision(7) << max_error << std::endl;
                                    std::cout << "Match index is " << catalog_index << std::endl;
                                    std::cout << "Match star values are ";
                                    for (auto& a : catalog_pattern.star_indices)
                                    {
                                        std::cout << a << ", ";
                                    }
                                    std::cout << std::endl;

                                    // Refine FOV estimate using the matched pattern
                                    float catalog_largest_edge = 0;
                                    for (int p = 0; p < 4; p++)
                                    {
                                        for (int q = p + 1; q < 4; q++)
                                        {
                                            float angle = vector_angle(catalog_vectors[p], catalog_vectors[q]);
                                            if (angle > catalog_largest_edge)
                                            {
                                                catalog_largest_edge = angle;
                                            }
                                        }
                                    }

                                    float image_largest_edge = 0;
                                    for (int p = 0; p < 4; p++)
                                    {
                                        for (int q = p + 1; q < 4; q++)
                                        {
                                            float angle = vector_angle(pattern_vectors[p], pattern_vectors[q]);
                                            if (angle > image_largest_edge)
                                            {
                                                image_largest_edge = angle;
                                            }
                                        }
                                    }

                                    std::cout << "Largest edges - catalog: " << catalog_largest_edge << ", image: " <<
                                        image_largest_edge << std::endl;

                                    if (image_largest_edge > 0.001f)
                                    {
                                        // Avoid division by zero
                                        fov = catalog_largest_edge / image_largest_edge * fov;
                                    }

                                    std::cout << "Refined FOV: " << fov << std::endl;


                                    // recompute based on FOV estimate
                                    pattern_vectors = compute_vectors(pattern_centroids, height, width, fov);

                                    // Sort both patterns by consistent criteria
                                    auto sorted_image_vectors = sort_pattern_by_centroid(pattern_vectors);
                                    auto sorted_catalog_vectors = catalog_vectors;
                                    // vectors are pre-sorted sort_pattern_by_centroid(catalog_vectors);

                                    // Calculate initial rotation matrix
                                    auto rotation_matrix_std_array = find_rotation_matrix(
                                        sorted_image_vectors, sorted_catalog_vectors);

                                    // Convert to Eigen::Matrix3f for easier operations
                                    Eigen::Matrix3f rotation_matrix_eigen;
                                    for (int row = 0; row < 3; ++row)
                                    {
                                        for (int col = 0; col < 3; ++col)
                                        {
                                            rotation_matrix_eigen(row, col) = rotation_matrix_std_array[row][col];
                                        }
                                    }

                                    std::cout << "Rotation matrix:\n";
                                    std::cout << rotation_matrix_eigen << std::endl; // Use Eigen's stream operator

                                    // fov_diagonal_rad = fov * np.sqrt(width**2 + height**2) / width
                                    float fov_diagonal_rad = fov * std::sqrt(
                                            static_cast<float>(width) * width + static_cast<float>(height) * height) /
                                        static_cast<float>(width);
                                    float search_radius_rad = fov_diagonal_rad / 2.0f;

                                    // Find all star vectors inside the (diagonal) field of view for matching
                                    // image_center_vector = rotation_matrix[0, :]
                                    // Assuming rotation_matrix.row(0) is the camera boresight vector in catalog frame
                                    std::array<float, 3> image_center_vector_arr = {
                                        rotation_matrix_eigen(0, 0), rotation_matrix_eigen(0, 1),
                                        rotation_matrix_eigen(0, 2)
                                    };

                                    // nearby_star_inds = self._get_nearby_stars(image_center_vector, fov_diagonal_rad/2)
                                    std::vector<int> nearby_star_inds = get_nearby_stars(
                                        image_center_vector_arr, search_radius_rad);

                                    // nearby_star_vectors = self.star_table[nearby_star_inds, 2:5]
                                    std::vector<std::array<double, 3>> nearby_star_vectors;
                                    nearby_star_vectors.reserve(nearby_star_inds.size());
                                    for (int idx : nearby_star_inds)
                                    {
                                        if (idx < star_catalog.size())
                                        {
                                            const auto& star = star_catalog[idx];
                                            nearby_star_vectors.push_back({star.x, star.y, star.z});
                                        }
                                    }

                                    // Derotate nearby stars and get their (undistorted) centroids using coarse fov
                                    // nearby_star_vectors_derot = np.dot(rotation_matrix, nearby_star_vectors.T).T
                                    // In Eigen: nearby_star_vectors_derot = (R_eigen * nearby_star_vectors_eigen.transpose()).transpose()
                                    // Or more directly: nearby_star_vectors_eigen * R_eigen.transpose()
                                    Eigen::MatrixXf nearby_star_vectors_eigen(nearby_star_vectors.size(), 3);
                                    for (size_t r = 0; r < nearby_star_vectors.size(); ++r)
                                    {
                                        for (int c = 0; c < 3; ++c)
                                        {
                                            nearby_star_vectors_eigen(r, c) = nearby_star_vectors[r][c];
                                        }
                                    }
                                    Eigen::MatrixXf nearby_star_vectors_derot_eigen = (
                                        rotation_matrix_eigen * nearby_star_vectors_eigen.transpose()).transpose();

                                    std::vector<std::array<float, 3>> nearby_star_vectors_derot_std(
                                        nearby_star_vectors_derot_eigen.rows());
                                    for (int r = 0; r < nearby_star_vectors_derot_eigen.rows(); ++r)
                                    {
                                        for (int c = 0; c < 3; ++c)
                                        {
                                            nearby_star_vectors_derot_std[r][c] = nearby_star_vectors_derot_eigen(r, c);
                                        }
                                    }

                                    // (nearby_star_centroids, kept) = _compute_centroids(nearby_star_vectors_derot, (height, width), fov)
                                    auto centroids_kept_pair = _compute_centroids_from_vectors(
                                        nearby_star_vectors_derot_std, height, width, fov);
                                    std::vector<Centroid> nearby_star_centroids = centroids_kept_pair.first;
                                    const std::vector<bool>& kept_mask = centroids_kept_pair.second;

                                    // Filter nearby_star_vectors and nearby_star_inds using the kept mask
                                    std::vector<std::array<double, 3>> filtered_nearby_star_vectors;
                                    std::vector<int> filtered_nearby_star_inds;
                                    for (size_t idx = 0; idx < kept_mask.size(); ++idx)
                                    {
                                        if (kept_mask[idx])
                                        {
                                            filtered_nearby_star_vectors.push_back(nearby_star_vectors[idx]);
                                            filtered_nearby_star_inds.push_back(nearby_star_inds[idx]);
                                        }
                                    }
                                    nearby_star_vectors = filtered_nearby_star_vectors;
                                    nearby_star_inds = filtered_nearby_star_inds;

                                    // Only keep as many as the centroids, they should ideally both be the num_stars brightest
                                    // nearby_star_centroids = nearby_star_centroids[:len(image_centroids)]
                                    // nearby_star_vectors = nearby_star_vectors[:len(image_centroids)]
                                    // nearby_star_inds = nearby_star_inds[:len(image_centroids)]
                                    int num_image_centroids = image_centroids_undist.size();
                                    // Use undistorted for matching
                                    if (nearby_star_centroids.size() > num_image_centroids)
                                    {
                                        nearby_star_centroids.resize(num_image_centroids);
                                        nearby_star_vectors.resize(num_image_centroids);
                                        nearby_star_inds.resize(num_image_centroids);
                                    }


                                    // Match these centroids to the image
                                    // matched_stars = _find_centroid_matches(image_centroids_undist, nearby_star_centroids, width*match_radius)
                                    float match_radius_pixels = static_cast<float>(width) * match_radius;
                                    std::vector<std::array<int, 2>> matched_stars = _find_centroid_matches(
                                        image_centroids_undist, nearby_star_centroids, match_radius_pixels);

                                    float boresight_x = rotation_matrix_eigen(0, 0);
                                    float boresight_y = rotation_matrix_eigen(1, 0);
                                    float boresight_z = rotation_matrix_eigen(2, 0);

                                    float ra = std::atan2(boresight_y, boresight_x);
                                    if (ra < 0) ra += 2.0f * M_PI; // Ensure positive RA

                                    float dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));

                                    // roll = np.rad2deg(np.arctan2(rotation_matrix[1, 2], rotation_matrix[2, 2])) % 360
                                    // This is a specific roll extraction. It maps to R_eigen(1,2) and R_eigen(2,2).
                                    float roll = std::atan2(rotation_matrix_eigen(1, 2),
                                                            rotation_matrix_eigen(2, 2));
                                    if (roll < 0) roll += 2.0f * M_PI; // Ensure positive roll

                                    std::cout << "RA: " << ra * 180.0f / M_PI << " deg, Dec: " << dec * 180.0f /
                                        M_PI << " deg, Roll: " << roll * 180.0f / M_PI << " deg" << std::endl;

                                    #ifdef PLOT
                                    matplot::cla();
                                    matplot::hold(matplot::on);
                                    plot_2d(centroids, "g");
                                    plot_2d(image_centroids_undist, "b");
                                    plot_2d(nearby_star_centroids, "r");
                                    // plot_2d(pattern_centroids, "k");
                                    matplot::hold(matplot::off);
                                    matplot::show();
                                    #endif

                                    int num_extracted_stars = image_centroids_undist.size();
                                    // Use count of undistorted centroids
                                    int num_nearby_catalog_stars = nearby_star_centroids.size();
                                    int num_star_matches = matched_stars.size();

                                    std::cout << "Number of nearby stars: " << num_nearby_catalog_stars
                                        << ", total matched: " << num_star_matches << std::endl;

                                    // Probability that a single star is a mismatch (fraction of area that are stars)
                                    // prob_single_star_mismatch = num_nearby_catalog_stars * match_radius**2
                                    long double prob_single_star_mismatch = static_cast<long double>(
                                            num_nearby_catalog_stars) * static_cast<long double>(match_radius) *
                                        static_cast
                                        <long double>(match_radius);

                                    // Probability that this rotation matrix's set of matches happen randomly
                                    // we subtract two degrees of freedom
                                    // prob_mismatch = scipy.stats.binom.cdf(num_extracted_stars - (num_star_matches - 2),
                                    //                                        num_extracted_stars,
                                    //                                        1 - prob_single_star_mismatch)
                                    int k_binom = num_extracted_stars - (num_star_matches - 2);
                                    int n_binom = num_extracted_stars;
                                    long double p_binom = 1.0L - prob_single_star_mismatch;

                                    long double prob_mismatch = calculate_binomial_cdf(k_binom, n_binom, p_binom);

                                    std::cout << "Mismatch probability = " << std::scientific << std::setprecision(2) <<
                                        prob_mismatch
                                        << ", at FOV = " << std::fixed << std::setprecision(5) << fov * 180.0f / M_PI <<
                                        "deg" << std::endl;

                                    if (prob_mismatch < match_threshold)
                                    {
                                        std::cout << "MATCH ACCEPTED" << std::endl;
                                        std::cout << "Prob: " << std::scientific << std::setprecision(4) <<
                                            prob_mismatch
                                            << ", corr: " << prob_mismatch * static_cast<long double>(
                                                num_patterns_in_catalog) << std::endl;

                                        // Get the vectors for all matches in the image using coarse fov
                                        // matched_image_centroids = image_centroids[matched_stars[:, 0], :]
                                        std::vector<Centroid> matched_image_centroids;
                                        matched_image_centroids.reserve(num_star_matches);
                                        for (const auto& match_pair : matched_stars)
                                        {
                                            matched_image_centroids.push_back(image_centroids_undist[match_pair[0]]);
                                            // Use undistorted
                                        }

                                        // matched_image_vectors = _compute_vectors(matched_image_centroids, (height, width), fov)
                                        std::vector<std::array<double, 3>> matched_image_vectors = compute_vectors(
                                            matched_image_centroids, height, width, fov);

                                        // matched_catalog_vectors = nearby_star_vectors[matched_stars[:, 1], :]
                                        std::vector<std::array<double, 3>> matched_catalog_vectors;
                                        matched_catalog_vectors.reserve(num_star_matches);
                                        for (const auto& match_pair : matched_stars)
                                        {
                                            matched_catalog_vectors.push_back(nearby_star_vectors[match_pair[1]]);
                                        }

                                        // Recompute rotation matrix for more accuracy
                                        // rotation_matrix = _find_rotation_matrix(matched_image_vectors, matched_catalog_vectors)
                                        rotation_matrix_std_array = find_rotation_matrix(
                                            matched_image_vectors, matched_catalog_vectors);
                                        for (int row = 0; row < 3; ++row)
                                        {
                                            // Convert back to Eigen for further ops
                                            for (int col = 0; col < 3; ++col)
                                            {
                                                rotation_matrix_eigen(row, col) = rotation_matrix_std_array[row][col];
                                            }
                                        }

                                        // Extract right ascension, declination, and roll from rotation matrix
                                        // Boresight is assumed to be the first COLUMN of the rotation matrix in C++
                                        // to align with standard attitude representations.
                                        // If R is the rotation from image to catalog frame, and camera X-axis is boresight,
                                        // then the boresight in catalog frame is R * [1,0,0]^T = R.col(0)
                                        float boresight_x = rotation_matrix_eigen(0, 0);
                                        float boresight_y = rotation_matrix_eigen(0, 1);
                                        float boresight_z = rotation_matrix_eigen(0, 2);

                                        float ra = std::atan2(boresight_y, boresight_x);
                                        if (ra < 0) ra += 2.0f * M_PI; // Ensure positive RA

                                        float dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));

                                        // roll = np.rad2deg(np.arctan2(rotation_matrix[1, 2], rotation_matrix[2, 2])) % 360
                                        // This is a specific roll extraction. It maps to R_eigen(1,2) and R_eigen(2,2).
                                        float roll = std::atan2(rotation_matrix_eigen(1, 2),
                                                                rotation_matrix_eigen(2, 2));
                                        if (roll < 0) roll += 2.0f * M_PI; // Ensure positive roll

                                        std::cout << "Boresight vector (catalog frame): " << boresight_x << ", " <<
                                            boresight_y << ", " << boresight_z << std::endl;
                                        std::cout << "RA: " << ra * 180.0f / M_PI << " deg, Dec: " << dec * 180.0f /
                                            M_PI << " deg, Roll: " << roll * 180.0f / M_PI << " deg" << std::endl;

                                        if (!distortion_coeff_in.has_value())
                                        {
                                            // No distortion model given initially
                                            // Compare mutual angles in catalogue to those with current
                                            // FOV estimate in order to scale accurately for fine FOV
                                            // angles_camera = 2 * np.arcsin(0.5 * pdist(matched_image_vectors))
                                            // angles_catalogue = 2 * np.arcsin(0.5 * pdist(matched_catalog_vectors))
                                            std::vector<float> angles_camera;
                                            std::vector<float> angles_catalogue;

                                            // pdist and 2 * np.arcsin(0.5 * dist) means calculating angular distance
                                            // between all pairs of vectors.
                                            for (size_t m1 = 0; m1 < matched_image_vectors.size(); ++m1)
                                            {
                                                for (size_t m2 = m1 + 1; m2 < matched_image_vectors.size(); ++m2)
                                                {
                                                    angles_camera.push_back(vector_angle(
                                                        matched_image_vectors[m1], matched_image_vectors[m2]));
                                                    angles_catalogue.push_back(vector_angle(
                                                        matched_catalog_vectors[m1], matched_catalog_vectors[m2]));
                                                }
                                            }

                                            if (!angles_camera.empty())
                                            {
                                                float mean_ratio = 0.0f;
                                                int count_ratios = 0;
                                                for (size_t ratio_idx = 0; ratio_idx < angles_camera.size(); ++
                                                     ratio_idx)
                                                {
                                                    if (angles_camera[ratio_idx] > 1e-9)
                                                    {
                                                        // Avoid division by zero for very small angles
                                                        mean_ratio += angles_catalogue[ratio_idx] / angles_camera[
                                                            ratio_idx];
                                                        count_ratios++;
                                                    }
                                                }
                                                if (count_ratios > 0)
                                                {
                                                    fov *= (mean_ratio / count_ratios); // Refine FOV
                                                }
                                            }
                                            k_distortion = 0.0f; // No distortion applied
                                            image_centroids_undist = matched_image_centroids;
                                            // Still undistorted, no further change
                                        }
                                        else
                                        {
                                            // Accurately calculate the FOV and distortion by looking at the angle from boresight
                                            // on all matched catalogue vectors and all matched image centroids

                                            // matched_catalog_vectors_derot = np.dot(rotation_matrix, matched_catalog_vectors.T).T
                                            Eigen::MatrixXf matched_catalog_vectors_eigen(
                                                matched_catalog_vectors.size(), 3);
                                            for (size_t r = 0; r < matched_catalog_vectors.size(); ++r)
                                            {
                                                for (int c = 0; c < 3; ++c)
                                                {
                                                    matched_catalog_vectors_eigen(r, c) = matched_catalog_vectors[r][c];
                                                }
                                            }
                                            Eigen::MatrixXf matched_catalog_vectors_derot_eigen =
                                                matched_catalog_vectors_eigen * rotation_matrix_eigen.transpose();

                                            // tangent_matched_catalog_vectors = norm(matched_catalog_vectors_derot[:, 1:], axis=1) / matched_catalog_vectors_derot[:, 0]
                                            // This is arctan(angle from boresight X-axis) if X is boresight
                                            Eigen::VectorXf tangent_matched_catalog_vectors(
                                                matched_catalog_vectors_derot_eigen.rows());
                                            for (int r = 0; r < matched_catalog_vectors_derot_eigen.rows(); ++r)
                                            {
                                                // sqrt(y^2 + z^2) / x
                                                float norm_yz = std::sqrt(
                                                    matched_catalog_vectors_derot_eigen(r, 1) *
                                                    matched_catalog_vectors_derot_eigen(r, 1) +
                                                    matched_catalog_vectors_derot_eigen(r, 2) *
                                                    matched_catalog_vectors_derot_eigen(r, 2));
                                                if (matched_catalog_vectors_derot_eigen(r, 0) > 1e-9)
                                                {
                                                    // Avoid division by zero
                                                    tangent_matched_catalog_vectors(r) = norm_yz /
                                                        matched_catalog_vectors_derot_eigen(r, 0);
                                                }
                                                else
                                                {
                                                    tangent_matched_catalog_vectors(r) = 0.0f;
                                                    // Or some large value indicating far off-axis
                                                }
                                            }

                                            // Get the (distorted) pixel distance from image centre for all matches
                                            // (scaled relative to width/2)
                                            // radius_matched_image_centroids = norm(matched_image_centroids - [height/2, width/2], axis=1)/width*2
                                            Eigen::VectorXf radius_matched_image_centroids(
                                                matched_image_centroids.size());
                                            float img_center_y_f = static_cast<float>(height) / 2.0f;
                                            float img_center_x_f = static_cast<float>(width) / 2.0f;
                                            float half_width_f = static_cast<float>(width) / 2.0f;
                                            for (size_t r = 0; r < matched_image_centroids.size(); ++r)
                                            {
                                                float dx_pix = matched_image_centroids[r].x - img_center_x_f;
                                                float dy_pix = matched_image_centroids[r].y - img_center_y_f;
                                                radius_matched_image_centroids(r) = std::sqrt(
                                                    dx_pix * dx_pix + dy_pix * dy_pix) / half_width_f;
                                            }


                                            // Solve system of equations in RMS sense for focal length f and distortion k
                                            // where f is focal length in units of image width/2
                                            // and k is distortion at width/2 (negative is barrel)
                                            // undistorted = distorted*(1 - k*(distorted*2/width)^2)
                                            // In form A * [f; k] = b
                                            // b = radius_matched_image_centroids
                                            // A = [tangent_matched_catalog_vectors, radius_matched_image_centroids^3]
                                            Eigen::MatrixXf A_ls(num_star_matches, 2);
                                            A_ls.col(0) = tangent_matched_catalog_vectors;
                                            A_ls.col(1) = radius_matched_image_centroids.array().pow(3).matrix();
                                            // element-wise power

                                            Eigen::VectorXf b_ls = radius_matched_image_centroids;

                                            // (f, k) = lstsq(A, b, rcond=None)[0].flatten()
                                            Eigen::VectorXf solution = A_ls.colPivHouseholderQr().solve(b_ls);
                                            float f_focal = solution(0); // focal length in units of image width/2
                                            k_distortion = solution(1); // distortion at width/2

                                            // Correct focal length to be at horizontal FOV
                                            // f = f/(1 - k)
                                            if (std::abs(1.0f - k_distortion) > 1e-9)
                                            {
                                                // Avoid division by zero
                                                f_focal = f_focal / (1.0f - k_distortion);
                                            }

                                            std::cout << "Calculated focal length to " << std::fixed <<
                                                std::setprecision(2) << f_focal
                                                << " and distortion to " << std::fixed << std::setprecision(3) <<
                                                k_distortion << std::endl;

                                            // Calculate (horizontal) true field of view
                                            // fov = 2*np.arctan(1/f)
                                            fov = 2.0f * std::atan(1.0f / f_focal);

                                            // Undistort centroids for final calculations
                                            // matched_image_centroids_undist = _undistort_centroids(matched_image_centroids, (height, width), k)
                                            image_centroids_undist = _undistort_centroids(
                                                matched_image_centroids, height, width, k_distortion);
                                        }

                                        // Get vectors
                                        // final_match_vectors = _compute_vectors(matched_image_centroids_undist, (height, width), fov)
                                        std::vector<std::array<double, 3>> final_match_vectors = compute_vectors(
                                            image_centroids_undist, height, width,
                                            fov); // Use newly computed FOV and undistorted centroids

                                        // Rotate to the sky (from image frame to catalog frame)
                                        // final_match_vectors = np.dot(rotation_matrix.T, final_match_vectors.T).T
                                        // This means: new_vec = R_transpose * old_vec
                                        Eigen::MatrixXf final_match_vectors_eigen(final_match_vectors.size(), 3);
                                        for (size_t r = 0; r < final_match_vectors.size(); ++r)
                                        {
                                            for (int c = 0; c < 3; ++c)
                                            {
                                                final_match_vectors_eigen(r, c) = final_match_vectors[r][c];
                                            }
                                        }
                                        Eigen::MatrixXf rotated_final_match_vectors_eigen = final_match_vectors_eigen *
                                            rotation_matrix_eigen; // Eigen R is image_to_catalog, so original_vec * R

                                        std::vector<std::array<double, 3>> rotated_final_match_vectors_std(
                                            rotated_final_match_vectors_eigen.rows());
                                        for (int r = 0; r < rotated_final_match_vectors_eigen.rows(); ++r)
                                        {
                                            for (int c = 0; c < 3; ++c)
                                            {
                                                rotated_final_match_vectors_std[r][c] =
                                                    rotated_final_match_vectors_eigen(r, c);
                                            }
                                        }

                                        // Calculate residual angles with more accurate formula
                                        // distance = norm(final_match_vectors - matched_catalog_vectors, axis=1)
                                        // angle = 2 * np.arcsin(.5 * distance)
                                        std::vector<float> angles_residual_rad;
                                        float sum_angle_sq = 0.0f;
                                        for (size_t res_idx = 0; res_idx < rotated_final_match_vectors_std.size(); ++
                                             res_idx)
                                        {
                                            float dist = vector_distance(
                                                rotated_final_match_vectors_std[res_idx],
                                                matched_catalog_vectors[res_idx]);
                                            float angle_rad = 2.0f * std::asin(std::min(1.0f, dist / 2.0f));
                                            // Clamp argument to asin to avoid NaN
                                            angles_residual_rad.push_back(angle_rad);
                                            sum_angle_sq += angle_rad * angle_rad;
                                        }

                                        // residual = np.rad2deg(np.sqrt(np.mean(angle**2))) * 3600
                                        float rmse_rad = 0.0f;
                                        if (!angles_residual_rad.empty())
                                        {
                                            rmse_rad = std::sqrt(
                                                sum_angle_sq / static_cast<float>(angles_residual_rad.size()));
                                        }
                                        float residual_arcsec = rmse_rad * 180.0f / M_PI * 3600.0f;

                                        std::cout << "Residual: " << std::fixed << std::setprecision(3) <<
                                            residual_arcsec << " arcsec" << std::endl;
                                        std::cout << "FOV estimate: " << fov_estimate_deg << " deg, refined: " << fov *
                                            180.0f / M_PI << " deg" << std::endl;
                                        std::cout << "Solution found!" << std::endl;
                                        std::cout << "RA: " << ra * 180.0f / M_PI << " degrees" << std::endl;
                                        std::cout << "Dec: " << dec * 180.0f / M_PI << " degrees" << std::endl;
                                        std::cout << "Roll: " << roll * 180.0f / M_PI << " degrees" << std::endl;
                                        std::cout << "FOV: " << fov * 180.0f / M_PI << " degrees" << std::endl;
                                        std::cout << "Matches: " << num_star_matches << std::endl;

                                        auto end_time = std::chrono::high_resolution_clock::now();
                                        std::chrono::duration<double, std::milli> solve_time = end_time - start_time;
                                        std::cout << "Solve time: " << std::fixed << std::setprecision(3) << solve_time.
                                            count() << " ms" << std::endl;

                                        result.solved = true;
                                        result.ra = ra;
                                        result.dec = dec;
                                        result.roll = roll;
                                        result.fov = fov;
                                        result.rmse = residual_arcsec;
                                        result.num_matches = num_star_matches;
                                        result.solve_time_ms = solve_time.count();
                                        result.distortion_k = k_distortion;
                                        return result; // Solution found, return
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
        return result; // No solution found;
    }

    // Get memory usage estimate
    size_t get_memory_usage() const
    {
        size_t star_memory = star_catalog.size() * sizeof(StarEntry);
        size_t pattern_memory = pattern_catalog.size() * sizeof(PatternEntry);
        return star_memory + pattern_memory;
    }

    void print_memory_usage() const
    {
        size_t total = get_memory_usage();
        std::cout << "Memory usage:\n";
        std::cout << "  Stars: " << star_catalog.size() << " entries, "
            << (star_catalog.size() * sizeof(StarEntry)) << " bytes\n";
        std::cout << "  Patterns: " << pattern_catalog.size() << " entries, "
            << (pattern_catalog.size() * sizeof(PatternEntry)) << " bytes\n";
        std::cout << "  Total: " << total << " bytes (" << total / 1024.0f << " KB)\n";
    }
};

class BinaryDatabaseLoader
{
private:
    std::vector<StarEntry> star_catalog;
    std::vector<PatternEntry> pattern_catalog;

public:
    bool load_database(const std::string& stars_file, const std::string& patterns_file)
    {
        return load_stars(stars_file) && load_patterns(patterns_file);
    }

    const std::vector<StarEntry>& get_stars() const { return star_catalog; }
    const std::vector<PatternEntry>& get_patterns() const { return pattern_catalog; }

    void print_info() const
    {
        size_t star_memory = star_catalog.size() * sizeof(StarEntry);
        size_t pattern_memory = pattern_catalog.size() * sizeof(PatternEntry);
        size_t total_memory = star_memory + pattern_memory;

        std::cout << "Loaded " << star_catalog.size() << " stars, "
            << pattern_catalog.size() << " patterns" << std::endl;
        std::cout << "Memory usage: " << total_memory / 1024.0f / 1024.0f << " MB" << std::endl;
    }

private:
    bool load_stars(const std::string& filename)
    {
        std::ifstream file(filename, std::ios::binary);
        if (!file) return false;

        uint32_t num_stars;
        file.read(reinterpret_cast<char*>(&num_stars), sizeof(uint32_t));
        if (!file) return false;

        star_catalog.resize(num_stars);
        file.read(reinterpret_cast<char*>(star_catalog.data()),
                  num_stars * sizeof(StarEntry));

        // Validate loaded data
        int invalid_count = 0;
        for (size_t i = 0; i < std::min(star_catalog.size(), (size_t)10); i++)
        {
            const auto& star = star_catalog[i];
            if (std::isnan(star.ra) || std::isnan(star.dec) ||
                std::isnan(star.x) || std::isnan(star.y) || std::isnan(star.z))
            {
                std::cout << "Invalid star " << i << ": RA=" << star.ra
                    << " Dec=" << star.dec << " xyz=(" << star.x
                    << "," << star.y << "," << star.z << ")" << std::endl;
                invalid_count++;
            }
        }

        if (invalid_count > 0)
        {
            std::cout << "Warning: Found " << invalid_count << " invalid stars in first 10" << std::endl;
        }

        return file.good() || file.eof();
    }

    bool load_patterns(const std::string& filename)
    {
        std::ifstream file(filename, std::ios::binary);
        if (!file) return false;

        uint32_t num_patterns;
        file.read(reinterpret_cast<char*>(&num_patterns), sizeof(uint32_t));
        if (!file) return false;

        pattern_catalog.resize(num_patterns);
        file.read(reinterpret_cast<char*>(pattern_catalog.data()),
                  num_patterns * sizeof(PatternEntry));

        return file.good() || file.eof();
    }
};

// Example usage
int main()
{
    BinaryDatabaseLoader loader;
    std::vector<StarEntry> stars;
    std::vector<PatternEntry> patterns;

    // Load the default tetra3 database
    std::string database_path =
        "/home/linus/git/OpenAstroExplorer/tracker/.venv/lib/python3.12/site-packages/tetra3/data/default_database.npz";
    // Adjust path as needed
    std::string stars_file =
        "/home/linus/git/OpenAstroExplorer/tracker/tetra3_db_stars.bin";
    std::string patterns_file =
        "/home/linus/git/OpenAstroExplorer/tracker/tetra3_db_patterns.bin";

    std::cout << "Loading tetra3 database: " << database_path << std::endl;

    std::vector<Centroid> centroids = {
        {692.4872f, 491.22507f},
        {946.2311f, 332.36243f},
        {12.089273f, 442.55933f},
        {943.6792f, 298.08932f},
        {734.17804f, 297.76077f},
        {1335.9685f, 518.3289f},
        {1374.8608f, 493.86456f},
        {905.30676f, 585.91943f},
        {1098.4895f, 316.72595f},
        {419.94678f, 1107.7972f},
        {989.1208f, 496.38684f},
        {850.3943f, 782.54486f},
        {1540.9384f, 747.3476f},
        {881.71277f, 99.38685f},
        {535.89136f, 267.47266f},
        {1441.419f, 593.14777f},
        {836.1878f, 981.1759f},
        {1409.5724f, 437.48322f},
        {618.7561f, 682.776f},
        {297.0253f, 1100.6609f},
        {307.83322f, 1160.8948f},
        {1111.1812f, 687.0469f},
        {307.18762f, 1002.8619f},
        {946.8009f, 63.03274f},
        {1386.8977f, 909.16187f},
        {827.1117f, 47.189194f},
        {1352.8467f, 956.95776f},
        {356.83426f, 704.9581f},
        {1309.0583f, 994.93353f},
        {1120.8969f, 1123.066f},
        {82.84231f, 652.91394f},
        {559.1041f, 1153.2173f},
        {458.92627f, 800.9313f},
        {1111.0157f, 873.15857f},
        {639.1404f, 724.95685f},
        {309.0613f, 1096.7659f},
        {425.03677f, 953.0256f},
        {568.9768f, 808.9929f},
        {94.97851f, 745.0461f},
        {1171.0835f, 492.99918f}
    };

    if (loader.load_database(stars_file, patterns_file))
    {
        std::cout << "Database loaded successfully!" << std::endl;

        // Now you can use these with your star solver
        SimpleStarSolver solver;
        solver.load_star_catalog(loader.get_stars());
        solver.load_pattern_catalog(loader.get_patterns());
        solver.print_memory_usage();
        auto result = solver.solve_from_centroids(centroids, 1200, 1600, 7.0f, 8);

        if (result.solved)
        {
            std::cout << "Solution found!\n";
            std::cout << "RA: " << result.ra << " degrees\n";
            std::cout << "Dec: " << result.dec << " degrees\n";
            std::cout << "Roll: " << result.roll << " degrees\n";
            std::cout << "FOV: " << result.fov << " degrees\n";
            std::cout << "Matches: " << result.num_matches << "\n";
            std::cout << "Solve time: " << result.solve_time_ms << " ms\n";
        }
        else
        {
            std::cout << "No solution found.\n";
            std::cout << "Solve time: " << result.solve_time_ms << " ms\n";
        }
        return 0;
    }
    else
    {
        std::cout << "Failed to load database!" << std::endl;
        return 1;
    }

    return 0;
    /// TARGET IS {'RA': 230.66761517248517, 'Dec': 11.03514732166332, 'Roll': 332.2796184441714, 'FOV': 11.424314488951488, 'distortion': 0.00070556672360746, 'RMSE': 5.408618080062131, 'Matches': 11, 'Prob': 5.2711548772017206e-14, 'epoch_equinox': 2000, 'epoch_proper_motion': 2025.0, 'T_solve': 7.733516002190299, 'T_extract': 81.35361300082877}
}
