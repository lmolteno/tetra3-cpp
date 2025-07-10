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

struct Centroid
{
    float y, x; // pixel coordinates (y=row, x=col)
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
};

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

    static constexpr uint64_t MAGIC_RAND = 2654435761ULL;


    // Helper function to calculate powers, similar to numpy's bin_factor**np.arange(...)
    // This avoids repeated multiplication inside the loop and can be more numerically stable for large powers.
    std::vector<uint64_t> calculate_powers(int base, int count)
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

        // Calculate powers once for all keys
        std::vector<uint64_t> powers = calculate_powers(bin_factor, key_length);

        hash_indices.reserve(keys.size()); // Reserve space for efficiency

        for (const auto& key : keys)
        {
            uint64_t hash_val = 0;
            for (size_t i = 0; i < key_length; ++i)
            {
                hash_val += static_cast<uint64_t>(key[i]) * powers[i];
            }
            hash_indices.push_back((hash_val * MAGIC_RAND) % max_index);
        }

        return hash_indices;
    }

    /**
     * @brief Get unit vectors from star centroids (pinhole camera model).
     * @param centroids A vector of Centroid structs, each containing (y, x) pixel coordinates.
     * @param height The height of the image in pixels.
     * @param width The width of the image in pixels.
     * @param fov The field-of-view in the x dimension (in radians).
     * @return A vector of std::array<float, 3>, where each array is a 3D unit vector [i, j, k].
     */
    std::vector<std::array<float, 3>> compute_vectors(
        const std::vector<Centroid>& centroids,
        int height, int width, float fov)
    {
        // Ensure width is not zero to avoid division by zero
        if (width == 0)
        {
            return {}; // Return empty vector if width is invalid
        }

        // Calculate the scaling factor based on FOV and image width
        // This is equivalent to (2 * tan(fov/2)) / width
        float scale_factor = std::tan(fov / 2.0f) / static_cast<float>(width) * 2.0f;

        // Initialize the vector to store the 3D star vectors
        std::vector<std::array<float, 3>> star_vectors(centroids.size());

        // Calculate the pixel center of the image
        float img_center_y = static_cast<float>(height) / 2.0f;
        float img_center_x = static_cast<float>(width) / 2.0f;

        // Iterate through each centroid to compute its corresponding 3D vector
        for (size_t i = 0; i < centroids.size(); ++i)
        {
            // The first component (index 0) is always 1.0 for the pinhole camera model
            star_vectors[i][0] = 1.0f;

            // Calculate the difference from the image center and apply the scale factor
            // The y-difference (row) maps to the 'k' component (index 2)
            // The x-difference (column) maps to the 'j' component (index 1)
            star_vectors[i][2] = (img_center_y - centroids[i].y) * scale_factor;
            star_vectors[i][1] = (img_center_x - centroids[i].x) * scale_factor;

            // Calculate the Euclidean norm (magnitude) of the current vector
            float norm = std::sqrt(
                star_vectors[i][0] * star_vectors[i][0] +
                star_vectors[i][1] * star_vectors[i][1] +
                star_vectors[i][2] * star_vectors[i][2]
            );

            // Normalize the vector to make it a unit vector
            // Avoid division by zero for very small norms
            if (norm > 1e-9)
            {
                // Using a small epsilon to check for near-zero norm
                star_vectors[i][0] /= norm;
                star_vectors[i][1] /= norm;
                star_vectors[i][2] /= norm;
            }
            else
            {
                // If norm is zero or very small, set vector to zero or handle as an error
                star_vectors[i][0] = 0.0f;
                star_vectors[i][1] = 0.0f;
                star_vectors[i][2] = 0.0f;
            }
        }

        return star_vectors;
    }

    // Calculate distance between two unit vectors
    float vector_distance(const std::array<float, 3>& v1, const std::array<float, 3>& v2)
    {
        float diff_x = v1[0] - v2[0];
        float diff_y = v1[1] - v2[1];
        float diff_z = v1[2] - v2[2];
        return std::sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
    }

    // Calculate angle between two unit vectors (more accurate than distance for small angles)
    float vector_angle(const std::array<float, 3>& v1, const std::array<float, 3>& v2)
    {
        float dist = vector_distance(v1, v2);
        return 2.0f * std::asin(dist / 2.0f);
    }

    // Calculate edge ratios for a pattern of 4 stars
    std::vector<float> calculate_edge_ratios(const std::vector<std::array<float, 3>>& pattern_vectors)
    {
        std::vector<float> edge_angles;

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

        std::vector<float> edge_ratios;
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
        const std::vector<float>& edge_ratios,
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
            float ratio_min_val = edge_ratios[i] * (1.0f - tolerance);
            float ratio_max_val = edge_ratios[i] * (1.0f + tolerance);

            hash_code_space_min[i] = std::max(0, static_cast<int>(ratio_min_val * pattern_bins));
            // Note: Python's `p_bins` is often the exclusive upper bound for bins (0 to p_bins-1)
            // np.minimum(p_bins, ...) implies max bin value is p_bins (which is out of range 0 to p_bins-1)
            // Let's assume p_bins corresponds to the count, so max bin index is p_bins - 1.
            // Adjust `pattern_bins - 1` if `pattern_bins` is intended as the exclusive upper limit.
            hash_code_space_max[i] = std::min(pattern_bins - 1, static_cast<int>(ratio_max_val * pattern_bins));
        }

        // Make an array of all combinations (Cartesian product)
        // This is the most complex part to replicate without a dedicated library like NumPy/itertools.
        // We'll use a recursive helper function or an iterative approach for Cartesian product.

        std::vector<std::vector<int>> hash_code_range_per_ratio(5);
        for (size_t i = 0; i < 5; ++i)
        {
            for (int bin = hash_code_space_min[i]; bin <= hash_code_space_max[i]; ++bin)
            {
                hash_code_range_per_ratio[i].push_back(bin);
            }
        }

        std::vector<std::vector<int>> all_combinations;
        if (hash_code_range_per_ratio.empty() || hash_code_range_per_ratio[0].empty())
        {
            return all_combinations; // No combinations if ranges are empty
        }

        // Recursive helper for Cartesian product
        std::function<void(int, std::vector<int>)> generate_product =
            [&](int k, std::vector<int> current_combination)
        {
            if (k == hash_code_range_per_ratio.size())
            {
                all_combinations.push_back(current_combination);
                return;
            }
            for (int val : hash_code_range_per_ratio[k])
            {
                current_combination.push_back(val);
                generate_product(k + 1, current_combination);
                current_combination.pop_back(); // Backtrack
            }
        };

        generate_product(0, {});

        // Make sure we have unique ascending codes
        // Python: hash_code_list = np.sort(hash_code_list, axis=1)
        // Python: hash_code_list = np.unique(hash_code_list, axis=0)

        // Sort each combination internally and then collect unique combinations
        std::set<std::vector<int>> unique_sorted_combinations;
        for (auto& combo : all_combinations)
        {
            std::sort(combo.begin(), combo.end()); // Sort each row
            unique_sorted_combinations.insert(combo); // Insert into set for uniqueness
        }

        return std::vector<std::vector<int>>(unique_sorted_combinations.begin(), unique_sorted_combinations.end());
    }

    // Find rotation matrix between two sets of vectors (simplified Kabsch algorithm)
    // You'll need to add a linear algebra library like Eigen for proper SVD
    // For now, here's a simplified but more correct approach:

    std::array<std::array<float, 3>, 3> find_rotation_matrix(
        const std::vector<std::array<float, 3>>& image_vectors,
        const std::vector<std::array<float, 3>>& catalog_vectors)
    {
        // Convert std::vector<std::array<float, 3>> to Eigen::MatrixXf
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

        // Centroiding (optional for unit vectors from origin, but good practice for general case)
        // For unit vectors from origin, centroids will be near zero, so centering won't change H.
        Eigen::RowVector3f centroid_A = A.colwise().mean();
        Eigen::RowVector3f centroid_B = B.colwise().mean();

        Eigen::MatrixXf A_centered = A.rowwise() - centroid_A;
        Eigen::MatrixXf B_centered = B.rowwise() - centroid_B;

        // Calculate the covariance matrix H = A_centered^T * B_centered
        // A is Nx3, B is Nx3. A_centered.transpose() is 3xN. B_centered is Nx3.
        // H will be 3x3.
        Eigen::Matrix3f H = A_centered.transpose() * B_centered;

        // Perform Singular Value Decomposition (SVD) on H
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Matrix3f V = svd.matrixV();

        // Calculate the rotation matrix R = V * U^T
        Eigen::Matrix3f R = V * U.transpose();

        // Handle potential reflection (det(R) == -1)
        // This happens if one set of vectors is a reflection of the other.
        // If the determinant is -1, flip the sign of the last column of V.
        if (R.determinant() < 0)
        {
            V.col(2) *= -1; // Flip the sign of the last column (corresponds to smallest singular value)
            R = V * U.transpose(); // Recompute R
        }

        // Convert Eigen::Matrix3f back to std::array<std::array<float, 3>, 3>
        std::array<std::array<float, 3>, 3> rotation_matrix_std;
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

    // Main solving function
    SolveResult solve_from_centroids(
        const std::vector<Centroid>& centroids,
        int height, int width,
        float fov_estimate_deg = 15.0f,
        int pattern_checking_stars = 8)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        SolveResult result = {false, 0, 0, 0, 0, 0, 0, 0};

        if (centroids.size() < 4)
        {
            return result;
        }

        float fov_estimate = fov_estimate_deg * M_PI / 180.0f;
        int max_stars = std::min((int)centroids.size(), pattern_checking_stars);

        for (int i = 0; i < max_stars - 3; i++)
        {
            for (int j = i + 1; j < max_stars - 2; j++)
            {
                for (int k = j + 1; k < max_stars - 1; k++)
                {
                    for (int l = k + 1; l < max_stars; l++)
                    {
                        std::vector<Centroid> pattern_centroids = {
                            centroids[i], centroids[j], centroids[k], centroids[l]
                        };

                        auto pattern_vectors = compute_vectors(pattern_centroids, height, width, fov_estimate);
                        auto edge_ratios = calculate_edge_ratios(pattern_vectors);
                        auto hash_code_list = generate_hash_code_combinations(
                            edge_ratios, pattern_bins, 0.005f);


                        // Convert to hash bins
                        std::vector<int> hash_bins;
                        for (float ratio : edge_ratios)
                        {
                            int bin = (int)(ratio * pattern_bins);
                            bin = std::max(0, std::min(pattern_bins - 1, bin));
                            hash_bins.push_back(bin);
                        }

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
                                const auto& catalog_pattern = pattern_catalog[catalog_index];

                                // Get catalog star vectors
                                std::vector<std::array<float, 3>> catalog_vectors;
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
                                        // std::cout << "Error is too great! Max error: " << max_error << std::endl;
                                        match = false;
                                        break;
                                    }
                                }

                                if (match)
                                {
                                    std::cout << "Pattern match found! Max error: " << max_error << std::endl;

                                    // Sort both patterns by distance from centroid (if not pre-sorted)
                                    auto sorted_image_vectors = sort_pattern_by_centroid(pattern_vectors);
                                    auto sorted_catalog_vectors = sort_pattern_by_centroid(catalog_vectors);

                                    // Calculate rotation matrix
                                    auto rotation_matrix = find_rotation_matrix(
                                        sorted_image_vectors, sorted_catalog_vectors);

                                    std::cout << "Rotation matrix:" << std::endl;
                                    for (int i = 0; i < 3; i++)
                                    {
                                        for (int j = 0; j < 3; j++)
                                        {
                                            std::cout << rotation_matrix[i][j] << " ";
                                        }
                                        std::cout << std::endl;
                                    }

                                    // Extract RA, Dec, Roll from rotation matrix
                                    // The first column of rotation matrix is the boresight direction
                                    float boresight_x = rotation_matrix[0][0];
                                    float boresight_y = rotation_matrix[0][1];
                                    float boresight_z = rotation_matrix[0][2];

                                    std::cout << "Boresight vector: " << boresight_x << ", " << boresight_y << ", " <<
                                        boresight_z << std::endl;

                                    // Convert to RA/Dec
                                    float ra = std::atan2(boresight_y, boresight_x);
                                    if (ra < 0) ra += 2.0f * M_PI; // Ensure positive

                                    float dec = std::asin(std::max(-1.0f, std::min(1.0f, boresight_z)));
                                    // Clamp for safety

                                    // Calculate roll (rotation around boresight)
                                    float roll = std::atan2(rotation_matrix[1][2], rotation_matrix[2][2]);
                                    if (roll < 0) roll += 2.0f * M_PI;

                                    std::cout << "RA: " << ra << " rad, Dec: " << dec << " rad, Roll: " << roll <<
                                        " rad" <<
                                        std::endl;

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

                                    float fov_refined = fov_estimate;
                                    if (image_largest_edge > 0.001f)
                                    {
                                        // Avoid division by zero
                                        fov_refined = catalog_largest_edge / image_largest_edge * fov_estimate;
                                    }

                                    std::cout << "FOV estimate: " << fov_estimate << ", refined: " << fov_refined <<
                                        std::endl;

                                    result.solved = true;
                                    result.ra = ra * 180.0f / M_PI; // Convert to degrees
                                    result.dec = dec * 180.0f / M_PI; // Convert to degrees
                                    result.roll = roll * 180.0f / M_PI; // Convert to degrees
                                    result.fov = fov_refined * 180.0f / M_PI; // Convert to degrees
                                    result.num_matches = 4; // Could expand this with verification

                                    auto end_time = std::chrono::high_resolution_clock::now();
                                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                                        end_time - start_time);
                                    result.solve_time_ms = duration.count() / 1000.0f;

                                    // return result; // Return successful result
                                }
                            }
                        }
                    }
                }
            }
        }

        // No solution found
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        result.solve_time_ms = duration.count() / 1000.0f;
        return result;
    }

    // Helper function to sort pattern by distance from centroid
    std::vector<std::array<float, 3>> sort_pattern_by_centroid(const std::vector<std::array<float, 3>>& vectors)
    {
        // Calculate centroid
        std::array<float, 3> centroid = {0, 0, 0};
        for (const auto& vec : vectors)
        {
            centroid[0] += vec[0];
            centroid[1] += vec[1];
            centroid[2] += vec[2];
        }
        centroid[0] /= vectors.size();
        centroid[1] /= vectors.size();
        centroid[2] /= vectors.size();

        // Calculate distances and create index array
        std::vector<std::pair<float, int>> distances;
        for (int i = 0; i < vectors.size(); i++)
        {
            float dist = std::sqrt(
                (vectors[i][0] - centroid[0]) * (vectors[i][0] - centroid[0]) +
                (vectors[i][1] - centroid[1]) * (vectors[i][1] - centroid[1]) +
                (vectors[i][2] - centroid[2]) * (vectors[i][2] - centroid[2])
            );
            distances.push_back({dist, i});
        }

        // Sort by distance
        std::sort(distances.begin(), distances.end());

        // Return sorted vectors
        std::vector<std::array<float, 3>> sorted_vectors;
        for (const auto& pair : distances)
        {
            sorted_vectors.push_back(vectors[pair.second]);
        }
        return sorted_vectors;
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
        "/home/linus/git/OpenAstroExplorer/tracker/.venv/lib/python3.12/site-packages/tetra3_db_stars.bin";
    std::string patterns_file =
        "/home/linus/git/OpenAstroExplorer/tracker/.venv/lib/python3.12/site-packages/tetra3_db_patterns.bin";

    std::cout << "Loading tetra3 database: " << database_path << std::endl;

    std::vector<Centroid> centroids = {
        {298.2792, 256.10538},
        {4.6912594, 635.3345},
        {322.22626, 200.69469},
        {43.070637, 219.5646},
        {347.67383, 870.07025},
        {229.61739, 265.70392},
        {265.86844, 581.27826},
        {122.80346, 216.62479},
        {510.49747, 690.5212},
        {455.73697, 774.504},
        {493.0147, 248.57089},
        {550.51825, 575.48016},
        {115.60284, 140.35214},
        {704.5442, 726.65936},
        {680.7061, 733.4033},
        {160.39578, 470.6078},
        {568.74475, 373.4044},
        {538.50354, 397.46884},
        {365.4988, 468.51767},
        {228.63228, 402.7744},
        {132.41939, 280.6431},
        {288.55194, 181.4497},
        {442.42972, 366.51596}
    };

    if (loader.load_database(stars_file, patterns_file))
    {
        std::cout << "Database loaded successfully!" << std::endl;

        // Now you can use these with your star solver
        SimpleStarSolver solver;
        solver.load_star_catalog(loader.get_stars());
        solver.load_pattern_catalog(loader.get_patterns());
        solver.print_memory_usage();
        auto result = solver.solve_from_centroids(centroids, 768, 1064, 15.0f, 8);

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
}
