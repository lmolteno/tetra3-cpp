#include "star_detector.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <stack>

StarDetector::StarDetector(const StarDetectorConfig &config) : config(config) {}

StarDetector::BackgroundStats StarDetector::estimate_background(const uint16_t *data, int width, int height) {
    // Subsample every 4th pixel for speed
    std::vector<double> samples;
    samples.reserve((width / 4) * (height / 4));

    for (int y = 0; y < height; y += 4) {
        for (int x = 0; x < width; x += 4) {
            samples.push_back(static_cast<double>(data[y * width + x]));
        }
    }

    double mean = 0, stddev = 0;

    // Iterative sigma-clipping
    for (int iter = 0; iter < config.sigma_clip_iterations; iter++) {
        if (samples.empty()) break;

        double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
        mean = sum / samples.size();

        double sq_sum = 0;
        for (double v : samples) {
            double diff = v - mean;
            sq_sum += diff * diff;
        }
        stddev = std::sqrt(sq_sum / samples.size());

        if (stddev < 1e-6) break;

        // Reject outliers
        double low = mean - config.sigma_clip_factor * stddev;
        double high = mean + config.sigma_clip_factor * stddev;

        std::vector<double> filtered;
        filtered.reserve(samples.size());
        for (double v : samples) {
            if (v >= low && v <= high) {
                filtered.push_back(v);
            }
        }
        samples = std::move(filtered);
    }

    // Final stats
    if (!samples.empty()) {
        double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
        mean = sum / samples.size();
        double sq_sum = 0;
        for (double v : samples) {
            double diff = v - mean;
            sq_sum += diff * diff;
        }
        stddev = std::sqrt(sq_sum / samples.size());
    }

    return {mean, stddev};
}

int StarDetector::flood_fill(const std::vector<uint8_t> &mask, int width, int height,
                              int start_y, int start_x, std::vector<uint8_t> &visited,
                              std::vector<Point> &cluster_points) {
    std::stack<Point> stack;
    stack.push({start_y, start_x});
    int count = 0;
    cluster_points.clear();

    while (!stack.empty()) {
        auto [y, x] = stack.top();
        stack.pop();

        if (y < 0 || y >= height || x < 0 || x >= width) continue;
        int idx = y * width + x;
        if (visited[idx] || !mask[idx]) continue;

        visited[idx] = 1;
        count++;

        if (count <= config.max_cluster_size) {
            cluster_points.push_back({y, x});
        }

        // 8-connected neighbors
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dy == 0 && dx == 0) continue;
                stack.push({y + dy, x + dx});
            }
        }
    }

    return count;
}

std::vector<DetectedStar> StarDetector::detect(const Frame &frame) {
    const uint16_t *data = frame.data.data();
    int width = frame.width;
    int height = frame.height;

    // Step 1: Background estimation via sigma-clipping
    auto [bg_mean, bg_stddev] = estimate_background(data, width, height);
    std::cout << "Background: mean=" << std::fixed << std::setprecision(1) << bg_mean
              << ", stddev=" << std::setprecision(1) << bg_stddev << std::endl;

    // Step 2: Create binary mask at detection threshold
    double threshold = bg_mean + config.detection_sigma * bg_stddev;
    std::vector<uint8_t> mask(width * height, 0);
    for (int i = 0; i < width * height; i++) {
        if (data[i] > threshold) {
            mask[i] = 1;
        }
    }

    // Step 3: Flood fill to find clusters
    std::vector<uint8_t> visited(width * height, 0);
    std::vector<DetectedStar> stars;
    std::vector<Point> cluster_points;

    for (int y = 0; y < height && static_cast<int>(stars.size()) < config.max_stars; y++) {
        for (int x = 0; x < width && static_cast<int>(stars.size()) < config.max_stars; x++) {
            int idx = y * width + x;
            if (mask[idx] && !visited[idx]) {
                int cluster_size = flood_fill(mask, width, height, y, x, visited, cluster_points);

                if (cluster_size < config.min_cluster_size || cluster_size > config.max_cluster_size) {
                    continue;
                }

                // Step 4: Intensity-weighted centroid
                double sum_iy = 0, sum_ix = 0, sum_i = 0;
                double peak = 0;

                for (const auto &pt : cluster_points) {
                    double intensity = static_cast<double>(data[pt.y * width + pt.x]) - bg_mean;
                    if (intensity < 0) intensity = 0;
                    sum_iy += pt.y * intensity;
                    sum_ix += pt.x * intensity;
                    sum_i += intensity;
                    if (data[pt.y * width + pt.x] > peak) {
                        peak = data[pt.y * width + pt.x];
                    }
                }

                if (sum_i < 1e-6) continue;

                double cy = sum_iy / sum_i;
                double cx = sum_ix / sum_i;

                // Step 5: SNR calculation (Poisson + read noise model)
                // SNR = total_flux / sqrt(total_flux + n_pixels * sigma_bg^2)
                double total_flux = sum_i;
                int n_pixels = static_cast<int>(cluster_points.size());
                double noise = std::sqrt(total_flux + n_pixels * bg_stddev * bg_stddev);
                double snr = (noise > 0) ? total_flux / noise : 0;

                DetectedStar star;
                star.centroid = {cy, cx};
                star.snr = snr;
                star.peak_value = peak;
                star.flux = total_flux;
                star.cluster_size = cluster_size;
                stars.push_back(star);
            }
        }
    }

    // Sort by descending SNR
    std::sort(stars.begin(), stars.end(),
              [](const DetectedStar &a, const DetectedStar &b) {
                  return a.snr > b.snr;
              });

    // Log each star
    std::cout << "Detected " << stars.size() << " stars (threshold=" << std::fixed
              << std::setprecision(0) << threshold << ")" << std::endl;
    for (size_t i = 0; i < stars.size(); i++) {
        const auto &s = stars[i];
        std::cout << "  Star " << std::setw(3) << i
                  << ": pos=(" << std::fixed << std::setprecision(1) << s.centroid.y
                  << ", " << std::setprecision(1) << s.centroid.x << ")"
                  << " SNR=" << std::setprecision(1) << s.snr
                  << " peak=" << std::setprecision(0) << s.peak_value
                  << " flux=" << std::setprecision(0) << s.flux
                  << " size=" << s.cluster_size << std::endl;
    }

    return stars;
}
