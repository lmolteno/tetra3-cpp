#include "solver/star_detector.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <stack>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

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

StarDetector::TiledBackground StarDetector::compute_local_background(
        const uint16_t *data, int width, int height) {
    int tile = config.bg_tile_size;
    int tiles_x = (width + tile - 1) / tile;
    int tiles_y = (height + tile - 1) / tile;

    TiledBackground bg;
    bg.tiles_x = tiles_x;
    bg.tile_size = tile;
    bg.tile_mean.resize(tiles_x * tiles_y);
    bg.tile_std.resize(tiles_x * tiles_y);
    std::vector<double> samples;

    for (int ty = 0; ty < tiles_y; ty++) {
        for (int tx = 0; tx < tiles_x; tx++) {
            int y0 = ty * tile, y1 = std::min(y0 + tile, height);
            int x0 = tx * tile, x1 = std::min(x0 + tile, width);

            samples.clear();
            for (int y = y0; y < y1; y += 2)
                for (int x = x0; x < x1; x += 2)
                    samples.push_back(data[y * width + x]);

            double mean = 0, stddev = 0;
            for (int iter = 0; iter < config.sigma_clip_iterations && !samples.empty(); iter++) {
                double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
                mean = sum / samples.size();
                double sq_sum = 0;
                for (double v : samples) { double d = v - mean; sq_sum += d * d; }
                stddev = std::sqrt(sq_sum / samples.size());
                if (stddev < 1e-6) break;

                double lo = mean - config.sigma_clip_factor * stddev;
                double hi = mean + config.sigma_clip_factor * stddev;
                std::erase_if(samples, [lo, hi](double v) { return v < lo || v > hi; });
            }
            if (!samples.empty()) {
                double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
                mean = sum / samples.size();
                double sq_sum = 0;
                for (double v : samples) { double d = v - mean; sq_sum += d * d; }
                stddev = std::sqrt(sq_sum / samples.size());
            }

            bg.tile_mean[ty * tiles_x + tx] = static_cast<float>(mean);
            bg.tile_std[ty * tiles_x + tx] = static_cast<float>(stddev);
        }
    }

    return bg;
}

std::vector<DetectedStar> StarDetector::detect(const Frame &frame,
                                                StarDetectorStats *stats) {
    const uint16_t *data = frame.data.data();
    int width = frame.width;
    int height = frame.height;

    using clock = std::chrono::steady_clock;
    auto ms_since = [](clock::time_point t) {
        return std::chrono::duration<double, std::milli>(clock::now() - t).count();
    };

    // Step 1: Local background estimation (tile grid only, no per-pixel map)
    auto t_bg = clock::now();
    auto bg = compute_local_background(data, width, height);
    double bg_ms = ms_since(t_bg);

    if (config.verbose) {
        auto [bg_mean, bg_stddev] = estimate_background(data, width, height);
        std::cout << "Background: mean=" << std::fixed << std::setprecision(1) << bg_mean
                  << ", stddev=" << std::setprecision(1) << bg_stddev << std::endl;
    }

    // Step 2: Create binary mask using local threshold.
    //
    // The bilinear thresh(x, y) is piecewise-linear in x once y is fixed —
    // within a span where (tx0, tx1) don't change it reduces to
    //   thresh(x) = thresh_at_x + delta * (x - x_start)
    // both terms constant for the span. We walk each row in tile-aligned
    // spans and let NEON fill 8 pixels per iter via running-add. The leading
    // [0, half) and trailing partial spans use delta=0 (clamped fx).
    auto t_mask = clock::now();
    int tile = bg.tile_size;
    int ntiles = static_cast<int>(bg.tile_mean.size());
    int tiles_x = bg.tiles_x;
    int tiles_y = ntiles / tiles_x;
    float half = tile * 0.5f;
    float inv_tile = 1.0f / tile;

    std::vector<float> tile_thresh(ntiles);
    for (int i = 0; i < ntiles; i++)
        tile_thresh[i] = bg.tile_mean[i] + config.detection_sigma * bg.tile_std[i];

    std::vector<uint8_t> mask(width * height, 0);
    for (int y = 0; y < height; y++) {
        float fy_raw = (static_cast<float>(y) - half) * inv_tile;
        int ty0 = std::max(0, static_cast<int>(fy_raw));
        int ty1 = std::min(tiles_y - 1, ty0 + 1);
        float fy = fy_raw - ty0;
        if (fy < 0.0f) fy = 0.0f;
        int row0 = ty0 * tiles_x;
        int row1 = ty1 * tiles_x;
        float w0y = 1.0f - fy, w1y = fy;

        const uint16_t *row_data = data + static_cast<size_t>(y) * width;
        uint8_t *row_mask = mask.data() + static_cast<size_t>(y) * width;

        int x = 0;
        while (x < width) {
            float fx_raw = (static_cast<float>(x) - half) * inv_tile;
            int tx0 = std::max(0, static_cast<int>(fx_raw));
            int tx1 = std::min(tiles_x - 1, tx0 + 1);
            bool left_clamp  = (fx_raw < 0.0f);            // x in [0, half)
            bool right_clamp = (tx0 == tx1);               // last partial span

            int span_end;
            if (left_clamp) {
                span_end = std::min(static_cast<int>(half), width);
            } else if (right_clamp) {
                span_end = width;
            } else {
                span_end = std::min(static_cast<int>((tx0 + 1) * tile + half), width);
            }

            float a = tile_thresh[row0 + tx0];
            float b = tile_thresh[row0 + tx1];
            float c = tile_thresh[row1 + tx0];
            float d = tile_thresh[row1 + tx1];

            float fx_at_x = left_clamp ? 0.0f : (fx_raw - tx0);
            float thresh_at_x = a*(1.0f - fx_at_x)*w0y + b*fx_at_x*w0y
                              + c*(1.0f - fx_at_x)*w1y + d*fx_at_x*w1y;
            float delta = (left_clamp || right_clamp)
                ? 0.0f
                : ((b - a) * w0y + (d - c) * w1y) * inv_tile;

#if defined(__ARM_NEON)
            if (x + 8 <= span_end) {
                static const float lane_off[8] = {0, 1, 2, 3, 4, 5, 6, 7};
                float32x4_t base    = vdupq_n_f32(thresh_at_x);
                float32x4_t delta_v = vdupq_n_f32(delta);
                float32x4_t step    = vdupq_n_f32(delta * 8.0f);
                float32x4_t off_lo  = vld1q_f32(lane_off);
                float32x4_t off_hi  = vld1q_f32(lane_off + 4);
                float32x4_t thresh_lo = vmlaq_f32(base, delta_v, off_lo);
                float32x4_t thresh_hi = vmlaq_f32(base, delta_v, off_hi);

                const uint8x8_t one = vdup_n_u8(1);
                for (; x + 8 <= span_end; x += 8) {
                    uint16x8_t  d_v  = vld1q_u16(row_data + x);
                    float32x4_t d_lo = vcvtq_f32_u32(vmovl_u16(vget_low_u16(d_v)));
                    float32x4_t d_hi = vcvtq_f32_u32(vmovl_u16(vget_high_u16(d_v)));

                    uint32x4_t cmp_lo = vcgtq_f32(d_lo, thresh_lo);
                    uint32x4_t cmp_hi = vcgtq_f32(d_hi, thresh_hi);
                    uint16x8_t cmp    = vcombine_u16(vmovn_u32(cmp_lo), vmovn_u32(cmp_hi));
                    uint8x8_t  out    = vand_u8(vmovn_u16(cmp), one);
                    vst1_u8(row_mask + x, out);

                    thresh_lo = vaddq_f32(thresh_lo, step);
                    thresh_hi = vaddq_f32(thresh_hi, step);
                }
            }
#endif
            for (; x < span_end; x++) {
                float fx_pix = left_clamp ? 0.0f
                                          : ((static_cast<float>(x) - half) * inv_tile - tx0);
                float thresh = a*(1.0f - fx_pix)*w0y + b*fx_pix*w0y
                             + c*(1.0f - fx_pix)*w1y + d*fx_pix*w1y;
                if (row_data[x] > thresh) row_mask[x] = 1;
            }
        }
    }

    double mask_ms = ms_since(t_mask);

    // Step 3: Flood fill to find clusters
    auto t_cluster = clock::now();
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

                // Reject elongated clusters (edges/structure) via bounding box aspect ratio
                int min_y = height, max_y = 0, min_x = width, max_x = 0;
                for (const auto &pt : cluster_points) {
                    if (pt.y < min_y) min_y = pt.y;
                    if (pt.y > max_y) max_y = pt.y;
                    if (pt.x < min_x) min_x = pt.x;
                    if (pt.x > max_x) max_x = pt.x;
                }
                int bbox_w = max_x - min_x + 1;
                int bbox_h = max_y - min_y + 1;
                float aspect = static_cast<float>(std::max(bbox_w, bbox_h)) /
                               std::max(1, std::min(bbox_w, bbox_h));
                if (aspect > 3.0f) continue;  // stars are roughly circular

                // Step 4: Intensity-weighted centroid using local background
                double sum_iy = 0, sum_ix = 0, sum_i = 0;
                double peak = 0;
                double local_bg_std_sum = 0;

                for (const auto &pt : cluster_points) {
                    int pidx = pt.y * width + pt.x;
                    double intensity = static_cast<double>(data[pidx]) - bg.mean_at(pt.y, pt.x);
                    if (intensity < 0) intensity = 0;
                    sum_iy += pt.y * intensity;
                    sum_ix += pt.x * intensity;
                    sum_i += intensity;
                    local_bg_std_sum += bg.std_at(pt.y, pt.x);
                    if (data[pidx] > peak) {
                        peak = data[pidx];
                    }
                }

                if (sum_i < 1e-6) continue;

                double cy = sum_iy / sum_i;
                double cx = sum_ix / sum_i;

                // Step 5: SNR using local background stddev
                double total_flux = sum_i;
                int n_pixels = static_cast<int>(cluster_points.size());
                double avg_bg_std = local_bg_std_sum / n_pixels;
                double noise = std::sqrt(total_flux + n_pixels * avg_bg_std * avg_bg_std);
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
    double cluster_ms = ms_since(t_cluster);

    if (stats) {
        stats->bg_ms      = bg_ms;
        stats->mask_ms    = mask_ms;
        stats->cluster_ms = cluster_ms;
    }

    if (config.verbose) {
        std::cout << "Detected " << stars.size() << " stars (local background, "
                  << config.detection_sigma << " sigma)" << std::endl;
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
    }

    return stars;
}
