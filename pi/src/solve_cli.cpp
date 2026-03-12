// CLI tool: reads a 16-bit grayscale image, detects stars, solves, outputs JSON.
// Two-pass approach: wide solve on full image, then tight solve on center crop.
// Usage: solve_cli --db-stars <path> --db-patterns <path> --fov <deg> <image>
// Output: single JSON line with solve result

#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

#include <tetra3/solver.h>
#include <tetra3/database_loader.h>

#include "frame_source.h"
#include "image_loader.h"
#include "star_detector.h"

int main(int argc, char *argv[]) {
    std::string image_path;
    std::string db_stars = "../tetra3_db_stars.bin";
    std::string db_patterns = "../tetra3_db_patterns.bin";
    float fov = 11.0f;
    float detection_sigma = 5.0f;
    int crop_size = 720;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--db-stars" && i + 1 < argc) {
            db_stars = argv[++i];
        } else if (arg == "--db-patterns" && i + 1 < argc) {
            db_patterns = argv[++i];
        } else if (arg == "--fov" && i + 1 < argc) {
            fov = std::stof(argv[++i]);
        } else if (arg == "--detection-sigma" && i + 1 < argc) {
            detection_sigma = std::stof(argv[++i]);
        } else if (arg == "--crop" && i + 1 < argc) {
            crop_size = std::stoi(argv[++i]);
        } else if (arg[0] != '-') {
            image_path = arg;
        }
    }

    if (image_path.empty()) {
        std::cerr << "Usage: solve_cli [options] <image>" << std::endl;
        return 1;
    }

    BinaryDatabaseLoader loader;
    if (!loader.load_database(db_stars, db_patterns)) {
        fprintf(stderr, "Failed to load database\n");
        printf("{\"solved\":false,\"error\":\"db_load_failed\"}\n");
        return 1;
    }

    SimpleStarSolver solver;
    solver.load_star_catalog(loader.get_stars());
    solver.load_pattern_catalog(loader.get_patterns());

    ImageFileSource source(image_path);
    if (!source.initialize()) {
        printf("{\"solved\":false,\"error\":\"image_load_failed\"}\n");
        return 1;
    }

    auto frame_opt = source.capture();
    if (!frame_opt) {
        printf("{\"solved\":false,\"error\":\"capture_failed\"}\n");
        return 1;
    }

    auto &frame = *frame_opt;

    StarDetectorConfig det_cfg;
    det_cfg.detection_sigma = detection_sigma;
    StarDetector detector(det_cfg);
    auto detected = detector.detect(frame);

    if (detected.empty()) {
        printf("{\"solved\":false,\"error\":\"no_stars\",\"num_detected\":0}\n");
        return 0;
    }

    int num_detected_full = static_cast<int>(detected.size());

    // --- Pass 1: tight solve on center crop ---
    int cw = std::min(crop_size, static_cast<int>(frame.width));
    int ch = std::min(crop_size, static_cast<int>(frame.height));
    int x0 = (frame.width - cw) / 2;
    int y0 = (frame.height - ch) / 2;

    // Crop FOV from the estimated full-image FOV
    float fov_rad = fov * static_cast<float>(M_PI) / 180.0f;
    float focal_px = frame.width / (2.0f * std::tan(fov_rad / 2.0f));
    float crop_fov_deg = 2.0f * std::atan(cw / (2.0f * focal_px)) * 180.0f / static_cast<float>(M_PI);

    Frame crop_frame;
    crop_frame.width = cw;
    crop_frame.height = ch;
    crop_frame.bit_depth = frame.bit_depth;
    crop_frame.data.resize(cw * ch);
    for (int row = 0; row < ch; row++) {
        for (int col = 0; col < cw; col++) {
            crop_frame.data[row * cw + col] = frame.data[(y0 + row) * frame.width + (x0 + col)];
        }
    }

    SolveResult result = {false, 0, 0, 0, 0, 0, 0, 0, 0.0f};

    auto crop_detected = detector.detect(crop_frame);
    if (crop_detected.size() >= 4) {
        std::vector<Centroid> crop_centroids;
        crop_centroids.reserve(crop_detected.size());
        for (const auto &s : crop_detected) {
            crop_centroids.push_back(s.centroid);
        }

        result = solver.solve_from_centroids(
            crop_centroids, ch, cw, crop_fov_deg, 16,
            0.01f, 1e-6f, std::nullopt, crop_fov_deg * 0.2f);

        // Accept crop solve if enough matches AND sub-pixel RMSE.
        // Max RMSE: 1 pixel worth of arcsec = crop_fov / crop_width * 3600.
        float max_crop_rmse = crop_fov_deg * 3600.0f / cw;
        if (result.solved && result.num_matches >= 6 && result.rmse < max_crop_rmse) {
            // Scale crop FOV back to full-image FOV
            result.fov = 2.0f * std::atan(
                std::tan(result.fov / 2.0f) * frame.width / cw);
        } else {
            result.solved = false;
        }
    }

    // --- Pass 2: fall back to full image if crop failed ---
    if (!result.solved) {
        std::vector<Centroid> centroids;
        centroids.reserve(detected.size());
        for (const auto &s : detected) {
            centroids.push_back(s.centroid);
        }

        result = solver.solve_from_centroids(
            centroids, frame.height, frame.width, fov, 16,
            0.01f, 0.001f, std::nullopt, fov * 0.2f);
    }

    // Output JSON
    if (result.solved) {
        printf("{\"solved\":true,"
               "\"ra_rad\":%.8f,\"dec_rad\":%.8f,\"roll_rad\":%.8f,"
               "\"fov_rad\":%.8f,\"rmse\":%.4f,"
               "\"num_matches\":%d,\"solve_time_ms\":%.1f,"
               "\"num_detected\":%d,"
               "\"ra_deg\":%.6f,\"dec_deg\":%.6f,\"fov_deg\":%.6f}\n",
               result.ra, result.dec, result.roll,
               result.fov, result.rmse,
               result.num_matches, result.solve_time_ms,
               num_detected_full,
               result.ra * 180.0 / M_PI,
               result.dec * 180.0 / M_PI,
               result.fov * 180.0 / M_PI);
    } else {
        printf("{\"solved\":false,\"solve_time_ms\":%.1f,\"num_detected\":%d}\n",
               result.solve_time_ms, num_detected_full);
    }

    return 0;
}
