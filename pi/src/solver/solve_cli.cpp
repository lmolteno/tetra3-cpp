// CLI tool: reads a 16-bit grayscale image, detects stars, solves, outputs JSON.
// Two-pass approach: wide solve on full image, then tight solve on center crop.
// Usage: solve_cli [options] <image>
//        solve_cli [options] --batch   (reads image paths from stdin, one per line)
// Output: one JSON line per image

#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <chrono>

#include <tetra3/solver.h>
#include <tetra3/database_loader.h>

#include "capture/frame_source.h"
#include "capture/image_loader.h"
#include "solver/star_detector.h"

static void solve_image(const std::string &image_path,
                        SimpleStarSolver &solver,
                        float fov, float detection_sigma, int crop_size,
                        int bg_tile_size, int max_cluster_size,
                        float min_sharpness, bool debug,
                        bool emit_matches,
                        std::optional<float> distortion_in,
                        std::optional<float> freeze_k,
                        float extended_match_radius,
                        const std::vector<double> &poly_undist) {
    using clock = std::chrono::steady_clock;
    auto ms_since = [](clock::time_point t) {
        return std::chrono::duration<double, std::milli>(clock::now() - t).count();
    };

    auto t_load = clock::now();
    ImageFileSource source(image_path, debug);
    if (!source.initialize()) {
        printf("{\"solved\":false,\"error\":\"image_load_failed\"}\n");
        fflush(stdout);
        return;
    }

    auto frame_opt = source.capture();
    if (!frame_opt) {
        printf("{\"solved\":false,\"error\":\"capture_failed\"}\n");
        fflush(stdout);
        return;
    }
    double load_ms = ms_since(t_load);

    auto &frame = *frame_opt;

    StarDetectorConfig det_cfg;
    det_cfg.detection_sigma = detection_sigma;
    det_cfg.bg_tile_size = bg_tile_size;
    det_cfg.max_cluster_size = max_cluster_size;
    det_cfg.min_sharpness = min_sharpness;
    det_cfg.verbose = debug;
    StarDetector detector(det_cfg);

    auto t_detect = clock::now();
    StarDetectorStats det_stats_full;
    auto detected = detector.detect(frame, &det_stats_full);
    double detect_ms = ms_since(t_detect);
    StarDetectorStats det_stats_crop;

    if (detected.empty()) {
        printf("{\"solved\":false,\"error\":\"no_stars\",\"num_detected\":0,"
               "\"load_ms\":%.1f,\"detect_ms\":%.1f,"
               "\"bg_ms\":%.1f,\"mask_ms\":%.1f,\"cluster_ms\":%.1f}\n",
               load_ms, detect_ms,
               det_stats_full.bg_ms, det_stats_full.mask_ms, det_stats_full.cluster_ms);
        fflush(stdout);
        return;
    }

    int num_detected_full = static_cast<int>(detected.size());

    // Lens calibration: poly_undist = [cx_off, cy_off, k1, k2, k3, ...].
    // We let the solver apply it internally (set_lens_calibration). For the
    // crop pass, k coefficients must be rescaled because the solver
    // normalizes radius by its own half_width.
    auto build_k = [&](float scale_squared) {
        std::vector<float> ks;
        if (poly_undist.size() <= 2) return ks;
        ks.reserve(poly_undist.size() - 2);
        float p = scale_squared;
        for (size_t i = 2; i < poly_undist.size(); ++i) {
            ks.push_back(static_cast<float>(poly_undist[i]) * p);
            p *= scale_squared;
        }
        return ks;
    };
    float lens_cx_off = poly_undist.size() >= 2 ? static_cast<float>(poly_undist[0]) : 0.0f;
    float lens_cy_off = poly_undist.size() >= 2 ? static_cast<float>(poly_undist[1]) : 0.0f;

    // --- Pass 1: tight solve on center crop ---
    int cw = std::min(crop_size, static_cast<int>(frame.width));
    int ch = std::min(crop_size, static_cast<int>(frame.height));
    int x0 = (frame.width - cw) / 2;
    int y0 = (frame.height - ch) / 2;

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
    std::vector<MatchInfo> matches;
    std::vector<MatchInfo> *matches_arg = emit_matches ? &matches : nullptr;
    bool used_crop = false;

    auto t_crop_detect = clock::now();
    auto crop_detected = detector.detect(crop_frame, &det_stats_crop);
    detect_ms += ms_since(t_crop_detect);
    if (crop_detected.size() >= 4) {
        std::vector<Centroid> crop_centroids;
        crop_centroids.reserve(crop_detected.size());
        for (const auto &s : crop_detected) {
            crop_centroids.push_back(s.centroid);
        }

        // The solver's distortion math normalizes radius by half-width, so
        // a single k or a polynomial calibrated against the full image must
        // be rescaled when the solver runs on the crop. For the crop:
        //   r_undist = r - k * r^3 / half_w²    (single-k legacy)
        //   r_norm   = r_pixels / (cw/2)        (polynomial argument)
        // so each polynomial coefficient k_i scales by (cw/W)^(2i).
        std::optional<float> dist_for_crop;
        if (freeze_k.has_value()) {
            float scale = static_cast<float>(cw) / static_cast<float>(frame.width);
            dist_for_crop = freeze_k.value() * scale * scale;
        } else {
            dist_for_crop = distortion_in;
        }

        // Set lens calibration on the solver (or clear it) for the crop pass.
        if (poly_undist.size() >= 3) {
            float ratio = static_cast<float>(cw) / static_cast<float>(frame.width);
            solver.set_lens_calibration(lens_cx_off, lens_cy_off,
                                         build_k(ratio * ratio));
        } else {
            solver.set_lens_calibration(0.0f, 0.0f, {});
        }

        result = solver.solve_from_centroids(
            crop_centroids, ch, cw, crop_fov_deg, 16,
            0.01f, 1e-6f, dist_for_crop, crop_fov_deg * 0.2f, matches_arg,
            /*freeze_distortion=*/freeze_k.has_value(),
            /*extended_match_radius=*/extended_match_radius);

        float max_crop_rmse = crop_fov_deg * 3600.0f / cw;
        if (result.solved && result.num_matches >= 6 && result.rmse < max_crop_rmse) {
            result.fov = 2.0f * std::atan(
                std::tan(result.fov / 2.0f) * frame.width / cw);
            used_crop = true;
        } else {
            result.solved = false;
            if (matches_arg) matches_arg->clear();
        }
    }

    // --- Pass 2: fall back to full image if crop failed ---
    if (!result.solved) {
        std::vector<Centroid> centroids;
        centroids.reserve(detected.size());
        for (const auto &s : detected) {
            centroids.push_back(s.centroid);
        }

        // Set lens calibration for the full-image solve (no rescaling).
        if (poly_undist.size() >= 3) {
            solver.set_lens_calibration(lens_cx_off, lens_cy_off, build_k(1.0f));
        } else {
            solver.set_lens_calibration(0.0f, 0.0f, {});
        }

        std::optional<float> dist_for_full = freeze_k.has_value() ? freeze_k : distortion_in;
        result = solver.solve_from_centroids(
            centroids, frame.height, frame.width, fov, 16,
            0.01f, 0.001f, dist_for_full, fov * 0.2f, matches_arg,
            /*freeze_distortion=*/freeze_k.has_value(),
            /*extended_match_radius=*/extended_match_radius);
    }

    // If the crop pass produced the solve, remap match coordinates from the
    // crop's frame into the full image. The center crop is just a translate;
    // focal length is identical because the crop FOV was derived from it.
    if (emit_matches && used_crop) {
        for (auto &m : matches) {
            m.img_x += x0;  m.img_y += y0;
            m.pred_x += x0; m.pred_y += y0;
        }
    }

    // Sum detector sub-stages across the full + crop passes so the breakdown
    // covers the same wall clock that detect_ms measures.
    double bg_ms      = det_stats_full.bg_ms      + det_stats_crop.bg_ms;
    double mask_ms    = det_stats_full.mask_ms    + det_stats_crop.mask_ms;
    double cluster_ms = det_stats_full.cluster_ms + det_stats_crop.cluster_ms;

    // Output JSON. solve_time_ms is the matcher's internal time (just pattern
    // search); load_ms / detect_ms cover the rest of solve_cli's wall clock.
    if (result.solved) {
        printf("{\"solved\":true,"
               "\"ra_rad\":%.8f,\"dec_rad\":%.8f,\"roll_rad\":%.8f,"
               "\"fov_rad\":%.8f,\"rmse\":%.4f,"
               "\"num_matches\":%d,\"solve_time_ms\":%.1f,"
               "\"load_ms\":%.1f,\"detect_ms\":%.1f,"
               "\"bg_ms\":%.1f,\"mask_ms\":%.1f,\"cluster_ms\":%.1f,"
               "\"num_detected\":%d,"
               "\"ra_deg\":%.6f,\"dec_deg\":%.6f,\"fov_deg\":%.6f,"
               "\"used_crop\":%s,\"width\":%d,\"height\":%d,\"distortion_k\":%.6f",
               result.ra, result.dec, result.roll,
               result.fov, result.rmse,
               result.num_matches, result.solve_time_ms,
               load_ms, detect_ms,
               bg_ms, mask_ms, cluster_ms,
               num_detected_full,
               result.ra * 180.0 / M_PI,
               result.dec * 180.0 / M_PI,
               result.fov * 180.0 / M_PI,
               used_crop ? "true" : "false",
               static_cast<int>(frame.width),
               static_cast<int>(frame.height),
               result.distortion_k);
        if (emit_matches) {
            printf(",\"matches\":[");
            for (size_t i = 0; i < matches.size(); ++i) {
                const auto &m = matches[i];
                printf("%s{\"cat_idx\":%d,\"ra\":%.8f,\"dec\":%.8f,\"mag\":%.3f,"
                       "\"img_x\":%.4f,\"img_y\":%.4f,"
                       "\"pred_x\":%.4f,\"pred_y\":%.4f,"
                       "\"residual_pix\":%.4f}",
                       i ? "," : "",
                       m.catalog_index, m.catalog_ra, m.catalog_dec, m.magnitude,
                       m.img_x, m.img_y, m.pred_x, m.pred_y, m.residual_pix);
            }
            printf("]");
        }
        printf("}\n");
    } else {
        printf("{\"solved\":false,\"solve_time_ms\":%.1f,"
               "\"load_ms\":%.1f,\"detect_ms\":%.1f,"
               "\"bg_ms\":%.1f,\"mask_ms\":%.1f,\"cluster_ms\":%.1f,"
               "\"num_detected\":%d}\n",
               result.solve_time_ms, load_ms, detect_ms,
               bg_ms, mask_ms, cluster_ms, num_detected_full);
    }
    fflush(stdout);
}

int main(int argc, char *argv[]) {
    std::string image_path;
    std::string db_stars = "../tetra3_db_stars.bin";
    std::string db_patterns = "../tetra3_db_patterns.bin";
    float fov = 70.0f;
    float detection_sigma = 5.0f;
    int crop_size = 720;
    int bg_tile_size = 128;
    int max_cluster_size = 30;
    float min_sharpness = 2.0f;
    bool batch_mode = false;
    bool debug = false;
    bool emit_matches = false;
    std::optional<float> distortion_in = std::nullopt;
    std::optional<float> freeze_k = std::nullopt;
    float extended_match_radius = 0.0f;
    std::vector<double> poly_undist;  // [cx_off, cy_off, k1, k2, k3, ...]

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
        } else if (arg == "--bg-tile" && i + 1 < argc) {
            bg_tile_size = std::stoi(argv[++i]);
        } else if (arg == "--max-cluster" && i + 1 < argc) {
            max_cluster_size = std::stoi(argv[++i]);
        } else if (arg == "--min-sharpness" && i + 1 < argc) {
            min_sharpness = std::stof(argv[++i]);
        } else if (arg == "--batch") {
            batch_mode = true;
        } else if (arg == "--emit-matches") {
            emit_matches = true;
        } else if (arg == "--distortion-fit") {
            // Enable solver's per-frame distortion fit (initial guess 0).
            distortion_in = 0.0f;
        } else if (arg == "--distortion-k" && i + 1 < argc) {
            // Pre-set a fixed distortion coefficient (also triggers
            // the solver's LS fit, which overwrites with its own value).
            distortion_in = std::stof(argv[++i]);
        } else if (arg == "--freeze-distortion-k" && i + 1 < argc) {
            // Undistort centroids in solve_cli using this k, then call the
            // solver in no-distortion mode so it doesn't refit f and k.
            freeze_k = std::stof(argv[++i]);
        } else if (arg == "--extended-match-radius" && i + 1 < argc) {
            // After solving, do an unrestricted projection+match pass with
            // this fractional radius (e.g. 0.02 = 2% of width) to harvest
            // many more pairs for distortion analysis. Replaces matches_out.
            extended_match_radius = std::stof(argv[++i]);
        } else if (arg == "--poly-undistort" && i + 1 < argc) {
            // Comma-separated list "cx_off,cy_off,k1,k2,k3,...". The model
            // is dx = (x-cx_t)*F(r²) where F = sum_i k_i * r^(2i).
            // Applied to all detections (in full-image coords) before solving.
            std::string s = argv[++i];
            size_t pos = 0;
            while (pos < s.size()) {
                size_t comma = s.find(',', pos);
                std::string tok = s.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
                if (!tok.empty()) poly_undist.push_back(std::stod(tok));
                if (comma == std::string::npos) break;
                pos = comma + 1;
            }
        } else if (arg == "--debug") {
            debug = true;
        } else if (arg[0] != '-') {
            image_path = arg;
        }
    }

    if (!batch_mode && image_path.empty()) {
        std::cerr << "Usage: solve_cli [options] <image>" << std::endl;
        std::cerr << "       solve_cli [options] --batch  (read paths from stdin)" << std::endl;
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

    if (batch_mode) {
        std::string line;
        while (std::getline(std::cin, line)) {
            // Trim whitespace
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r' || line.back() == ' '))
                line.pop_back();
            if (line.empty()) continue;
            solve_image(line, solver, fov, detection_sigma, crop_size, bg_tile_size,
                        max_cluster_size, min_sharpness, debug, emit_matches,
                        distortion_in, freeze_k, extended_match_radius,
                        poly_undist);
        }
    } else {
        solve_image(image_path, solver, fov, detection_sigma, crop_size, bg_tile_size,
                    max_cluster_size, min_sharpness, debug, emit_matches,
                    distortion_in, freeze_k, extended_match_radius,
                    poly_undist);
    }

    return 0;
}
