#include <iostream>
#include <string>
#include <memory>
#include <cstring>
#include <cmath>
#include <chrono>
#include <cstdlib>

#include <tetra3/solver.h>
#include <tetra3/database_loader.h>

#include "frame_source.h"
#include "image_loader.h"
#include "star_detector.h"
#include "display.h"
#include "display_null.h"
#include "star_map.h"
#include "polar_align.h"

#ifdef HAS_SDL_DISPLAY
#include "display_sdl.h"
#include <SDL2/SDL.h>
#endif

#ifdef HAS_LIBCAMERA
#include "camera_capture.h"
#endif

#include "display_i2c.h"
#include "display_buffer.h"

struct AppConfig {
    std::string image_path;
    int width = 2028;
    int height = 1520;
    std::string display_type = "sdl";
    float fov = 11.0f;
    std::string db_stars = "../tetra3_db_stars.bin";
    std::string db_patterns = "../tetra3_db_patterns.bin";
    float detection_sigma = 5.0f;
    bool loop = false;
    bool stdin_json = false;
};

static AppConfig parse_args(int argc, char *argv[]) {
    AppConfig cfg;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--image" && i + 1 < argc) {
            cfg.image_path = argv[++i];
        } else if (arg == "--width" && i + 1 < argc) {
            cfg.width = std::stoi(argv[++i]);
        } else if (arg == "--height" && i + 1 < argc) {
            cfg.height = std::stoi(argv[++i]);
        } else if (arg == "--display" && i + 1 < argc) {
            cfg.display_type = argv[++i];
        } else if (arg == "--fov" && i + 1 < argc) {
            cfg.fov = std::stof(argv[++i]);
        } else if (arg == "--db-stars" && i + 1 < argc) {
            cfg.db_stars = argv[++i];
        } else if (arg == "--db-patterns" && i + 1 < argc) {
            cfg.db_patterns = argv[++i];
        } else if (arg == "--detection-sigma" && i + 1 < argc) {
            cfg.detection_sigma = std::stof(argv[++i]);
        } else if (arg == "--loop") {
            cfg.loop = true;
        } else if (arg == "--stdin-json") {
            cfg.stdin_json = true;
            cfg.loop = true; // stdin mode is always a loop
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: pi_tracker [options]\n"
                      << "  --image <path>           Load 16-bit PNG/TIFF for testing\n"
                      << "  --width <n>              Camera resolution width (default: 2028)\n"
                      << "  --height <n>             Camera resolution height (default: 1520)\n"
                      << "  --display sdl|i2c|none   Display backend (default: sdl)\n"
                      << "  --fov <degrees>          FOV estimate (default: 11.0)\n"
                      << "  --db-stars <path>        Star catalog path\n"
                      << "  --db-patterns <path>     Pattern catalog path\n"
                      << "  --detection-sigma <n>    Star detection threshold (default: 5.0)\n"
                      << "  --loop                   Continuous capture mode\n"
                      << "  --stdin-json             Read JSON solve results from stdin\n"
                      << "\nStdin JSON format (one per line):\n"
                      << "  {\"solved\":true,\"ra_rad\":1.5,\"dec_rad\":0.3,\"roll_rad\":0,\"fov_rad\":0.2}\n"
                      << "  {\"solved\":false}\n"
                      << "  {\"mode\":\"pa_sampling\"}  -- switch to PA sampling mode\n"
                      << "  {\"mode\":\"pa_fix\"}       -- switch to PA fix mode\n"
                      << "  {\"mode\":\"tracking\"}     -- switch to tracking mode\n"
                      << "\nKeyboard (SDL):\n"
                      << "  SPACE  Cycle: Tracking -> PA Sampling -> PA Fix -> Tracking\n"
                      << "  ESC    Return to Tracking mode\n";
            exit(0);
        }
    }
    return cfg;
}

// Simple input: poll keyboard for mode changes
enum class InputEvent { None, ModeButton, BackButton };

static InputEvent poll_input([[maybe_unused]] const std::string &display_type) {
#ifdef HAS_SDL_DISPLAY
    if (display_type == "sdl") {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) return InputEvent::BackButton;
            if (ev.type == SDL_KEYDOWN) {
                if (ev.key.keysym.sym == SDLK_SPACE) return InputEvent::ModeButton;
                if (ev.key.keysym.sym == SDLK_ESCAPE) return InputEvent::BackButton;
            }
        }
    }
#endif
    return InputEvent::None;
}

static double steady_now() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Minimal JSON value extraction (avoids pulling in a JSON library).
// Handles: "key":1.5  "key":true  "key":false  "key":"string"
static bool json_has_key(const std::string &json, const std::string &key) {
    return json.find("\"" + key + "\"") != std::string::npos;
}

static double json_float(const std::string &json, const std::string &key, double def = 0) {
    auto pos = json.find("\"" + key + "\"");
    if (pos == std::string::npos) return def;
    pos = json.find(':', pos);
    if (pos == std::string::npos) return def;
    return std::strtod(json.c_str() + pos + 1, nullptr);
}

static bool json_bool(const std::string &json, const std::string &key, bool def = false) {
    auto pos = json.find("\"" + key + "\"");
    if (pos == std::string::npos) return def;
    pos = json.find(':', pos);
    if (pos == std::string::npos) return def;
    // skip whitespace
    pos++;
    while (pos < json.size() && json[pos] == ' ') pos++;
    return json.compare(pos, 4, "true") == 0;
}

static std::string json_string(const std::string &json, const std::string &key) {
    auto pos = json.find("\"" + key + "\"");
    if (pos == std::string::npos) return "";
    pos = json.find(':', pos);
    if (pos == std::string::npos) return "";
    auto q1 = json.find('"', pos + 1);
    if (q1 == std::string::npos) return "";
    auto q2 = json.find('"', q1 + 1);
    if (q2 == std::string::npos) return "";
    return json.substr(q1 + 1, q2 - q1 - 1);
}

static const char B64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static std::string base64_encode(const uint8_t *data, size_t len) {
    std::string out;
    out.reserve((len + 2) / 3 * 4);
    for (size_t i = 0; i < len; i += 3) {
        uint32_t n = static_cast<uint32_t>(data[i]) << 16;
        if (i + 1 < len) n |= static_cast<uint32_t>(data[i + 1]) << 8;
        if (i + 2 < len) n |= data[i + 2];
        out += B64[(n >> 18) & 63];
        out += B64[(n >> 12) & 63];
        out += (i + 1 < len) ? B64[(n >> 6) & 63] : '=';
        out += (i + 2 < len) ? B64[n & 63] : '=';
    }
    return out;
}

static void render_display(IDisplay &display, StarMapRenderer &star_map,
                            PolarAligner &aligner, AppMode mode,
                            const SolveResult &result) {
    display.clear();

    switch (mode) {
    case AppMode::Tracking: {
        if (result.solved) {
            // Show the central half of the camera FOV on the map
            float map_fov = result.fov * 0.5f;
            star_map.render(display, result.ra, result.dec, map_fov,
                            0.0f, IDisplay::MAP_Y0, IDisplay::MAP_Y1);
        }
        display.show_solve_result(result);
        break;
    }

    case AppMode::PASampling: {
        auto pole = aligner.estimate_pole();
        display.show_pa_sampling(aligner.num_samples(), pole, result.solved);
        break;
    }

    case AppMode::PAFix: {
        auto pole = aligner.estimate_pole();
        auto offset = aligner.pole_error(pole);
        display.show_pa_fix(offset);
        break;
    }
    }

    display.update();
}

static void handle_mode_switch(InputEvent event, AppMode &mode,
                                PolarAligner &aligner) {
    if (event == InputEvent::ModeButton) {
        switch (mode) {
            case AppMode::Tracking:
                mode = AppMode::PASampling;
                aligner.clear_samples();
                std::cerr << "Mode: PA Sampling" << std::endl;
                break;
            case AppMode::PASampling:
                if (aligner.num_samples() >= 3) {
                    mode = AppMode::PAFix;
                    std::cerr << "Mode: PA Fix" << std::endl;
                }
                break;
            case AppMode::PAFix:
                mode = AppMode::Tracking;
                std::cerr << "Mode: Tracking" << std::endl;
                break;
        }
    } else if (event == InputEvent::BackButton) {
        mode = AppMode::Tracking;
        std::cerr << "Mode: Tracking" << std::endl;
    }
}

// --stdin-json main loop: reads JSON lines from stdin, drives UI
static int run_stdin_json(AppConfig &cfg, IDisplay &display) {
    StarMapRenderer star_map;
    PolarAligner aligner(OBSERVER_LAT_DEG, OBSERVER_LON_DEG);
    AppMode mode = AppMode::Tracking;

    display.show_status("Waiting...");

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;

        // Poll SDL events to keep window alive
        auto event = poll_input(cfg.display_type);
        if (event != InputEvent::None) {
            handle_mode_switch(event, mode, aligner);
        }

        // Check for mode command
        auto mode_str = json_string(line, "mode");
        if (!mode_str.empty()) {
            if (mode_str == "pa_sampling") {
                mode = AppMode::PASampling;
                aligner.clear_samples();
                std::cerr << "Mode: PA Sampling" << std::endl;
            } else if (mode_str == "pa_fix") {
                if (aligner.num_samples() >= 3) {
                    mode = AppMode::PAFix;
                    std::cerr << "Mode: PA Fix" << std::endl;
                }
            } else if (mode_str == "tracking") {
                mode = AppMode::Tracking;
                std::cerr << "Mode: Tracking" << std::endl;
            }

            // Re-render with current state
            SolveResult dummy{};
            render_display(display, star_map, aligner, mode, dummy);
            // Ack the command with framebuffer
            std::string fb_str;
            if (auto *fb = display.get_framebuffer())
                fb_str = base64_encode(fb, IDisplay::FB_SIZE);
            std::cout << "{\"ack\":\"" << mode_str << "\",\"fb\":\"" << fb_str << "\"}" << std::endl;
            continue;
        }

        // Parse solve result
        SolveResult result{};
        result.solved = json_bool(line, "solved");
        if (result.solved) {
            result.ra = static_cast<float>(json_float(line, "ra_rad"));
            result.dec = static_cast<float>(json_float(line, "dec_rad"));
            result.roll = static_cast<float>(json_float(line, "roll_rad"));
            result.fov = static_cast<float>(json_float(line, "fov_rad", 0.2));
            result.rmse = static_cast<float>(json_float(line, "rmse"));
            result.num_matches = static_cast<int>(json_float(line, "num_matches"));
            result.solve_time_ms = static_cast<float>(json_float(line, "solve_time_ms"));
        }

        // Feed to polar alignment if in sampling/fix mode
        if (result.solved && (mode == AppMode::PASampling || mode == AppMode::PAFix)) {
            aligner.add_sample(result.ra, result.dec, steady_now());
        }

        render_display(display, star_map, aligner, mode, result);

        // Echo back status + framebuffer
        std::string fb_str;
        if (auto *fb = display.get_framebuffer())
            fb_str = base64_encode(fb, IDisplay::FB_SIZE);

        if (mode == AppMode::PAFix) {
            auto pole = aligner.estimate_pole();
            auto offset = aligner.pole_error(pole);
            std::cout << "{\"mode\":\"pa_fix\",\"alt_arcmin\":"
                      << offset.alt_arcmin << ",\"az_arcmin\":"
                      << offset.az_arcmin << ",\"total_arcmin\":"
                      << offset.total_arcmin
                      << ",\"fb\":\"" << fb_str << "\"}" << std::endl;
        } else if (mode == AppMode::PASampling) {
            std::cout << "{\"mode\":\"pa_sampling\",\"num_samples\":"
                      << aligner.num_samples()
                      << ",\"fb\":\"" << fb_str << "\"}" << std::endl;
        } else {
            std::cout << "{\"mode\":\"tracking\",\"solved\":"
                      << (result.solved ? "true" : "false")
                      << ",\"fb\":\"" << fb_str << "\"}" << std::endl;
        }
    }

    return 0;
}

// Normal main loop: capture, detect, solve
static int run_normal(AppConfig &cfg, IDisplay &display) {
    // Load catalogs
    BinaryDatabaseLoader loader;
    if (!loader.load_database(cfg.db_stars, cfg.db_patterns)) {
        std::cerr << "Failed to load database!" << std::endl;
        display.show_status("DB load fail");
        return 1;
    }
    loader.print_info();

    SimpleStarSolver solver;
    solver.load_star_catalog(loader.get_stars());
    solver.load_pattern_catalog(loader.get_patterns());
    solver.print_memory_usage();

    // Create frame source
    std::unique_ptr<IFrameSource> source;
    if (!cfg.image_path.empty()) {
        source = std::make_unique<ImageFileSource>(cfg.image_path);
    } else {
#ifdef HAS_LIBCAMERA
        source = std::make_unique<LibcameraCapture>(cfg.width, cfg.height);
#else
        std::cerr << "No image specified and libcamera not available. Use --image <path>." << std::endl;
        return 1;
#endif
    }

    if (!source->initialize()) {
        std::cerr << "Failed to initialize frame source" << std::endl;
        display.show_status("Source fail");
        return 1;
    }

    StarDetectorConfig det_cfg;
    det_cfg.detection_sigma = cfg.detection_sigma;
    StarDetector detector(det_cfg);
    StarMapRenderer star_map;
    PolarAligner aligner(OBSERVER_LAT_DEG, OBSERVER_LON_DEG);
    AppMode mode = AppMode::Tracking;

    do {
        auto event = poll_input(cfg.display_type);
        if (event == InputEvent::ModeButton) {
            handle_mode_switch(event, mode, aligner);
        } else if (event == InputEvent::BackButton) {
            if (mode != AppMode::Tracking) {
                handle_mode_switch(event, mode, aligner);
            } else {
                break;
            }
        }

        display.show_status("Capturing...");

        auto frame_opt = source->capture();
        if (!frame_opt) {
            std::cerr << "Failed to capture frame" << std::endl;
            display.show_status("Capture fail");
            continue;
        }

        auto &frame = *frame_opt;
        std::cerr << "Frame: " << frame.width << "x" << frame.height
                  << " " << frame.bit_depth << "-bit" << std::endl;

        display.show_status("Detecting...");
        auto detected_stars = detector.detect(frame);

        if (detected_stars.empty()) {
            std::cerr << "No stars detected" << std::endl;
            display.show_status("No stars");
            continue;
        }

        std::vector<Centroid> centroids;
        centroids.reserve(detected_stars.size());
        for (const auto &s : detected_stars) {
            centroids.push_back(s.centroid);
        }

        display.show_status("Solving...");
        auto result = solver.solve_from_centroids(
            centroids, frame.height, frame.width, cfg.fov, 16,
            0.01f, 0.001f, std::nullopt, cfg.fov * 0.2f);

        if (result.solved && (mode == AppMode::PASampling || mode == AppMode::PAFix)) {
            aligner.add_sample(result.ra, result.dec, steady_now());
        }

        render_display(display, star_map, aligner, mode, result);

        if (result.solved) {
            std::cerr << "RA: " << result.ra * 180.0f / M_PI
                      << " Dec: " << result.dec * 180.0f / M_PI
                      << " FOV: " << result.fov * 180.0f / M_PI
                      << " RMSE: " << result.rmse << "\"" << std::endl;
        } else {
            std::cerr << "No solution (" << result.solve_time_ms << " ms)" << std::endl;
        }

        if (!cfg.loop) {
#ifdef HAS_SDL_DISPLAY
            if (cfg.display_type == "sdl") {
                std::cerr << "Press Ctrl+C to exit (SDL window open)..." << std::endl;
                bool running = true;
                while (running) {
                    SDL_Event ev;
                    while (SDL_PollEvent(&ev)) {
                        if (ev.type == SDL_QUIT || ev.type == SDL_KEYDOWN) {
                            running = false;
                        }
                    }
                    SDL_Delay(100);
                }
            }
#endif
        }

    } while (cfg.loop);

    return 0;
}

int main(int argc, char *argv[]) {
    AppConfig cfg = parse_args(argc, argv);

    // Create display
    std::unique_ptr<IDisplay> display;
    if (cfg.stdin_json && cfg.display_type != "sdl") {
        // For stdin-json mode, default to buffer display (has framebuffer, no GUI)
        display = std::make_unique<BufferDisplay>();
    } else if (cfg.display_type == "none") {
        display = std::make_unique<NullDisplay>();
    } else if (cfg.display_type == "i2c") {
        display = std::make_unique<I2cSh1106Display>();
#ifdef HAS_SDL_DISPLAY
    } else if (cfg.display_type == "sdl") {
        display = std::make_unique<SdlDisplay>();
#endif
    } else if (cfg.stdin_json) {
        display = std::make_unique<BufferDisplay>();
    } else {
        std::cerr << "Unknown display type: " << cfg.display_type << std::endl;
        return 1;
    }

    display->show_status("Loading...");

    if (cfg.stdin_json) {
        return run_stdin_json(cfg, *display);
    } else {
        return run_normal(cfg, *display);
    }
}
