// pi_tracker daemon: capture -> archive -> solve (out-of-process) -> render -> publish.
//
//   stdin  : NDJSON commands  ({"mode":...} | {"lens_pos":...} | {"exposure_us":...} | {"gain":...})
//   stdout : NDJSON status    (tick state + base64 framebuffer)
//   tmpfs  : rolling 10-frame archive (defaults to /dev/shm)
//   state  : lens_position, exposure_us, gain (persistent, on the SD card)
//
// Solving runs as a `solve_cli --batch` subprocess; one solve in flight at a
// time, with frames dropped while the solver is busy so the publish loop
// stays at frame rate.

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <thread>

#include "capture/frame_source.h"
#include "capture/image_loader.h"
#include "render/buffer_canvas.h"
#include "render/i2c_ssd1306.h"
#include "render/tracker_view.h"
#include "app/polar_align.h"
#include "io/json_sink.h"
#include "capture/frame_writer.h"
#include "solver/solver_process.h"
#include "io/control_input.h"
#include "app/cam_state.h"

#ifdef HAS_LIBCAMERA
#include "capture/camera_capture.h"
#endif

struct AppConfig {
    std::string image_path;
    int cam_width  = 2028;
    int cam_height = 1520;
    float lat = OBSERVER_LAT_DEG;
    float lon = OBSERVER_LON_DEG;
    std::string frames_dir = "/dev/shm/pi_tracker/frames";
    std::string state_dir  = "./state";
    float max_fps = 5.0f;

    // Camera overrides (optional — saved state takes precedence otherwise).
    bool     have_lens_pos    = false;  float    lens_pos    = 0.0f;
    bool     have_exposure_us = false;  uint32_t exposure_us = 0;
    bool     have_gain        = false;  float    gain        = 0.0f;

    // Solver subprocess
    std::string solver_path = "solve_cli";
    std::string db_stars;
    std::string db_patterns;
    float fov_deg = 11.0f;
    float detection_sigma = 5.0f;
    int crop_size = 720;
    int bg_tile = 128;
    // Lens distortion correction (Camera Module 3 NoIR Standard, calibrated).
    // Pass NaN to disable. See SolverProcess::Config.distortion_k.
    float distortion_k = 0.024f;
    bool no_solver = false;
    bool oled = false;
    std::string oled_bus = "/dev/i2c-1";
};

static AppConfig parse_args(int argc, char *argv[]) {
    AppConfig cfg;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if      (arg == "--image"      && i + 1 < argc) cfg.image_path  = argv[++i];
        else if (arg == "--width"      && i + 1 < argc) cfg.cam_width   = std::stoi(argv[++i]);
        else if (arg == "--height"     && i + 1 < argc) cfg.cam_height  = std::stoi(argv[++i]);
        else if (arg == "--lat"        && i + 1 < argc) cfg.lat         = std::stof(argv[++i]);
        else if (arg == "--lon"        && i + 1 < argc) cfg.lon         = std::stof(argv[++i]);
        else if (arg == "--frames-dir" && i + 1 < argc) cfg.frames_dir  = argv[++i];
        else if (arg == "--state-dir"  && i + 1 < argc) cfg.state_dir   = argv[++i];
        else if (arg == "--max-fps"    && i + 1 < argc) cfg.max_fps     = std::stof(argv[++i]);
        else if (arg == "--lens-pos"    && i + 1 < argc) { cfg.lens_pos    = std::stof(argv[++i]); cfg.have_lens_pos    = true; }
        else if (arg == "--exposure-us" && i + 1 < argc) { cfg.exposure_us = static_cast<uint32_t>(std::stoul(argv[++i])); cfg.have_exposure_us = true; }
        else if (arg == "--gain"        && i + 1 < argc) { cfg.gain        = std::stof(argv[++i]); cfg.have_gain        = true; }
        else if (arg == "--solver"     && i + 1 < argc) cfg.solver_path = argv[++i];
        else if (arg == "--db-stars"   && i + 1 < argc) cfg.db_stars    = argv[++i];
        else if (arg == "--db-patterns"&& i + 1 < argc) cfg.db_patterns = argv[++i];
        else if (arg == "--fov"        && i + 1 < argc) cfg.fov_deg     = std::stof(argv[++i]);
        else if (arg == "--detection-sigma" && i + 1 < argc) cfg.detection_sigma = std::stof(argv[++i]);
        else if (arg == "--crop"       && i + 1 < argc) cfg.crop_size   = std::stoi(argv[++i]);
        else if (arg == "--bg-tile"    && i + 1 < argc) cfg.bg_tile     = std::stoi(argv[++i]);
        else if (arg == "--distortion-k"     && i + 1 < argc) cfg.distortion_k = std::stof(argv[++i]);
        else if (arg == "--no-distortion-k")            cfg.distortion_k = std::nanf("");
        else if (arg == "--no-solver")                  cfg.no_solver   = true;
        else if (arg == "--oled")                       cfg.oled        = true;
        else if (arg == "--oled-bus"   && i + 1 < argc) cfg.oled_bus   = argv[++i];
        else if (arg == "--help" || arg == "-h") {
            std::cerr <<
                "Usage: pi_tracker [options]\n"
                "  --image <path>            Use a test image instead of libcamera\n"
                "  --width/--height <n>      Camera capture dims (default 2028x1520)\n"
                "  --lat/--lon <deg>         Observer location (default compile-time)\n"
                "  --frames-dir <path>       Rolling frame archive (default /dev/shm/pi_tracker/frames)\n"
                "  --state-dir <path>        Persistent state dir  (default ./state)\n"
                "  --max-fps <n>             Output rate cap       (default 5)\n"
                "  --lens-pos <diopters>     Manual focus (overrides saved state)\n"
                "  --exposure-us <us>        Exposure time (overrides saved state)\n"
                "  --gain <x>                Analogue gain  (overrides saved state)\n"
                "  --solver <path>           Path to solve_cli     (default 'solve_cli' on PATH)\n"
                "  --db-stars/--db-patterns  Star catalog files (passed through to solve_cli)\n"
                "  --fov <deg>               FOV hint for the solver (default 11)\n"
                "  --detection-sigma <n>     Star detection threshold (default 5)\n"
                "  --crop <n>                Solver crop size (default 720)\n"
                "  --bg-tile <n>             Solver bg-tile size (default 128)\n"
                "  --distortion-k <k>        Frozen radial distortion coefficient passed to\n"
                "                            solve_cli --freeze-distortion-k (default 0.024,\n"
                "                            calibrated for Camera Module 3 NoIR Standard)\n"
                "  --no-distortion-k         Disable distortion correction (fall back to k=0)\n"
                "  --no-solver               Don't spawn solve_cli; emit no-solve frames only\n"
                "  --oled                    Push framebuffer to SSD1306 OLED via I2C (0x3C)\n"
                "  --oled-bus <dev>          I2C bus device (default /dev/i2c-1)\n";
            std::exit(0);
        }
    }
    return cfg;
}

static std::unique_ptr<IFrameSource> make_source(const AppConfig &cfg,
                                                  const cam_state::Settings &s) {
    if (!cfg.image_path.empty()) {
        return std::make_unique<ImageFileSource>(cfg.image_path);
    }
#ifdef HAS_LIBCAMERA
    auto cam = std::make_unique<LibcameraCapture>(cfg.cam_width, cfg.cam_height);
    cam->set_lens_position(s.lens_pos);
    cam->set_exposure_us(s.exposure_us);
    cam->set_gain(s.gain);
    return cam;
#else
    (void)s;
    std::cerr << "No --image given and libcamera not available at build time\n";
    return nullptr;
#endif
}

// Allow Ctrl-C to break the loop cleanly.
static std::atomic<bool> g_running{true};
static void on_sigint(int) { g_running = false; }

int main(int argc, char *argv[]) {
    AppConfig cfg = parse_args(argc, argv);

    // Quiet libcamera's INFO / WARN chatter to stderr by default. Otherwise
    // running over `ssh -tt host "./pi_tracker ..." | view.py` floods the
    // viewer's stdin (pty merges stderr into stdout under -tt). User can
    // still override via env if they want the noise.
    setenv("LIBCAMERA_LOG_LEVELS", "*:ERROR", /*overwrite=*/0);

    std::signal(SIGINT, on_sigint);
    std::signal(SIGTERM, on_sigint);
    std::signal(SIGPIPE, SIG_IGN);  // writing to a closed stdout shouldn't kill us

    BufferCanvas canvas;
    TrackerView view;

    std::unique_ptr<I2cSsd1306> oled;
    if (cfg.oled) {
        try {
            oled = std::make_unique<I2cSsd1306>(cfg.oled_bus.c_str());
        } catch (const std::exception &e) {
            std::cerr << "OLED init failed: " << e.what() << '\n';
        }
    }
    PolarAligner aligner(cfg.lat, cfg.lon);
    AppMode mode = AppMode::Tracking;
    ControlInput input;

    cam_state::Settings cs = cam_state::load(cfg.state_dir);
    // CLI flags override saved state (and persist on use).
    if (cfg.have_lens_pos)    { cs.lens_pos    = cfg.lens_pos;    cam_state::save_lens_pos   (cfg.state_dir, cs.lens_pos);    }
    if (cfg.have_exposure_us) { cs.exposure_us = cfg.exposure_us; cam_state::save_exposure_us(cfg.state_dir, cs.exposure_us); }
    if (cfg.have_gain)        { cs.gain        = cfg.gain;        cam_state::save_gain       (cfg.state_dir, cs.gain);        }

    auto fill_state = [&](json_sink::State &s) {
        s.mode = mode;
        s.lens_pos    = cs.lens_pos;
        s.exposure_us = cs.exposure_us;
        s.gain        = cs.gain;
    };

    auto announce = [&](const std::string &msg) {
        view.render_status(canvas, msg);
        json_sink::State s;
        fill_state(s);
        s.fb = canvas.framebuffer();
        s.fb_size = canvas.framebuffer_size();
        json_sink::publish(std::cout, s);
    };

    announce("Starting");

    auto source = make_source(cfg, cs);
    if (!source) return 1;
    if (!source->initialize()) {
        announce("Source fail");
        return 1;
    }

    FrameWriter writer(cfg.frames_dir);
    if (!writer.ensure_dir()) {
        announce("Frames dir fail");
        return 1;
    }

    std::unique_ptr<SolverProcess> solver;
    if (!cfg.no_solver) {
        SolverProcess::Config sc;
        sc.binary_path     = cfg.solver_path;
        sc.db_stars        = cfg.db_stars;
        sc.db_patterns     = cfg.db_patterns;
        sc.fov_deg         = cfg.fov_deg;
        sc.detection_sigma = cfg.detection_sigma;
        sc.crop_size       = cfg.crop_size;
        sc.bg_tile         = cfg.bg_tile;
        sc.distortion_k    = cfg.distortion_k;
        solver = std::make_unique<SolverProcess>(std::move(sc));
        if (!solver->start()) {
            announce("Solver fail");
            return 1;
        }
    }

    const auto frame_period = std::chrono::nanoseconds(
        static_cast<int64_t>(1e9 / cfg.max_fps));
    auto next_tick = std::chrono::steady_clock::now();

    SolveResult last_result{};
    last_result.solved = false;
    std::string last_frame_path;
    std::string last_archived_path;
    json_sink::Timing last_solve_timing;          // load/detect/match/solve_wall from latest result
    std::chrono::steady_clock::time_point last_result_at{};
    std::chrono::steady_clock::time_point last_frame_at{};

    auto now_steady = []() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
    };

    // Apply a camera tweak if we're actually driving a libcamera source. Falls
    // through to a no-op for the file-source case (or when libcamera wasn't
    // available at build time).
    auto apply_to_camera = [&](auto fn) {
#ifdef HAS_LIBCAMERA
        if (auto *cam = dynamic_cast<LibcameraCapture *>(source.get())) fn(cam);
#else
        (void)fn;
#endif
    };

    while (g_running) {
        next_tick += frame_period;

        // 1. Drain any pending stdin commands.
        for (const auto &cmd : input.poll()) {
            if (cmd.mode) {
                AppMode want = *cmd.mode;
                if (want == AppMode::PAFix && aligner.num_samples() < 3) {
                    // Need samples first; ignore promotion.
                } else {
                    if (want == AppMode::PASampling && mode != AppMode::PASampling) {
                        aligner.clear_samples();
                    }
                    mode = want;
                }
            }
            if (cmd.lens_pos) {
                cs.lens_pos = *cmd.lens_pos;
                cam_state::save_lens_pos(cfg.state_dir, cs.lens_pos);
                apply_to_camera([&](auto *cam) { cam->set_lens_position(cs.lens_pos); });
            }
            if (cmd.exposure_us) {
                cs.exposure_us = *cmd.exposure_us;
                cam_state::save_exposure_us(cfg.state_dir, cs.exposure_us);
                apply_to_camera([&](auto *cam) { cam->set_exposure_us(cs.exposure_us); });
            }
            if (cmd.gain) {
                cs.gain = *cmd.gain;
                cam_state::save_gain(cfg.state_dir, cs.gain);
                apply_to_camera([&](auto *cam) { cam->set_gain(cs.gain); });
            }
        }

        using clock = std::chrono::steady_clock;
        auto ms_since = [](clock::time_point t) {
            return std::chrono::duration<float, std::milli>(clock::now() - t).count();
        };

        // 2. Try to grab a fresh frame. capture() is non-blocking now — the
        //    libcamera path runs an internal capture thread and returns the
        //    most recent unread frame, or nullopt if nothing new since last
        //    call. Image-source returns the loaded frame each call.
        auto t_capture = clock::now();
        auto frame = source->capture();
        float capture_ms = ms_since(t_capture);

        std::string frame_path;
        float archive_ms = -1.0f;
        if (frame) {
            auto t_archive = clock::now();
            frame_path = writer.write(*frame);
            archive_ms = ms_since(t_archive);
            if (!frame_path.empty()) {
                last_archived_path = frame_path;
                last_frame_at = clock::now();
            }
        }

        // 3. Submit to the solver if we have a new frame and it's free.
        if (frame && solver && !frame_path.empty()) {
            solver->submit(frame_path);
        }

        // 4. Pick up any fresh result.
        if (solver) {
            if (auto rs = solver->poll_fresh()) {
                last_result = rs->result;
                last_frame_path = rs->frame_path;
                last_solve_timing = json_sink::Timing{};
                last_solve_timing.load       = rs->load_ms;
                last_solve_timing.detect     = rs->detect_ms;
                last_solve_timing.bg         = rs->bg_ms;
                last_solve_timing.mask       = rs->mask_ms;
                last_solve_timing.cluster    = rs->cluster_ms;
                last_solve_timing.match      = last_result.solve_time_ms;
                last_solve_timing.solve_wall = rs->wall_ms;
                last_result_at = rs->received_at;
                if (last_result.solved &&
                    (mode == AppMode::PASampling || mode == AppMode::PAFix)) {
                    aligner.add_sample(last_result.ra, last_result.dec, now_steady());
                }
            }
        }

        // 5. Render + publish. We render every tick (even if no fresh frame
        //    arrived) so commands and solver results show up at publish rate.
        view.render(canvas, mode, last_result, aligner);
        if (oled) {
            try { oled->flush(canvas.framebuffer()); }
            catch (const std::exception &e) {
                std::cerr << "OLED flush: " << e.what() << '\n';
                oled.reset();  // stop trying after a hard I2C error
            }
        }

        json_sink::State s;
        fill_state(s);
        s.result = last_result;
        s.pa_samples = aligner.num_samples();
        s.pole = aligner.estimate_pole();
        s.offset = aligner.pole_error(s.pole);
        // Always point at the most recently archived frame so the web UI can
        // render even on stale ticks.
        s.frame_path = frame_path.empty() ? last_archived_path : frame_path;
        if (s.frame_path.empty()) s.frame_path = last_frame_path;
        s.timing = last_solve_timing;
        s.timing.capture = capture_ms;
        s.timing.archive = archive_ms;
        apply_to_camera([&](auto *cam) {
            s.timing.wait = cam->last_wait_ms();
            s.timing.bin  = cam->last_bin_ms();
        });
        if (last_frame_at != clock::time_point{}) {
            s.timing.frame_age = ms_since(last_frame_at);
        }
        if (last_result_at != clock::time_point{}) {
            s.timing.result_age = ms_since(last_result_at);
        }
        s.fb = canvas.framebuffer();
        s.fb_size = canvas.framebuffer_size();
        json_sink::publish(std::cout, s);

        std::this_thread::sleep_until(next_tick);
    }

    return 0;
}
