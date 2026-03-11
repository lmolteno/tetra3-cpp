#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <string>
#include <thread>
#include <chrono>

#include <tetra3/types.h>
#include "display.h"
#include "display_null.h"
#include "star_map.h"
#include "polar_align.h"

#ifdef HAS_SDL_DISPLAY
#include "display_sdl.h"
#include <SDL2/SDL.h>
#endif

static bool wait_or_key(int ms) {
#ifdef HAS_SDL_DISPLAY
    auto end = std::chrono::steady_clock::now() + std::chrono::milliseconds(ms);
    while (std::chrono::steady_clock::now() < end) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) exit(0);
            if (ev.type == SDL_KEYDOWN) return true;
        }
        SDL_Delay(30);
    }
#else
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
#endif
    return false;
}

static void wait_for_exit() {
#ifdef HAS_SDL_DISPLAY
    std::cout << "Press any key or close window to continue." << std::endl;
    bool running = true;
    while (running) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT || ev.type == SDL_KEYDOWN)
                running = false;
        }
        SDL_Delay(50);
    }
#endif
}

// ---------------------------------------------------------------------------
// Large-format SDL canvas for debugging star map rendering
// ---------------------------------------------------------------------------
#ifdef HAS_SDL_DISPLAY
class BigSdlCanvas : public IStarMapCanvas {
public:
    int width, height;

    BigSdlCanvas(int w, int h, int scale = 2)
        : width(w), height(h), scale_(scale), pixels_(w * h, 0) {
        SDL_Init(SDL_INIT_VIDEO);
        window_ = SDL_CreateWindow("Star Map (big)",
                                    SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                                    w * scale, h * scale, SDL_WINDOW_SHOWN);
        renderer_ = SDL_CreateRenderer(window_, -1, SDL_RENDERER_ACCELERATED);
        texture_ = SDL_CreateTexture(renderer_, SDL_PIXELFORMAT_RGB332,
                                      SDL_TEXTUREACCESS_STREAMING, w, h);
    }

    ~BigSdlCanvas() {
        if (texture_)  SDL_DestroyTexture(texture_);
        if (renderer_) SDL_DestroyRenderer(renderer_);
        if (window_)   SDL_DestroyWindow(window_);
    }

    void draw_pixel(int x, int y) override {
        if (x >= 0 && x < width && y >= 0 && y < height)
            pixels_[y * width + x] = 0xFF;
    }

    void clear() { std::fill(pixels_.begin(), pixels_.end(), 0); }

    void update() {
        SDL_UpdateTexture(texture_, nullptr, pixels_.data(), width);
        SDL_RenderClear(renderer_);
        SDL_RenderCopy(renderer_, texture_, nullptr, nullptr);
        SDL_RenderPresent(renderer_);
    }

    // Draw text at pixel position (simple, using SDL -- just label it)
    void draw_label(int x, int y, const char *text) {
        // No font available in this canvas, just skip
        // Labels go to stdout instead
    }

private:
    int scale_;
    std::vector<uint8_t> pixels_;
    SDL_Window *window_ = nullptr;
    SDL_Renderer *renderer_ = nullptr;
    SDL_Texture *texture_ = nullptr;
};
#endif

// ---------------------------------------------------------------------------
// Star map demos
// ---------------------------------------------------------------------------

static void run_starmap_demo(IDisplay &display) {
    StarMapRenderer map;

    display.show_status("Loading...");
    wait_or_key(1000);

    struct Scene { const char *name; float ra_deg, dec_deg, fov_deg; };
    Scene scenes[] = {
        {"Orion",          84.05f,  -1.2f, 40.0f},
        {"Southern Cross", 186.5f, -63.0f, 30.0f},
        {"Ursa Major",     180.0f,  55.0f, 35.0f},
        {"Scorpius",       247.5f, -30.0f, 50.0f},
        {"Cassiopeia",      15.0f,  60.0f, 40.0f},
    };

    for (auto &s : scenes) {
        float ra = s.ra_deg * static_cast<float>(M_PI) / 180.0f;
        float dec = s.dec_deg * static_cast<float>(M_PI) / 180.0f;
        float fov = s.fov_deg * static_cast<float>(M_PI) / 180.0f;

        display.clear();
        map.render(display, ra, dec, fov, 0.0f, IDisplay::MAP_Y0, IDisplay::MAP_Y1);
        display.show_coordinates(ra, dec);
        display.update();
        std::cout << "  " << s.name << std::endl;
        wait_or_key(3000);
    }

    display.clear();
    SolveResult unsolved{};
    unsolved.solved = false;
    display.show_solve_result(unsolved);
    display.draw_text(2, 24, "Waiting for solve...");
    display.update();

    wait_for_exit();
}

#ifdef HAS_SDL_DISPLAY
static void run_big_starmap_demo(int size) {
    StarMapRenderer map;
    BigSdlCanvas canvas(size, size, 2);

    struct Scene { const char *name; float ra_deg, dec_deg, fov_deg; };
    Scene scenes[] = {
        {"Orion",          84.05f,  -1.2f, 75.0f},
        {"Southern Cross", 186.5f, -63.0f, 75.0f},
        {"Ursa Major",     180.0f,  55.0f, 75.0f},
        {"Scorpius",       247.5f, -30.0f, 75.0f},
        {"Cassiopeia",      15.0f,  60.0f, 75.0f},
    };

    for (auto &s : scenes) {
        float ra = s.ra_deg * static_cast<float>(M_PI) / 180.0f;
        float dec = s.dec_deg * static_cast<float>(M_PI) / 180.0f;
        float fov = s.fov_deg * static_cast<float>(M_PI) / 180.0f;

        canvas.clear();
        map.render(canvas, size, size, ra, dec, fov);
        canvas.update();
        std::cout << "  " << s.name
                  << " (RA=" << s.ra_deg << " Dec=" << s.dec_deg
                  << " FOV=" << s.fov_deg << ")" << std::endl;

        if (wait_or_key(5000)) {
            // key pressed, advance immediately
        }
    }

    std::cout << "Done. Press any key to exit." << std::endl;
    wait_for_exit();
}
#endif

// ---------------------------------------------------------------------------
// Polar alignment demo
// ---------------------------------------------------------------------------

static void run_polar_align_demo(IDisplay &display) {
    PolarAligner aligner(OBSERVER_LAT_DEG, OBSERVER_LON_DEG);
    bool southern = OBSERVER_LAT_DEG < 0;

    float pole_dec_deg = southern ? -88.5f : 88.5f;
    float pole_ra_deg  = southern ? 170.0f : 10.0f;
    float pole_ra  = pole_ra_deg * static_cast<float>(M_PI) / 180.0f;
    float pole_dec = pole_dec_deg * static_cast<float>(M_PI) / 180.0f;

    float circle_radius = 20.0f * static_cast<float>(M_PI) / 180.0f;

    constexpr int NUM_SAMPLES = 6;
    float rotation_angles[NUM_SAMPLES];
    for (int i = 0; i < NUM_SAMPLES; i++) {
        rotation_angles[i] = static_cast<float>(i) * 15.0f * static_cast<float>(M_PI) / 180.0f;
    }

    float px = std::cos(pole_dec) * std::cos(pole_ra);
    float py = std::cos(pole_dec) * std::sin(pole_ra);
    float pz = std::sin(pole_dec);

    float e1x = -std::sin(pole_ra);
    float e1y =  std::cos(pole_ra);
    float e1z =  0.0f;
    float e2x = py * e1z - pz * e1y;
    float e2y = pz * e1x - px * e1z;
    float e2z = px * e1y - py * e1x;

    double t0 = std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();

    struct SynthSample { float ra, dec; double timestamp; };
    SynthSample samples[NUM_SAMPLES];

    float sin_r = std::sin(circle_radius);
    float cos_r = std::cos(circle_radius);
    for (int i = 0; i < NUM_SAMPLES; i++) {
        float angle = rotation_angles[i];
        float sx = cos_r * px + sin_r * (std::cos(angle) * e1x + std::sin(angle) * e2x);
        float sy = cos_r * py + sin_r * (std::cos(angle) * e1y + std::sin(angle) * e2y);
        float sz = cos_r * pz + sin_r * (std::cos(angle) * e1z + std::sin(angle) * e2z);

        float dec = std::asin(std::clamp(sz, -1.0f, 1.0f));
        float ra  = std::atan2(sy, sx);
        if (ra < 0) ra += 2.0f * static_cast<float>(M_PI);

        double dt = i * 30.0;
        ra += static_cast<float>(7.2921159e-5 * dt);

        samples[i] = {ra, dec, t0 + dt};
    }

    std::cout << "Polar Alignment Demo" << std::endl;
    std::cout << "  Simulated mount pole: RA=" << pole_ra_deg
              << " Dec=" << pole_dec_deg << std::endl;
    std::cout << "  Observer: lat=" << OBSERVER_LAT_DEG
              << " lon=" << OBSERVER_LON_DEG << std::endl;

    for (int i = 0; i < NUM_SAMPLES; i++) {
        aligner.add_sample(samples[i].ra, samples[i].dec, samples[i].timestamp);

        display.clear();
        auto pole = aligner.estimate_pole();
        display.show_pa_sampling(aligner.num_samples(), pole, true);
        display.update();

        std::cout << "  Sample " << (i + 1) << "/" << NUM_SAMPLES
                  << ": RA=" << samples[i].ra * 180.0f / M_PI
                  << " Dec=" << samples[i].dec * 180.0f / M_PI;
        if (pole.valid) {
            std::cout << "  -> Pole est: RA=" << pole.ra * 180.0f / M_PI
                      << " Dec=" << pole.dec * 180.0f / M_PI
                      << " arc=" << pole.arc_deg << "deg";
        }
        std::cout << std::endl;

        wait_or_key(2000);
    }

    std::cout << "  Entering fix mode..." << std::endl;

    auto pole = aligner.estimate_pole();
    auto offset = aligner.pole_error(pole);

    display.clear();
    display.show_pa_fix(offset);
    display.update();

    std::cout << "  Alt offset: " << offset.alt_arcmin << "'" << std::endl;
    std::cout << "  Az  offset: " << offset.az_arcmin << "'" << std::endl;
    std::cout << "  Total:      " << offset.total_arcmin << "'" << std::endl;

    if (pole.valid) {
        std::cout << "  Estimated pole: RA=" << pole.ra * 180.0f / M_PI
                  << " Dec=" << pole.dec * 180.0f / M_PI << std::endl;
    }

    wait_for_exit();
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char *argv[]) {
#ifndef HAS_SDL_DISPLAY
    std::cerr << "SDL not available" << std::endl;
    return 1;
#else
    std::string mode = "starmap";
    int big_size = 0;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--polar-align" || arg == "-p") {
            mode = "polar";
        } else if (arg == "--starmap" || arg == "-s") {
            mode = "starmap";
        } else if (arg == "--big" || arg == "-b") {
            mode = "big";
            big_size = 256;
            if (i + 1 < argc) {
                try { big_size = std::stoi(argv[i + 1]); i++; }
                catch (...) {} // not a number, use default
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: display_test [options]\n"
                      << "  --starmap, -s         Star map demo on 128x64 OLED (default)\n"
                      << "  --big [N], -b [N]     Large NxN star map (default 256)\n"
                      << "  --polar-align, -p     Polar alignment demo\n";
            return 0;
        }
    }

    if (mode == "big") {
        std::cout << "Big star map: " << big_size << "x" << big_size << std::endl;
        run_big_starmap_demo(big_size);
    } else if (mode == "polar") {
        auto display = std::make_unique<SdlDisplay>(4);
        run_polar_align_demo(*display);
    } else {
        auto display = std::make_unique<SdlDisplay>(4);
        run_starmap_demo(*display);
    }

    return 0;
#endif
}
