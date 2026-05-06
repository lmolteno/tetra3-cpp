// Renders every TrackerView UI state to stdout as ASCII art.
// With --oled it also pushes each frame to a real SSD1306 (2-second hold).
// No camera, solver, or live data needed.
//
// Usage: oled_test [--oled] [--oled-bus /dev/i2c-N]

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <thread>
#include <chrono>

#include "render/buffer_canvas.h"
#include "render/i2c_ssd1306.h"
#include "render/tracker_view.h"
#include "app/polar_align.h"
#include <tetra3/types.h>

static void print_fb(const uint8_t *fb) {
    // Each terminal row covers 2 pixel rows using Unicode half-blocks:
    //   ' ' = neither  '▀' = top only  '▄' = bottom only  '█' = both
    static const char *const kGlyph[4] = {" ", "\xe2\x96\x80", "\xe2\x96\x84", "\xe2\x96\x88"};

    fputs("+", stdout);
    for (int i = 0; i < BufferCanvas::W; i++) fputs("-", stdout);
    fputs("+\n", stdout);

    for (int y = 0; y < BufferCanvas::H; y += 2) {
        fputc('|', stdout);
        for (int x = 0; x < BufferCanvas::W; x++) {
            int top    = (fb[x + (y / 8)       * BufferCanvas::W] >> (y % 8))       & 1;
            int bottom = (fb[x + ((y + 1) / 8) * BufferCanvas::W] >> ((y + 1) % 8)) & 1;
            fputs(kGlyph[top | (bottom << 1)], stdout);
        }
        fputs("|\n", stdout);
    }

    fputs("+", stdout);
    for (int i = 0; i < BufferCanvas::W; i++) fputs("-", stdout);
    fputs("+\n", stdout);
}

static void show(BufferCanvas &canvas, I2cSsd1306 *oled, const char *label) {
    printf("\n=== %s ===\n", label);
    print_fb(canvas.framebuffer());
    if (oled) {
        try { oled->flush(canvas.framebuffer()); }
        catch (const std::exception &e) { fprintf(stderr, "OLED flush: %s\n", e.what()); }
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }
}

// Southern Cross center: RA ~12h26m (~186.5°), Dec ~-60°.
static constexpr float SC_RA  = 186.5f * float(M_PI) / 180.0f;
static constexpr float SC_DEC = -60.0f * float(M_PI) / 180.0f;
static constexpr float FOV    = 70.0f  * float(M_PI) / 180.0f;

// Build a PolarAligner with `n` samples evenly spread in RA at `dec_deg`
// from the south celestial pole, 1-second timestamps so sidereal derotation
// is negligible.  n=3 at -85° gives arc ~8.7° (<10°, "Rotate more!").
//                 n=3 at -80° gives arc ~17°  (>10°, valid estimate).
static PolarAligner make_pa(float dec_deg, int n) {
    PolarAligner pa(OBSERVER_LAT_DEG, OBSERVER_LON_DEG);
    float dec = dec_deg * float(M_PI) / 180.0f;
    for (int i = 0; i < n; i++) {
        float ra = i * 2.0f * float(M_PI) / float(n);
        pa.add_sample(ra, dec, double(i));
    }
    return pa;
}

int main(int argc, char *argv[]) {
    bool use_oled = false;
    std::string oled_bus = "/dev/i2c-1";
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--oled")                         use_oled  = true;
        if (std::string(argv[i]) == "--oled-bus" && i + 1 < argc) oled_bus = argv[++i];
    }

    BufferCanvas canvas;
    TrackerView  view;

    std::unique_ptr<I2cSsd1306> oled;
    if (use_oled) {
        try { oled = std::make_unique<I2cSsd1306>(oled_bus.c_str()); }
        catch (const std::exception &e) { fprintf(stderr, "OLED init: %s\n", e.what()); }
    }

    SolveResult solved{};
    solved.solved = true;
    solved.ra     = SC_RA;
    solved.dec    = SC_DEC;
    solved.fov    = FOV;
    solved.roll   = 0.0f;

    SolveResult no_solve{};
    no_solve.solved = false;

    PolarAligner pa0(OBSERVER_LAT_DEG, OBSERVER_LON_DEG);  // empty

    // Status screens
    view.render_status(canvas, "Starting");
    show(canvas, oled.get(), "status: Starting");

    view.render_status(canvas, "Source fail");
    show(canvas, oled.get(), "status: Source fail");

    // Tracking
    view.render(canvas, AppMode::Tracking, no_solve, pa0);
    show(canvas, oled.get(), "tracking: no solution");

    view.render(canvas, AppMode::Tracking, solved, pa0);
    show(canvas, oled.get(), "tracking: Southern Cross solved, FOV=70");

    // Polar-alignment sampling
    view.render(canvas, AppMode::PASampling, no_solve, pa0);
    show(canvas, oled.get(), "PA sampling: 0 samples, no solve");

    {
        auto pa = make_pa(-85.0f, 1);
        view.render(canvas, AppMode::PASampling, solved, pa);
        show(canvas, oled.get(), "PA sampling: 1 sample, solved");
    }
    {
        auto pa = make_pa(-85.0f, 3);  // arc ~8.7° -> "Rotate more!"
        view.render(canvas, AppMode::PASampling, solved, pa);
        show(canvas, oled.get(), "PA sampling: 3 samples, small arc (<10 deg)");
    }
    {
        auto pa = make_pa(-80.0f, 3);  // arc ~17° -> valid estimate shown
        view.render(canvas, AppMode::PASampling, solved, pa);
        show(canvas, oled.get(), "PA sampling: 3 samples, good arc (>10 deg)");
    }

    // Polar-alignment fix
    {
        auto pa = make_pa(-80.0f, 4);
        view.render(canvas, AppMode::PAFix, solved, pa);
        show(canvas, oled.get(), "PA fix: alt/az correction");
    }

    return 0;
}
