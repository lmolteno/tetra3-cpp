#include "render/tracker_view.h"
#include "render/text.h"

#include <cmath>
#include <cstdio>

namespace {
// Adapts an ICanvas to the IStarMapCanvas single-method interface.
struct StarMapAdapter : IStarMapCanvas {
    ICanvas &canvas;
    explicit StarMapAdapter(ICanvas &c) : canvas(c) {}
    void draw_pixel(int x, int y) override { canvas.draw_pixel(x, y); }
};
}

void TrackerView::render(ICanvas &canvas, AppMode mode,
                         const SolveResult &result,
                         const PolarAligner &aligner,
                         const ImuHint *imu) {
    canvas.clear();

    switch (mode) {
    case AppMode::Tracking:
        draw_tracking(canvas, result, imu);
        break;
    case AppMode::PASampling: {
        auto pole = aligner.estimate_pole();
        draw_pa_sampling(canvas, aligner.num_samples(), pole, result.solved);
        break;
    }
    case AppMode::PAFix: {
        auto pole = aligner.estimate_pole();
        auto offset = aligner.pole_error(pole);
        draw_pa_fix(canvas, offset);
        break;
    }
    }

    canvas.update();
}

void TrackerView::render_status(ICanvas &canvas, const std::string &msg) {
    canvas.clear();
    text::draw(canvas, 2, 28, msg, /*large=*/true);
    canvas.update();
}

void TrackerView::draw_tracking(ICanvas &canvas, const SolveResult &result,
                                const ImuHint *imu) {
    if (result.solved) {
        // Show the central half of the camera FOV on the map.
        float map_fov = result.fov * 0.5f;
        StarMapAdapter adapter(canvas);
        star_map_.render(adapter, WIDTH, MAP_Y1 + 1,
                         result.ra, result.dec, map_fov, result.roll);
        draw_coord_bar(canvas, result.ra, result.dec);
    } else if (imu) {
        StarMapAdapter adapter(canvas);
        star_map_.render(adapter, WIDTH, MAP_Y1 + 1,
                         imu->ra, imu->dec, imu->fov * 0.5f, imu->roll);
        draw_imu_bar(canvas, *imu);
    } else {
        draw_no_solution(canvas);
    }
}

void TrackerView::draw_no_solution(ICanvas &canvas) {
    text::draw(canvas, 1, COORD_Y, "No solution");
}

void TrackerView::draw_imu_bar(ICanvas &canvas, const ImuHint &imu) {
    // Pre-calibration: the IMU pointing is a rough accel+mag fix without the
    // benefit of a plate solve, so don't pretend it's a precise RA/Dec. After
    // calibration we still tag it as IMU so the user knows it isn't fresh.
    const char *tag = imu.calibrated ? "IMU" : "IMU?";
    char buf[24];
    float ra_h = imu.ra * 12.0f / static_cast<float>(M_PI);
    if (ra_h < 0)      ra_h += 24.0f;
    if (ra_h >= 24.0f) ra_h -= 24.0f;
    float dec_deg = imu.dec * 180.0f / static_cast<float>(M_PI);
    snprintf(buf, sizeof(buf), "%s %4.1fh %+5.1f\xB0", tag, ra_h, dec_deg);
    text::draw(canvas, 1, COORD_Y, buf);
}

void TrackerView::draw_coord_bar(ICanvas &canvas, float ra_rad, float dec_rad) {
    float ra_deg  = ra_rad  * 180.0f / static_cast<float>(M_PI);
    float dec_deg = dec_rad * 180.0f / static_cast<float>(M_PI);

    float ra_hours = ra_deg / 15.0f;
    if (ra_hours < 0) ra_hours += 24.0f;
    int ra_h = static_cast<int>(ra_hours) % 24;
    int ra_m = static_cast<int>((ra_hours - static_cast<int>(ra_hours)) * 60.0f);
    int ra_s = static_cast<int>(((ra_hours - static_cast<int>(ra_hours)) * 60.0f - ra_m) * 60.0f);

    char dec_sign = dec_deg >= 0 ? '+' : '-';
    float abs_dec = std::fabs(dec_deg);
    int dec_d = static_cast<int>(abs_dec);
    int dec_m = static_cast<int>((abs_dec - dec_d) * 60.0f);
    int dec_s = static_cast<int>(((abs_dec - dec_d) * 60.0f - dec_m) * 60.0f);

    char buf[24];
    snprintf(buf, sizeof(buf), "%02dh%02dm%02ds %c%02d\xB0%02d'%02d\"",
             ra_h, ra_m, ra_s, dec_sign, dec_d, dec_m, dec_s);
    text::draw(canvas, 1, COORD_Y, buf);
}

void TrackerView::draw_pa_sampling(ICanvas &canvas, int num_samples,
                                    const PoleEstimate &pole, bool solved) {
    text::draw(canvas, 1, 0, "POLAR ALIGN");
    canvas.draw_hline(0, WIDTH - 1, 10);

    char buf[24];
    snprintf(buf, sizeof(buf), "Samples: %d", num_samples);
    text::draw(canvas, 1, 14, buf);

    if (!solved) {
        text::draw(canvas, 1, 24, "No solution");
    } else if (num_samples < 3) {
        text::draw(canvas, 1, 24, "Rotate mount in RA");
        text::draw(canvas, 1, 34, "then press SPACE");
    } else if (pole.valid) {
        float dec_deg = pole.dec * 180.0f / static_cast<float>(M_PI);
        char sign = dec_deg >= 0 ? '+' : '-';
        float ad = std::fabs(dec_deg);
        int d = static_cast<int>(ad);
        int m = static_cast<int>((ad - d) * 60.0f);
        snprintf(buf, sizeof(buf), "Pole: %c%d\xB0%02d'", sign, d, m);
        text::draw(canvas, 1, 24, buf);

        snprintf(buf, sizeof(buf), "Arc: %.1f\xB0", pole.arc_deg);
        text::draw(canvas, 1, 34, buf);

        if (pole.arc_deg < 10.0f) {
            text::draw(canvas, 1, 44, "Rotate more!");
        }
    }

    canvas.draw_hline(0, WIDTH - 1, 55);
    text::draw(canvas, 1, COORD_Y, "SPC:sample ESC:done");
}

void TrackerView::draw_pa_fix(ICanvas &canvas, const AltAzOffset &offset) {
    text::draw(canvas, 1, 0, "ADJUST MOUNT");
    canvas.draw_hline(0, WIDTH - 1, 10);

    char buf[24];
    char alt_sign = offset.alt_arcmin >= 0 ? '+' : '-';
    float alt_abs = std::fabs(offset.alt_arcmin);
    snprintf(buf, sizeof(buf), "Alt: %c%.1f'", alt_sign, alt_abs);
    text::draw(canvas, 1, 16, buf, /*large=*/true);

    char az_sign = offset.az_arcmin >= 0 ? '+' : '-';
    float az_abs = std::fabs(offset.az_arcmin);
    snprintf(buf, sizeof(buf), "Az:  %c%.1f'", az_sign, az_abs);
    text::draw(canvas, 1, 32, buf, /*large=*/true);

    canvas.draw_hline(0, WIDTH - 1, 46);
    snprintf(buf, sizeof(buf), "Total: %.1f'", offset.total_arcmin);
    text::draw(canvas, 1, 48, buf);

    canvas.draw_hline(0, WIDTH - 1, 55);
    text::draw(canvas, 1, COORD_Y, "SPC:resample ESC:ok");
}
