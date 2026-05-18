#include "render/tracker_view.h"
#include "render/text.h"

#include <algorithm>
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
                         const ImuHint *imu,
                         const PAFixState *pa_fix) {
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
        // Prefer the live (tracker-updated) offset; fall back to the frozen
        // fit if the tracker hasn't been initialised yet.
        AltAzOffset offset = (pa_fix && pa_fix->live_pole.valid)
                           ? pa_fix->live_offset
                           : aligner.pole_error(aligner.estimate_pole());
        draw_pa_fix(canvas, offset, pa_fix);
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

namespace {

// Adapter used to overlay markers on top of a star map rendered into an
// ICanvas (we only need draw_pixel from the marker side).
struct PaCanvasAdapter : IStarMapCanvas {
    ICanvas &canvas;
    explicit PaCanvasAdapter(ICanvas &c) : canvas(c) {}
    void draw_pixel(int x, int y) override { canvas.draw_pixel(x, y); }
};

// Gnomonic projection — identical to StarMapRenderer::project, replicated
// here so we can reuse the *same* scale + roll the star map was just drawn
// with and overlay a marker at the matching pixel.
bool gnomonic(float ra, float dec,
              float center_ra, float center_dec,
              float roll, float scale,
              int cx, int cy, int &px, int &py) {
    float dra = ra - center_ra;
    float sd  = std::sin(dec),     cd  = std::cos(dec);
    float sd0 = std::sin(center_dec), cd0 = std::cos(center_dec);
    float cdra = std::cos(dra),    sdra = std::sin(dra);
    float cosc = sd0 * sd + cd0 * cd * cdra;
    if (cosc <= 0.01f) return false;
    float x = (cd * sdra) / cosc;
    float y = (cd0 * sd - sd0 * cd * cdra) / cosc;
    if (roll != 0.0f) {
        float cr = std::cos(roll), sr = std::sin(roll);
        float rx = x * cr - y * sr;
        float ry = x * sr + y * cr;
        x = rx; y = ry;
    }
    px = cx - static_cast<int>(x * scale + 0.5f);
    py = cy - static_cast<int>(y * scale + 0.5f);
    return true;
}

// Continuous autoscale for the pole-centered chart. We want the live
// mount-pole marker to sit at roughly 1/4 of the chart radius from centre
// regardless of how big the offset is, so FOV is set proportional to the
// offset itself (FOV ≈ 8·offset → marker at 64/8 ≈ 8 px from centre, about
// 36% of the 22-px vertical chart half-height). Clamped at the small end so
// the star map still has stars to anchor against (Octans / UMi are sparse
// below ~5'), and at the large end so a wildly-bad initial fit doesn't blow
// the projection up.
float autoscale_fov_rad(float total_arcmin) {
    constexpr float ARCMIN = static_cast<float>(M_PI) / (180.0f * 60.0f);
    float f = 8.0f * total_arcmin;
    if (f <    5.0f) f =    5.0f;   //  min 5'   FOV — keep stars in frame
    if (f > 3600.0f) f = 3600.0f;   //  max 60°  FOV
    return f * ARCMIN;
}

void format_fov_label(float fov_rad, char *out, size_t n) {
    float arcmin = fov_rad * (180.0f * 60.0f / static_cast<float>(M_PI));
    if (arcmin < 60.0f) snprintf(out, n, "%.0f'", arcmin);
    else                snprintf(out, n, "%.1f\xB0", arcmin / 60.0f);
}

}  // namespace

void TrackerView::draw_pa_fix(ICanvas &canvas, const AltAzOffset &offset,
                              const PAFixState *pa_fix) {
    // Layout: pole-centered star map fills the top of the screen, the live
    // mount-pole marker sits on top of it, and the numerical offset is shown
    // below. The map uses the SAME projection (gnomonic via StarMapRenderer)
    // and the SAME camera roll, so the orientation matches what the operator
    // sees through the camera — they can correlate stars on-screen with what
    // they see in the camera view, regardless of where the camera is pointed.
    const int map_h = 45;             // y: 0 .. 44
    const int map_y0 = 0, map_y1 = map_h - 1;

    // Pick the true celestial pole — north or south based on observer lat.
    // Fall back to south if we have no PAFixState (compile-time default at
    // OBSERVER_LAT_DEG = -33.86 → southern hemisphere).
    float observer_lat = pa_fix ? pa_fix->observer_lat
                                : static_cast<float>(OBSERVER_LAT_DEG)
                                  * static_cast<float>(M_PI) / 180.0f;
    float true_pole_dec = (observer_lat >= 0.0f)
                        ?  static_cast<float>(M_PI / 2)
                        : -static_cast<float>(M_PI / 2);

    // Autoscale FOV so the mount-pole marker doesn't get buried. As the
    // operator drives the offset down, the chart zooms in.
    float fov = autoscale_fov_rad(offset.total_arcmin);
    float cam_roll = pa_fix ? pa_fix->cam_roll : 0.0f;

    // Render stars + lines + center crosshair around the celestial pole. The
    // crosshair the renderer drops at the centre IS the true pole.
    PaCanvasAdapter adapter(canvas);
    star_map_.render(adapter, WIDTH, map_h,
                     /*ra=*/0.0f, true_pole_dec, fov, cam_roll);

    // Overlay the live mount-pole marker using the same scale the star map
    // used internally. scale_pix_per_rad = width / (2 * tan(fov/2)).
    if (pa_fix && pa_fix->live_pole.valid) {
        float scale = static_cast<float>(WIDTH) /
                      (2.0f * std::tan(fov / 2.0f));
        int cx = WIDTH / 2;
        int cy = map_y0 + map_h / 2;
        int mx, my;
        if (gnomonic(pa_fix->live_pole.ra, pa_fix->live_pole.dec,
                     0.0f, true_pole_dec, cam_roll, scale,
                     cx, cy, mx, my)) {
            if (mx >= 0 && mx < WIDTH && my >= map_y0 && my <= map_y1) {
                // Hollow square (3x3) so the marker stands out from the
                // single-pixel star dots: corners + middle of each side.
                canvas.draw_pixel(mx - 1, my - 1);
                canvas.draw_pixel(mx + 1, my - 1);
                canvas.draw_pixel(mx - 1, my + 1);
                canvas.draw_pixel(mx + 1, my + 1);
                canvas.draw_pixel(mx,     my - 1);
                canvas.draw_pixel(mx,     my + 1);
                canvas.draw_pixel(mx - 1, my);
                canvas.draw_pixel(mx + 1, my);
            } else {
                // Off-screen at this zoom — draw an edge arrow pointing at it.
                int ex = std::clamp(mx, 1, WIDTH - 2);
                int ey = std::clamp(my, map_y0 + 1, map_y1 - 1);
                canvas.draw_pixel(ex, ey);
                canvas.draw_pixel(ex - 1, ey);
                canvas.draw_pixel(ex + 1, ey);
                canvas.draw_pixel(ex, ey - 1);
                canvas.draw_pixel(ex, ey + 1);
            }
        }
    }

    // Scale label (left) + numerical offset (right).
    char buf[24];
    char lbl[12];
    format_fov_label(fov, lbl, sizeof(lbl));
    snprintf(buf, sizeof(buf), "FOV %s", lbl);
    text::draw(canvas, 1, 47, buf);
    snprintf(buf, sizeof(buf), "%+.1f' %+.1f' T%.1f'",
             offset.alt_arcmin, offset.az_arcmin, offset.total_arcmin);
    // 22 chars max at 6 px per glyph = 132 px → too wide. Compact form:
    snprintf(buf, sizeof(buf), "%+.1f'/%+.1f'", offset.alt_arcmin,
             offset.az_arcmin);
    text::draw(canvas, 56, 47, buf);
    snprintf(buf, sizeof(buf), "T%.1f'", offset.total_arcmin);
    text::draw(canvas, 1, 56, buf, /*large=*/true);

    text::draw(canvas, 80, COORD_Y, "SPC/ESC");
}
