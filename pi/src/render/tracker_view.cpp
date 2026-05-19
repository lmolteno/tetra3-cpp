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

// Bresenham line + small arrowhead. Used to draw "move the mount this way"
// from the chart centre to wherever the real pole is projecting.
void draw_line(ICanvas &c, int x0, int y0, int x1, int y1) {
    int dx = std::abs(x1 - x0), dy = std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
    int err = dx - dy;
    while (true) {
        c.draw_pixel(x0, y0);
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x0 += sx; }
        if (e2 <  dx) { err += dx; y0 += sy; }
    }
}

// 5-pixel arrowhead at (tx, ty), pointing along unit direction (ux, uy).
// (tx, ty) is the tip. The two "fins" sit ~3 px back along -direction and
// 2 px perpendicular.
void draw_arrowhead(ICanvas &c, int tx, int ty, float ux, float uy) {
    auto round = [](float v) { return static_cast<int>(v + (v >= 0 ? 0.5f : -0.5f)); };
    float px = -uy, py = ux;          // perpendicular
    int bx = tx - round(3.0f * ux);
    int by = ty - round(3.0f * uy);
    int lx = bx + round(2.0f * px), ly = by + round(2.0f * py);
    int rx = bx - round(2.0f * px), ry = by - round(2.0f * py);
    draw_line(c, tx, ty, lx, ly);
    draw_line(c, tx, ty, rx, ry);
}

}  // namespace

void TrackerView::draw_pa_fix(ICanvas &canvas, const AltAzOffset &offset,
                              const PAFixState *pa_fix) {
    // Camera-centred chart at a fixed 5° FOV. The star field is drawn around
    // the camera's actual pointing (with the camera's actual roll), so the
    // operator sees the patch of sky the camera is on. The real celestial
    // pole is drawn as a target marker; an arrow goes from chart centre
    // toward it. As alignment improves, the mount slews the camera toward
    // the pole, the marker walks in toward the centre, and the arrow
    // shortens. When the marker is at the centre the mount points at the
    // real pole — i.e. it's polar-aligned.
    constexpr int CHART_H = 48;                  // y: 0 .. 47
    constexpr int CHART_Y0 = 0;
    constexpr int CHART_Y1 = CHART_H - 1;
    constexpr float FOV_RAD =
        5.0f * static_cast<float>(M_PI) / 180.0f;

    // Default to south if no state — matches the OBSERVER_LAT_DEG compile
    // default. With a PAFixState this is unused except for picking the pole.
    float observer_lat = pa_fix ? pa_fix->observer_lat
                                : static_cast<float>(OBSERVER_LAT_DEG)
                                  * static_cast<float>(M_PI) / 180.0f;
    float true_pole_dec = (observer_lat >= 0.0f)
                        ?  static_cast<float>(M_PI / 2)
                        : -static_cast<float>(M_PI / 2);

    float cam_ra   = pa_fix ? pa_fix->cam_ra   : 0.0f;
    float cam_dec  = pa_fix ? pa_fix->cam_dec  : 0.0f;
    float cam_roll = pa_fix ? pa_fix->cam_roll : 0.0f;

    // Star field at the camera's current pointing. Camera roll = chart roll
    // here (chart is centred on the camera, no parallel-transport needed).
    PaCanvasAdapter adapter(canvas);
    star_map_.render(adapter, WIDTH, CHART_H,
                     cam_ra, cam_dec, FOV_RAD, cam_roll);

    // Project the real celestial pole into chart space.
    int cx = WIDTH / 2;
    int cy = CHART_Y0 + CHART_H / 2;
    const float scale =
        static_cast<float>(WIDTH) / (2.0f * std::tan(FOV_RAD / 2.0f));
    int px, py;
    bool projected = gnomonic(0.0f, true_pole_dec,
                              cam_ra, cam_dec, cam_roll, scale,
                              cx, cy, px, py);

    if (projected) {
        const int MARGIN = 4;
        bool on_chart = (px >= MARGIN && px < WIDTH - MARGIN &&
                         py >= CHART_Y0 + MARGIN && py <= CHART_Y1 - MARGIN);
        int tip_x = px, tip_y = py;
        if (!on_chart) {
            // Clamp to a point near the edge, on the line from centre to
            // the pole's projected position. Arrow lives at the edge.
            float dx = static_cast<float>(px - cx);
            float dy = static_cast<float>(py - cy);
            float len = std::sqrt(dx * dx + dy * dy);
            if (len < 1.0f) { dx = 1.0f; dy = 0.0f; len = 1.0f; }
            float ux = dx / len, uy = dy / len;
            int max_x = WIDTH / 2 - MARGIN;
            int max_y = CHART_H / 2 - MARGIN;
            float r = std::min(static_cast<float>(max_x) / std::fabs(ux + 1e-6f),
                                static_cast<float>(max_y) / std::fabs(uy + 1e-6f));
            tip_x = cx + static_cast<int>(ux * r + 0.5f);
            tip_y = cy + static_cast<int>(uy * r + 0.5f);
        }

        // Shaft from centre to tip + arrowhead at the tip. As alignment
        // improves the shaft gets shorter.
        if (std::abs(tip_x - cx) > 1 || std::abs(tip_y - cy) > 1) {
            draw_line(canvas, cx, cy, tip_x, tip_y);
            float dx = static_cast<float>(tip_x - cx);
            float dy = static_cast<float>(tip_y - cy);
            float len = std::sqrt(dx * dx + dy * dy);
            draw_arrowhead(canvas, tip_x, tip_y, dx / len, dy / len);
        }

        if (on_chart) {
            // Target marker: small hollow square at the projected pole.
            // (When the camera is exactly on the real pole, this lands on
            // top of the star-map crosshair and reinforces "you're here".)
            canvas.draw_pixel(px - 1, py - 1);
            canvas.draw_pixel(px,     py - 1);
            canvas.draw_pixel(px + 1, py - 1);
            canvas.draw_pixel(px - 1, py);
            canvas.draw_pixel(px + 1, py);
            canvas.draw_pixel(px - 1, py + 1);
            canvas.draw_pixel(px,     py + 1);
            canvas.draw_pixel(px + 1, py + 1);
        }
    }

    // Bottom strip: alt/az breakdown on one line, total on the next. All
    // small text so nothing clips off the bottom of the 64-px screen.
    char buf[32];
    snprintf(buf, sizeof(buf), "%+5.1f' %+5.1f'",
             offset.alt_arcmin, offset.az_arcmin);
    text::draw(canvas, 1, 49, buf);
    // Cap the total at the visible cell width so absurd values from a bad
    // fit don't spill into "SPC/ESC".
    if (offset.total_arcmin >= 999.5f)
        snprintf(buf, sizeof(buf), "T>999'");
    else
        snprintf(buf, sizeof(buf), "T%5.1f'", offset.total_arcmin);
    text::draw(canvas, 1, 57, buf);
    text::draw(canvas, 80, 57, "SPC/ESC");
}
