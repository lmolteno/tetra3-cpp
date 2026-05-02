#include "star_map.h"
#include "star_map_data.h"
#include <cmath>
#include <algorithm>

bool StarMapRenderer::project(float star_ra, float star_dec,
                               float center_ra, float center_dec,
                               float roll, float scale,
                               int cx, int cy, int map_h,
                               int &px, int &py) const {
    float dra = star_ra - center_ra;

    float sin_dec = std::sin(star_dec);
    float cos_dec = std::cos(star_dec);
    float sin_dec0 = std::sin(center_dec);
    float cos_dec0 = std::cos(center_dec);
    float cos_dra = std::cos(dra);
    float sin_dra = std::sin(dra);

    float cos_c = sin_dec0 * sin_dec + cos_dec0 * cos_dec * cos_dra;
    if (cos_c <= 0.01f) return false;

    float x = (cos_dec * sin_dra) / cos_c;
    float y = (cos_dec0 * sin_dec - sin_dec0 * cos_dec * cos_dra) / cos_c;

    if (roll != 0.0f) {
        float cr = std::cos(roll);
        float sr = std::sin(roll);
        float rx = x * cr - y * sr;
        float ry = x * sr + y * cr;
        x = rx;
        y = ry;
    }

    px = cx - static_cast<int>(x * scale + 0.5f);
    py = cy - static_cast<int>(y * scale + 0.5f);

    return true;
}

void StarMapRenderer::draw_line(IStarMapCanvas &canvas, int w,
                                 int x0, int y0, int x1, int y1,
                                 int clip_y0, int clip_y1, bool dashed) {
    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;
    int step = 0;

    while (true) {
        if (x0 >= 0 && x0 < w && y0 >= clip_y0 && y0 <= clip_y1) {
            if (!dashed || (step % 4 < 2)) {
                canvas.draw_pixel(x0, y0);
            }
        }

        if (x0 == x1 && y0 == y1) break;

        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x0 += sx; }
        if (e2 < dx) { err += dx; y0 += sy; }
        step++;
    }
}

void StarMapRenderer::render_impl(IStarMapCanvas &canvas, int w,
                                    float ra, float dec, float fov_rad,
                                    float roll, int map_y0, int map_y1) {
    int map_h = map_y1 - map_y0 + 1;
    int cx = w / 2;
    int cy = map_y0 + map_h / 2;

    float scale = static_cast<float>(w) / (2.0f * std::tan(fov_rad / 2.0f));

    struct VisibleStar {
        int idx;
        int px, py;
        int8_t mag10;
    };
    static VisibleStar visible[STAR_MAP_NUM_STARS];
    int num_visible = 0;

    static int star_visible_idx[STAR_MAP_NUM_STARS];

    // Use a wide margin so off-screen stars are still tracked for
    // constellation line drawing (lines are clipped during rasterisation).
    int margin = w * 2;

    for (int i = 0; i < STAR_MAP_NUM_STARS; i++) {
        star_visible_idx[i] = -1;
        int px, py;
        if (project(STAR_MAP_RA[i], STAR_MAP_DEC[i], ra, dec, roll, scale,
                    cx, cy, map_h, px, py)) {
            if (px >= -margin && px < w + margin &&
                py >= map_y0 - margin && py <= map_y1 + margin) {
                star_visible_idx[i] = num_visible;
                visible[num_visible] = {i, px, py, STAR_MAP_MAG10[i]};
                num_visible++;
            }
        }
    }

    // Constellation lines
    for (int i = 0; i < STAR_MAP_NUM_LINES; i++) {
        int i1 = STAR_MAP_LINES[i][0];
        int i2 = STAR_MAP_LINES[i][1];
        int vi1 = star_visible_idx[i1];
        int vi2 = star_visible_idx[i2];
        if (vi1 >= 0 && vi2 >= 0) {
            draw_line(canvas,  w,
                      visible[vi1].px, visible[vi1].py,
                      visible[vi2].px, visible[vi2].py,
                      map_y0, map_y1, true);
        }
    }

    // Stars: bright (mag < 3.5) = cross, dim = single pixel
    for (int i = 0; i < num_visible; i++) {
        int px = visible[i].px;
        int py = visible[i].py;
        int8_t m10 = visible[i].mag10;

        if (px < 0 || px >= w || py < map_y0 || py > map_y1) continue;

        canvas.draw_pixel(px, py);

        if (m10 < 35) {
            if (px > 0) canvas.draw_pixel(px - 1, py);
            if (px < w - 1) canvas.draw_pixel(px + 1, py);
            if (py > map_y0) canvas.draw_pixel(px, py - 1);
            if (py <= map_y1) canvas.draw_pixel(px, py + 1);
        }
    }

    // Center crosshair
    int ch_size = 3;
    for (int d = 1; d <= ch_size; d++) {
        if (cx - d >= 0)      canvas.draw_pixel(cx - d, cy);
        if (cx + d < w)       canvas.draw_pixel(cx + d, cy);
        if (cy - d >= map_y0) canvas.draw_pixel(cx, cy - d);
        if (cy + d <= map_y1) canvas.draw_pixel(cx, cy + d);
    }
}

void StarMapRenderer::render(IStarMapCanvas &canvas, int width, int height,
                              float ra, float dec, float fov_rad, float roll) {
    render_impl(canvas, width, ra, dec, fov_rad, roll, 0, height - 1);
}
