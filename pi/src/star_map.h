#pragma once

#include "display.h"
#include <cmath>

// Renders a star field onto an IDisplay using gnomonic projection.
// Shows constellation stick figure lines and star dots scaled by magnitude.
// Callback interface for star map rendering -- allows any resolution.
struct IStarMapCanvas {
    virtual ~IStarMapCanvas() = default;
    virtual void draw_pixel(int x, int y) = 0;
};

class StarMapRenderer {
public:
    // Render onto IDisplay (128x64, uses IDisplay::WIDTH)
    void render(IDisplay &disp, float ra, float dec, float fov_rad,
                float roll = 0.0f, int map_y0 = 0, int map_y1 = 63);

    // Render onto arbitrary-size canvas
    void render(IStarMapCanvas &canvas, int width, int height,
                float ra, float dec, float fov_rad, float roll = 0.0f);

private:
    bool project(float star_ra, float star_dec,
                 float center_ra, float center_dec,
                 float roll, float scale,
                 int cx, int cy, int map_h,
                 int &px, int &py) const;

    void draw_line(IStarMapCanvas &canvas, int w, int x0, int y0, int x1, int y1,
                   int clip_y0, int clip_y1, bool dashed = false);

    void render_impl(IStarMapCanvas &canvas, int w, float ra, float dec,
                     float fov_rad, float roll, int map_y0, int map_y1);
};
