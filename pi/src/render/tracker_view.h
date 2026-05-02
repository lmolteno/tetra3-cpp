#pragma once

#include "render/canvas.h"
#include "app/polar_align.h"
#include "render/star_map.h"
#include <tetra3/types.h>
#include <string>

// Composes the tracker UI on top of an ICanvas. Knows nothing about how
// pixels reach hardware. `render()` clears, draws, then calls update().
class TrackerView {
public:
    static constexpr int WIDTH  = 128;
    static constexpr int HEIGHT = 64;

    // Star map fills y 0..54; coordinate bar lives at y 55..63.
    static constexpr int MAP_Y0 = 0;
    static constexpr int MAP_Y1 = 54;
    static constexpr int COORD_Y = 56;

    void render(ICanvas &canvas,
                AppMode mode,
                const SolveResult &result,
                const PolarAligner &aligner);

    // Full-screen centered status (loading, errors).
    void render_status(ICanvas &canvas, const std::string &text);

private:
    StarMapRenderer star_map_;

    void draw_tracking(ICanvas &canvas, const SolveResult &result);
    void draw_pa_sampling(ICanvas &canvas, int num_samples,
                          const PoleEstimate &pole, bool solved);
    void draw_pa_fix(ICanvas &canvas, const AltAzOffset &offset);

    void draw_coord_bar(ICanvas &canvas, float ra_rad, float dec_rad);
    void draw_no_solution(ICanvas &canvas);
};
