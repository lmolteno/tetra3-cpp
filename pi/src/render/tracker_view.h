#pragma once

#include "render/canvas.h"
#include "app/polar_align.h"
#include "render/star_map.h"
#include <tetra3/types.h>
#include <string>

// IMU-derived pointing hint, passed in by main.cpp. Decoupled from the imu/
// module so oled_test (which doesn't link the IMU code) can still compile.
struct ImuHint {
    float ra;     // radians
    float dec;    // radians
    float roll;   // radians
    float fov;    // radians — what FOV to draw the star map at
    bool  calibrated;
};

// Live PAFix state. main.cpp fills this in while in PAFix mode; the view
// uses it to draw the chart and to display live-updating alt/az numbers as
// the operator turns the alt/az knobs (no new fit is required).
//
//   live_pole    : live mount-pole estimate, from PAFixTracker
//   live_offset  : alt/az error of the live pole vs the true celestial pole
//   cam_ra/dec   : current camera direction (from latest plate solve)
//   cam_roll     : current camera roll — rotates the chart so the alt/az
//                  motion matches what the operator sees through the camera
//   observer_lat : observer latitude (radians); only used to pick north vs
//                  south celestial pole for the chart target
struct PAFixState {
    PoleEstimate live_pole;
    AltAzOffset  live_offset;
    float cam_ra      = 0;
    float cam_dec     = 0;
    float cam_roll    = 0;
    float observer_lat = 0;
};

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

    // imu: optional IMU-derived pointing. If `result` is unsolved we fall
    // back to drawing the star map at the IMU's pointing (and skip the
    // coord bar — the IMU estimate is rough). Pass nullptr to disable.
    // pa_fix: live PAFix state (only consulted in PAFix mode). If absent,
    // PAFix mode falls back to drawing the (frozen-since-PASampling) fit.
    void render(ICanvas &canvas,
                AppMode mode,
                const SolveResult &result,
                const PolarAligner &aligner,
                const ImuHint *imu = nullptr,
                const PAFixState *pa_fix = nullptr);

    // Full-screen centered status (loading, errors).
    void render_status(ICanvas &canvas, const std::string &text);

private:
    StarMapRenderer star_map_;

    void draw_tracking(ICanvas &canvas, const SolveResult &result,
                       const ImuHint *imu);
    void draw_pa_sampling(ICanvas &canvas, int num_samples,
                          const PoleEstimate &pole, bool solved);
    void draw_pa_fix(ICanvas &canvas, const AltAzOffset &offset,
                     const PoleEstimate &mount_pole,
                     const PAFixState *pa_fix);

    void draw_coord_bar(ICanvas &canvas, float ra_rad, float dec_rad);
    void draw_imu_bar(ICanvas &canvas, const ImuHint &imu);
    void draw_no_solution(ICanvas &canvas);
};
