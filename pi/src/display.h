#pragma once

#include <tetra3/types.h>
#include "polar_align.h"
#include <string>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cstdio>

class IDisplay {
public:
    static constexpr int WIDTH = 128;
    static constexpr int HEIGHT = 64;
    static constexpr int FB_SIZE = WIDTH * HEIGHT / 8; // 1024 bytes

    // Star map area: y 0..54 (55 pixels)
    static constexpr int MAP_Y0 = 0;
    static constexpr int MAP_Y1 = 54;

    // Coordinate bar: y 55..63 (9 pixels, fits one line of 6x8 text)
    static constexpr int COORD_Y = 56;

    virtual ~IDisplay() = default;
    virtual void clear() = 0;
    virtual void draw_text(int x, int y, const std::string &text, bool large = false) = 0;
    virtual void draw_pixel(int x, int y) = 0;
    virtual void update() = 0;

    // Return raw 1-bit framebuffer (1024 bytes). nullptr if not available.
    virtual const uint8_t *get_framebuffer() const { return nullptr; }

    // Draw a horizontal line
    void draw_hline(int x0, int x1, int y) {
        for (int x = x0; x <= x1; x++) draw_pixel(x, y);
    }

    // Full-screen centered status message
    void show_status(const std::string &text) {
        clear();
        draw_text(2, 28, text, true);
        update();
    }

    // Draw the coordinate bar at the bottom.
    // Format: "05h36m12s  -01\xB012'00\""  (RA in HMS, Dec in DMS with degree symbol)
    void show_coordinates(float ra_rad, float dec_rad) {
        float ra_deg = ra_rad * 180.0f / M_PI;
        float dec_deg = dec_rad * 180.0f / M_PI;

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

        // "05h36m12s -01°12'00"" = 21 chars * 6px = 126px — fits exactly
        char buf[24];
        snprintf(buf, sizeof(buf), "%02dh%02dm%02ds %c%02d\xB0%02d'%02d\"",
                 ra_h, ra_m, ra_s, dec_sign, dec_d, dec_m, dec_s);
        draw_text(1, COORD_Y, buf);
    }

    // Show solve result: just the coordinate bar
    void show_solve_result(const SolveResult &r) {
        if (!r.solved) {
            draw_text(1, COORD_Y, "No solution");
            return;
        }
        show_coordinates(r.ra, r.dec);
    }

    // Polar alignment sampling screen
    void show_pa_sampling(int num_samples, const PoleEstimate &pole, bool solved) {
        draw_text(1, 0, "POLAR ALIGN");
        draw_hline(0, WIDTH - 1, 10);

        char buf[24];
        snprintf(buf, sizeof(buf), "Samples: %d", num_samples);
        draw_text(1, 14, buf);

        if (!solved) {
            draw_text(1, 24, "No solution");
        } else if (num_samples < 3) {
            draw_text(1, 24, "Rotate mount in RA");
            draw_text(1, 34, "then press SPACE");
        } else if (pole.valid) {
            float dec_deg = pole.dec * 180.0f / static_cast<float>(M_PI);
            char sign = dec_deg >= 0 ? '+' : '-';
            float ad = std::fabs(dec_deg);
            int d = static_cast<int>(ad);
            int m = static_cast<int>((ad - d) * 60.0f);
            snprintf(buf, sizeof(buf), "Pole: %c%d\xB0%02d'", sign, d, m);
            draw_text(1, 24, buf);

            snprintf(buf, sizeof(buf), "Arc: %.1f\xB0", pole.arc_deg);
            draw_text(1, 34, buf);

            if (pole.arc_deg < 10.0f) {
                draw_text(1, 44, "Rotate more!");
            }
        }

        draw_hline(0, WIDTH - 1, 55);
        draw_text(1, COORD_Y, "SPC:sample ESC:done");
    }

    // Polar alignment fix/correction screen
    void show_pa_fix(const AltAzOffset &offset) {
        draw_text(1, 0, "ADJUST MOUNT");
        draw_hline(0, WIDTH - 1, 10);

        char buf[24];
        char alt_sign = offset.alt_arcmin >= 0 ? '+' : '-';
        float alt_abs = std::fabs(offset.alt_arcmin);
        snprintf(buf, sizeof(buf), "Alt: %c%.1f'", alt_sign, alt_abs);
        draw_text(1, 16, buf, true);

        char az_sign = offset.az_arcmin >= 0 ? '+' : '-';
        float az_abs = std::fabs(offset.az_arcmin);
        snprintf(buf, sizeof(buf), "Az:  %c%.1f'", az_sign, az_abs);
        draw_text(1, 32, buf, true);

        draw_hline(0, WIDTH - 1, 46);
        snprintf(buf, sizeof(buf), "Total: %.1f'", offset.total_arcmin);
        draw_text(1, 48, buf);

        draw_hline(0, WIDTH - 1, 55);
        draw_text(1, COORD_Y, "SPC:resample ESC:ok");
    }
};
