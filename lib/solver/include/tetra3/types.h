#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <cstdint>

struct Centroid {
    double y, x; // pixel coordinates (y=row, x=col)
};

struct StarEntry {
    float ra, dec; // radians
    float x, y, z; // unit vector components
    float magnitude;
};

struct PatternEntry {
    uint16_t star_indices[4]; // indices into star catalog
};

struct SolveResult {
    bool solved;
    float ra, dec, roll; // radians
    float fov; // horizontal FOV in radians
    float rmse; // RMS error in arcseconds
    int num_matches;
    float solve_time_ms;
    float distortion_k;
};

// Per-pair diagnostic emitted on demand by SimpleStarSolver::solve_from_centroids
// (pass a non-null vector*). Coordinates are in the original (un-cropped) image
// the caller passed in; the solver does NOT remap them.
struct MatchInfo {
    int catalog_index;             // index into the loaded star catalog
    double img_x, img_y;           // observed centroid (post-undistortion)
    double pred_x, pred_y;         // catalog star projected with the final solve
    double residual_pix;           // sqrt((img-pred)^2)
    float catalog_ra, catalog_dec; // radians, copied for convenience
    float magnitude;
};

inline long double combinations(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    if (k > n / 2) {
        k = n - k;
    }
    long double res = 1.0L;
    for (int i = 1; i <= k; ++i) {
        res = res * (n - i + 1) / i;
    }
    return res;
}

// Binomial CDF: P(X <= k) = sum_{i=0 to k} C(n, i) * p^i * (1-p)^(n-i)
inline long double calculate_binomial_cdf(int k, int n, long double p) {
    if (k < 0)
        return 0.0L;
    if (k >= n)
        return 1.0L;
    if (p < 0.0L || p > 1.0L)
        return 0.0L;

    long double cdf = 0.0L;
    long double one_minus_p = 1.0L - p;

    for (int i = 0; i <= k; ++i) {
        long double term = combinations(n, i);
        term *= std::pow(p, i);
        term *= std::pow(one_minus_p, n - i);
        cdf += term;
    }
    return cdf;
}
