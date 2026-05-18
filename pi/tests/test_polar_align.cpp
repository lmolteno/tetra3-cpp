// Synthetic-sample tests for the polar-alignment pipeline.
//
// Generates samples around a known fake mount pole (with optional sidereal
// drift and per-sample noise), feeds them to PolarAligner, and checks that:
//   - estimate_pole() recovers the true pole tightly
//   - pole_error() reports near-zero when the mount is well-aligned
//   - pole_error() reports a non-trivial offset when the mount is misaligned
//   - the estimator stays sane under realistic per-sample noise
//
// No test framework — just asserts with a small CHECK macro.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "app/polar_align.h"

namespace {

constexpr double DEG = M_PI / 180.0;
constexpr double SIDEREAL_RATE = 7.2921159e-5;  // rad/s — must match polar_align.cpp

struct V3 { double x, y, z; };

V3 sph_to_vec(double ra, double dec) {
    double cd = std::cos(dec);
    return { cd * std::cos(ra), cd * std::sin(ra), std::sin(dec) };
}

void vec_to_sph(const V3 &v, double &ra, double &dec) {
    dec = std::asin(std::clamp(v.z, -1.0, 1.0));
    ra  = std::atan2(v.y, v.x);
    if (ra < 0) ra += 2.0 * M_PI;
}

// Rotate v around unit axis k by angle theta (Rodrigues).
V3 rotate_axis(const V3 &v, const V3 &k, double theta) {
    double c = std::cos(theta), s = std::sin(theta);
    double dot = k.x * v.x + k.y * v.y + k.z * v.z;
    V3 cross = { k.y * v.z - k.z * v.y,
                 k.z * v.x - k.x * v.z,
                 k.x * v.y - k.y * v.x };
    return {
        v.x * c + cross.x * s + k.x * dot * (1.0 - c),
        v.y * c + cross.y * s + k.y * dot * (1.0 - c),
        v.z * c + cross.z * s + k.z * dot * (1.0 - c),
    };
}

double angle_between(const V3 &a, const V3 &b) {
    double d = std::clamp(a.x * b.x + a.y * b.y + a.z * b.z, -1.0, 1.0);
    return std::acos(d);
}

// Synthesize a run of (ra, dec, t) samples that look like a star field
// observed as the mount rotates around `pole_*` by `total_arc_rad`, while the
// sky drifts at the sidereal rate (PolarAligner is expected to de-rotate it).
std::vector<PASample> synth_samples(
    double target_ra, double target_dec,
    double pole_ra,   double pole_dec,
    double total_arc_rad,
    int n,
    double dt_seconds,
    double noise_arcsec = 0.0,
    unsigned seed = 42)
{
    std::vector<PASample> out;
    out.reserve(n);
    V3 target = sph_to_vec(target_ra, target_dec);
    V3 pole   = sph_to_vec(pole_ra, pole_dec);

    std::mt19937 rng(seed);
    const double noise_rad = noise_arcsec * (M_PI / (180.0 * 3600.0));
    std::normal_distribution<double> nz(0.0, noise_rad);

    const double t0 = 1700000000.0;  // arbitrary unix epoch
    for (int i = 0; i < n; ++i) {
        double frac  = (n == 1) ? 0.0 : double(i) / double(n - 1);
        double theta = frac * total_arc_rad;
        V3 rotated   = rotate_axis(target, pole, theta);

        double ra, dec;
        vec_to_sph(rotated, ra, dec);

        double t = t0 + i * dt_seconds;
        // Add sidereal drift — PolarAligner removes this internally.
        ra += SIDEREAL_RATE * (t - t0);

        if (noise_arcsec > 0) {
            ra  += nz(rng);
            dec += nz(rng);
        }
        out.push_back({ static_cast<float>(ra),
                        static_cast<float>(dec),
                        t });
    }
    return out;
}

int failures = 0;

void report_ok(const char *msg) {
    std::printf("  ok    %s\n", msg);
}
void report_fail(const char *msg, const char *file, int line) {
    std::printf("  FAIL  %s  (%s:%d)\n", msg, file, line);
    ++failures;
}

#define CHECK(cond, msg) do { \
    if (cond) report_ok(msg); else report_fail(msg, __FILE__, __LINE__); \
} while (0)

#define CHECK_NEAR(actual, expect, tol, msg) do { \
    double _a = (actual), _e = (expect), _t = (tol); \
    char buf[256]; \
    std::snprintf(buf, sizeof(buf), "%s [%.4g ~ %.4g ± %.4g]", msg, _a, _e, _t); \
    if (std::fabs(_a - _e) <= _t) report_ok(buf); \
    else report_fail(buf, __FILE__, __LINE__); \
} while (0)

void feed(PolarAligner &a, const std::vector<PASample> &s) {
    for (auto &x : s) a.add_sample(x.ra, x.dec, x.timestamp);
}

double pole_error_arcsec(float true_ra, float true_dec, const PoleEstimate &est) {
    return angle_between(sph_to_vec(true_ra, true_dec),
                         sph_to_vec(est.ra, est.dec))
         * (180.0 * 3600.0 / M_PI);
}

// ----------------------------------------------------------------------------
// Tests

void test_insufficient_samples_invalid() {
    std::puts("test: <3 samples returns invalid estimate");
    PolarAligner a(-33.86f, 151.21f);
    a.add_sample(0.0f, -1.5f, 0.0);
    a.add_sample(0.1f, -1.5f, 10.0);
    auto est = a.estimate_pole();
    CHECK(!est.valid, "estimate marked invalid");
}

void test_perfect_pole_recovery_southern() {
    std::puts("test: perfect pole recovery, southern hemisphere, 30' from SCP");
    const double pole_ra  =  90.0 * DEG;
    const double pole_dec = -89.5 * DEG;   // 30' from south celestial pole
    auto samples = synth_samples(80*DEG, -85*DEG, pole_ra, pole_dec,
                                 /*arc=*/90*DEG, /*n=*/8, /*dt=*/15.0);
    PolarAligner a(-33.86f, 151.21f);
    feed(a, samples);

    auto est = a.estimate_pole();
    CHECK(est.valid, "estimate valid");
    // arc_deg is the on-sky span, not the rotation angle. For target ~5°
    // from the pole rotated by 90°, the chord is acos(cos²5° + sin²5°·0) ≈ 10°.
    CHECK(est.arc_deg > 5.0 && est.arc_deg < 15.0, "on-sky arc in [5°, 15°]");
    CHECK_NEAR(pole_error_arcsec(pole_ra, pole_dec, est), 0.0, 1.0,
               "pole recovered within 1 arcsec");
}

void test_perfect_pole_recovery_northern() {
    std::puts("test: perfect pole recovery, northern hemisphere");
    const double pole_ra  =  60.0 * DEG;
    const double pole_dec =  89.2 * DEG;   // ~48' from NCP
    auto samples = synth_samples(45*DEG, 80*DEG, pole_ra, pole_dec,
                                 /*arc=*/60*DEG, /*n=*/6, /*dt=*/20.0);
    PolarAligner a(40.0f, -74.0f);   // approximately NYC
    feed(a, samples);

    auto est = a.estimate_pole();
    CHECK(est.valid, "estimate valid");
    CHECK_NEAR(pole_error_arcsec(pole_ra, pole_dec, est), 0.0, 1.0,
               "pole recovered within 1 arcsec");
}

void test_well_aligned_reports_near_zero() {
    std::puts("test: mount pole == true SCP -> pole_error near zero");
    auto samples = synth_samples(45*DEG, -80*DEG,
                                 /*pole=*/0.0, -90*DEG,
                                 60*DEG, 8, 15.0);
    PolarAligner a(-33.86f, 151.21f);
    feed(a, samples);
    auto est = a.estimate_pole();
    auto err = a.pole_error(est);
    CHECK(est.valid, "estimate valid");
    CHECK_NEAR(err.total_arcmin, 0.0, 1.0, "alt/az error < 1 arcmin");
}

void test_misaligned_reports_offset() {
    std::puts("test: 30' mount offset -> pole_error reports non-trivial alt/az");
    const double pole_ra  =  90.0 * DEG;
    const double pole_dec = -89.5 * DEG;  // 30' great-circle from SCP
    auto samples = synth_samples(80*DEG, -85*DEG, pole_ra, pole_dec,
                                 90*DEG, 8, 15.0);
    PolarAligner a(-33.86f, 151.21f);
    feed(a, samples);
    auto est = a.estimate_pole();
    auto err = a.pole_error(est);
    // pole_error's "total" uses raw (alt, az) deltas (az is not scaled by
    // cos(alt)). At observer_lat=-33.86 the az-direction component is
    // inflated by up to 1/cos(33.86°) ≈ 1.20 — so a 30' great-circle offset
    // can land anywhere in roughly [30', 36']. Be permissive.
    CHECK(err.total_arcmin > 20.0 && err.total_arcmin < 45.0,
          "total alt/az error in [20, 45] arcmin for a 30' great-circle offset");
}

// ----------------------------------------------------------------------------
// PAFix-mode regression tests.
//
// In `pi_tracker`'s state machine, samples accumulate in BOTH PASampling and
// PAFix (see pi/src/main.cpp around the `aligner.add_sample(...)` call). The
// tests below exercise what happens to the alt/az feedback while in PAFix:
//   1. The intended path: user clears, rotates, adjusts mount, clears again,
//      rotates again -> reported total_arcmin shrinks. This is the loop the
//      operator is supposed to drive by hand.
//   2. Failure mode reported in field testing: while sitting in PAFix with
//      a stationary mount, the camera keeps emitting solves that all land
//      at the same (de-rotated) sky point. The spherical circle fit becomes
//      degenerate (D^T*D has a near-zero smallest eigenvalue) and the
//      reported alt/az numbers grow unboundedly.
//   3. Mixed-arc corruption: continuing to add samples after the operator
//      physically adjusts the mount mixes two small-circles in one fit.

void test_pafix_clean_cycle_reduces_error() {
    std::puts("test: clear+resample after mount adjustment reduces alt/az error");
    PolarAligner a(-33.86f, 151.21f);

    // Round 1: mount 60' from SCP.
    {
        auto s = synth_samples(80*DEG, -85*DEG,
                               /*pole=*/45*DEG, -89.0*DEG,
                               90*DEG, 8, 15.0);
        feed(a, s);
    }
    auto est1 = a.estimate_pole();
    auto err1 = a.pole_error(est1);

    // Operator physically adjusts the mount toward SCP; aligner is reset
    // by re-entering PASampling. Round 2: now 10' from SCP.
    a.clear_samples();
    {
        auto s = synth_samples(80*DEG, -85*DEG,
                               /*pole=*/45*DEG, -89.833*DEG,  // ~10'
                               90*DEG, 8, 15.0);
        feed(a, s);
    }
    auto est2 = a.estimate_pole();
    auto err2 = a.pole_error(est2);

    CHECK(est1.valid && est2.valid, "both rounds produce valid estimates");
    char buf[160];
    std::snprintf(buf, sizeof(buf),
                  "feedback shrinks: %.2f' -> %.2f'",
                  err1.total_arcmin, err2.total_arcmin);
    CHECK(err2.total_arcmin < err1.total_arcmin * 0.3, buf);
}

void test_pafix_stationary_camera_is_harmless() {
    std::puts("test: stationary camera in PAFix does NOT degrade the fit");
    PolarAligner a(-33.86f, 151.21f);

    // Good initial arc — 8 samples spanning 90° around a mount pole 30'
    // from SCP. This is what the user collected in PASampling.
    auto good = synth_samples(80*DEG, -85*DEG, 90*DEG, -89.5*DEG,
                              90*DEG, 8, 15.0);
    feed(a, good);
    auto err_good = a.pole_error(a.estimate_pole());

    // Camera then sits still for ~80 seconds. The mount isn't rotating, so
    // the same field is solved repeatedly; the only motion in (ra, dec) is
    // sidereal, which the aligner de-rotates out. Add 80 such samples with
    // ~5″ per-sample plate-solve jitter — the magnitude the field test saw.
    std::mt19937 rng(123);
    const double jitter_rad = 5.0 * (M_PI / (180.0 * 3600.0));
    std::normal_distribution<double> nz(0.0, jitter_rad);
    PASample last{ good.back().ra, good.back().dec, good.back().timestamp };
    for (int i = 1; i <= 80; ++i) {
        double t  = last.timestamp + i * 1.0;
        double ra = last.ra + SIDEREAL_RATE * (t - last.timestamp) + nz(rng);
        double dec = last.dec + nz(rng);
        a.add_sample(static_cast<float>(ra), static_cast<float>(dec), t);
    }
    auto err_after = a.pole_error(a.estimate_pole());

    char buf[200];
    std::snprintf(buf, sizeof(buf),
                  "stationary PAFix is stable: total %.2f' -> %.2f' (Δ=%.2f')",
                  err_good.total_arcmin, err_after.total_arcmin,
                  std::fabs(err_after.total_arcmin - err_good.total_arcmin));
    // Stationary clustered samples contribute near-zero difference vectors
    // to D^T*D, so the original good samples continue to dominate the fit.
    // We assert <1' of drift over the 80-sample run.
    CHECK(std::fabs(err_after.total_arcmin - err_good.total_arcmin) < 1.0, buf);
}

void test_pafix_single_outlier_solve_corrupts_fit() {
    std::puts("test: a single bad solve in PAFix shifts the pole significantly");
    PolarAligner a(-33.86f, 151.21f);

    auto good = synth_samples(80*DEG, -85*DEG, 90*DEG, -89.5*DEG,
                              90*DEG, 8, 15.0);
    feed(a, good);
    auto est_good = a.estimate_pole();
    double err_good_arcsec = pole_error_arcsec(90*DEG, -89.5*DEG, est_good);

    // Inject ONE bogus solve — wildly off, the kind of false-positive a
    // plate solver occasionally returns from a noisy frame.
    double t_bad = good.back().timestamp + 30.0;
    a.add_sample(static_cast<float>(200*DEG),
                 static_cast<float>(-70*DEG),
                 t_bad);

    auto est_bad = a.estimate_pole();
    double err_bad_arcsec = pole_error_arcsec(90*DEG, -89.5*DEG, est_bad);

    char buf[200];
    std::snprintf(buf, sizeof(buf),
                  "one outlier shifts pole: %.2f' -> %.2f'",
                  err_good_arcsec / 60.0, err_bad_arcsec / 60.0);
    // A single false-positive moves the estimate by hundreds of arcminutes;
    // this is what an operator would experience as "the numbers blowing out".
    CHECK(err_bad_arcsec > err_good_arcsec + 60.0 * 60.0, buf);  // >1° shift
}

void test_pafix_marginal_initial_samples_are_unstable() {
    std::puts("test: PAFix entered with just 3 close-spaced samples is unstable");
    PolarAligner a(-33.86f, 151.21f);

    // Enter PAFix with the minimum 3 samples and a very short rotation arc
    // (operator was impatient). The circle fit is barely determined.
    auto sparse = synth_samples(80*DEG, -85*DEG, 90*DEG, -89.5*DEG,
                                /*arc=*/5*DEG, /*n=*/3, /*dt=*/15.0);
    feed(a, sparse);
    auto est_init = a.estimate_pole();
    double init_arcmin = pole_error_arcsec(90*DEG, -89.5*DEG, est_init) / 60.0;

    // Now sit still in PAFix with 5″ plate-solve jitter. Because the
    // initial fit is barely constrained, jitter has more leverage.
    std::mt19937 rng(7);
    const double jitter_rad = 5.0 * (M_PI / (180.0 * 3600.0));
    std::normal_distribution<double> nz(0.0, jitter_rad);
    PASample last{ sparse.back().ra, sparse.back().dec, sparse.back().timestamp };
    for (int i = 1; i <= 50; ++i) {
        double t  = last.timestamp + i * 1.0;
        double ra = last.ra + SIDEREAL_RATE * (t - last.timestamp) + nz(rng);
        double dec = last.dec + nz(rng);
        a.add_sample(static_cast<float>(ra), static_cast<float>(dec), t);
    }
    auto est_final = a.estimate_pole();
    double final_arcmin = pole_error_arcsec(90*DEG, -89.5*DEG, est_final) / 60.0;

    char buf[200];
    std::snprintf(buf, sizeof(buf),
                  "marginal init + sitting: %.2f' -> %.2f' (informational)",
                  init_arcmin, final_arcmin);
    // Informational — we don't assert a hard bound, just that the fit
    // moves. The result depends on the RNG seed.
    CHECK(est_init.valid && est_final.valid, buf);
}

void test_pafix_adjust_mount_does_not_reduce_numbers() {
    std::puts("test: adjusting mount while sitting in PAFix does NOT shrink numbers");
    PolarAligner a(-33.86f, 151.21f);

    // 1. Tracking -> PASampling: get a clean 60° arc around the misaligned
    //    mount pole (60' from SCP). Enter PAFix at this point.
    const double old_pole_ra  = 90*DEG;
    const double old_pole_dec = -89.0*DEG;          // 60' from SCP
    auto arc = synth_samples(80*DEG, -85*DEG, old_pole_ra, old_pole_dec,
                             /*arc=*/60*DEG, /*n=*/8, /*dt=*/15.0);
    feed(a, arc);
    auto err_before = a.pole_error(a.estimate_pole());

    // 2. In PAFix the operator turns the alt/az knobs to correct the mount.
    //    Mechanically this both moves the mount pole toward SCP AND drags
    //    the camera pointing by roughly the same angular amount on the sky.
    //    Model: the camera now sits at the last arc sample shifted by ~50'
    //    in the alt direction (a representative "halfway-corrected" state).
    //    Sidereal drift continues; the aligner de-rotates it out.
    const double new_ra  = arc.back().ra;
    const double new_dec = arc.back().dec + (50.0 / 60.0) * DEG;
    std::mt19937 rng(11);
    const double jitter_rad = 5.0 * (M_PI / (180.0 * 3600.0));
    std::normal_distribution<double> nz(0.0, jitter_rad);
    const double t0 = arc.back().timestamp;
    for (int i = 1; i <= 30; ++i) {
        double t   = t0 + i * 1.0;                  // 1 Hz solves
        double ra  = new_ra  + SIDEREAL_RATE * (t - t0) + nz(rng);
        double dec = new_dec + nz(rng);
        a.add_sample(static_cast<float>(ra), static_cast<float>(dec), t);
    }
    auto err_after = a.pole_error(a.estimate_pole());

    char buf[200];
    std::snprintf(buf, sizeof(buf),
                  "PAFix-then-adjust does not improve numbers: %.2f' -> %.2f'",
                  err_before.total_arcmin, err_after.total_arcmin);
    // The post-adjustment cluster lies OFF the original small circle, so it
    // pulls the fit toward itself instead of refining the same arc. We
    // assert the reported number did NOT meaningfully drop.
    CHECK(err_after.total_arcmin >= err_before.total_arcmin * 0.7, buf);
}

void test_pafix_mixed_arcs_corrupts_fit() {
    std::puts("test: mixing samples from two mount poles in one fit is bogus");
    PolarAligner a(-33.86f, 151.21f);

    // Operator gets a clean half-arc around the misaligned pole...
    auto arc1 = synth_samples(80*DEG, -85*DEG, 90*DEG, -89.5*DEG,
                              45*DEG, 4, 15.0);  // 4 samples, 45°
    feed(a, arc1);
    auto est_partial = a.estimate_pole();
    double err_partial_arcmin = pole_error_arcsec(90*DEG, -89.5*DEG, est_partial) / 60.0;

    // ...adjusts the mount (new pole much closer to SCP)...
    // ...and *without clearing*, keeps rotating: 4 more samples around the
    // corrected pole get appended to the same fit.
    auto arc2 = synth_samples(80*DEG, -85*DEG, 90*DEG, -89.95*DEG,
                              45*DEG, 4, 15.0);
    feed(a, arc2);
    auto est_mixed = a.estimate_pole();
    double err_mixed_old_arcmin = pole_error_arcsec(90*DEG, -89.5*DEG,  est_mixed) / 60.0;
    double err_mixed_new_arcmin = pole_error_arcsec(90*DEG, -89.95*DEG, est_mixed) / 60.0;

    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "mixed fit is far from both poles: partial=%.2f' "
                  "vs old=%.2f' new=%.2f'",
                  err_partial_arcmin, err_mixed_old_arcmin, err_mixed_new_arcmin);
    // The mixed-arc estimate should be measurably worse than the
    // partial-arc estimate against either of the two true poles. We pick a
    // lenient threshold (5x worse than the partial result) so the test
    // isn't sensitive to which side of the mixed fit lands closer.
    bool bad_against_old = err_mixed_old_arcmin > err_partial_arcmin * 5.0 + 1.0;
    bool bad_against_new = err_mixed_new_arcmin > err_partial_arcmin * 5.0 + 1.0;
    CHECK(bad_against_old || bad_against_new, buf);
}

void test_noisy_samples_still_close() {
    std::puts("test: 20 arcsec per-sample noise -> pole within 5'");
    const double pole_ra  =  45.0 * DEG;
    const double pole_dec = -89.0 * DEG;  // ~60' from SCP
    auto samples = synth_samples(70*DEG, -80*DEG, pole_ra, pole_dec,
                                 120*DEG, 12, 10.0, /*noise_arcsec=*/20.0);
    PolarAligner a(-33.86f, 151.21f);
    feed(a, samples);
    auto est = a.estimate_pole();
    CHECK(est.valid, "estimate valid");
    double err_arcmin = pole_error_arcsec(pole_ra, pole_dec, est) / 60.0;
    char buf[128];
    std::snprintf(buf, sizeof(buf), "pole within 5' (got %.2f')", err_arcmin);
    CHECK(err_arcmin < 5.0, buf);
}

void test_small_arc_recovers_pole() {
    std::puts("test: short 20° rotation arc still recovers the pole");
    const double pole_ra  =  90.0 * DEG;
    const double pole_dec = -89.5 * DEG;
    auto samples = synth_samples(80*DEG, -85*DEG, pole_ra, pole_dec,
                                 20*DEG, 6, 10.0);
    PolarAligner a(-33.86f, 151.21f);
    feed(a, samples);
    auto est = a.estimate_pole();
    CHECK(est.valid, "estimate valid");
    CHECK(est.arc_deg > 0.5 && est.arc_deg < 5.0,
          "on-sky arc small but nonzero for short rotation");
    // Tighter than the noisy test because we have no noise here, but looser
    // than the 90° case because short arcs are geometrically less stable.
    CHECK_NEAR(pole_error_arcsec(pole_ra, pole_dec, est), 0.0, 60.0,
               "pole recovered within 1 arcmin from a 20° arc");
}

// ----------------------------------------------------------------------------
// PAFixTracker — live pole update with no new fit.

void test_pafix_tracker_idle_returns_drifted_pole() {
    std::puts("test: PAFixTracker idle (no mount motion) drifts pole by sidereal");
    PAFixTracker pa;
    const float pole_ra = 100*DEG, pole_dec = -89.5f*DEG;
    const float cam_ra  = 80*DEG,  cam_dec  = -85*DEG;
    pa.start(pole_ra, pole_dec, cam_ra, cam_dec, /*t=*/0.0);

    // 120 s later, mount untouched: camera (ra, dec) appears to have advanced
    // in RA at sidereal rate. Feed that as "current camera".
    const double dt = 120.0;
    const double sid = PAFixTracker::SIDEREAL_RATE;
    auto est = pa.live_pole(cam_ra + static_cast<float>(sid * dt),
                            cam_dec, dt);
    CHECK(est.valid, "estimate valid");
    // Live pole should be the original pole, sidereally drifted by sid*dt.
    double expected_ra = pole_ra + sid * dt;
    double dra = std::fmod(static_cast<double>(est.ra) - expected_ra
                           + 4 * M_PI, 2 * M_PI);
    if (dra > M_PI) dra -= 2 * M_PI;
    CHECK_NEAR(dra * (180.0 * 3600.0 / M_PI), 0.0, 1.0,
               "live pole RA drifts by sidereal amount (within 1″)");
    CHECK_NEAR((est.dec - pole_dec) * (180.0 * 3600.0 / M_PI), 0.0, 1.0,
               "live pole Dec unchanged");
}

void test_pafix_tracker_detects_mount_adjustment() {
    std::puts("test: PAFixTracker recovers a small mount rotation");
    PAFixTracker pa;

    // Snapshot: mount pole 30' from SCP, camera 5° from it.
    const double pole_ra  =  60*DEG;
    const double pole_dec = -89.5*DEG;
    const double cam_ra   =  80*DEG;
    const double cam_dec  = -85*DEG;
    pa.start(static_cast<float>(pole_ra), static_cast<float>(pole_dec),
             static_cast<float>(cam_ra),  static_cast<float>(cam_dec), 0.0);

    // Operator rotates the whole mount by 30' about an arbitrary 3D axis
    // (simulating a combined alt/az knob adjustment). Apply that same
    // rotation to the camera and feed the new camera direction.
    const double dt = 30.0;
    // Choose an axis in the celestial frame that resembles a real horizontal-
    // axis rotation: e.g., axis = (1, 0, 0). Magnitude 30 arcmin.
    V3 axis = {1, 0, 0};
    double angle = 30.0 / 60.0 * DEG;

    auto rotate = [&](const V3 &v) {
        double c = std::cos(angle), s = std::sin(angle);
        V3 cr = { axis.y*v.z - axis.z*v.y,
                  axis.z*v.x - axis.x*v.z,
                  axis.x*v.y - axis.y*v.x };
        double d = axis.x*v.x + axis.y*v.y + axis.z*v.z;
        return V3{ v.x*c + cr.x*s + axis.x*d*(1-c),
                   v.y*c + cr.y*s + axis.y*d*(1-c),
                   v.z*c + cr.z*s + axis.z*d*(1-c) };
    };

    // Drift the entry frame to t=dt and then apply the operator rotation.
    const double sid = PAFixTracker::SIDEREAL_RATE;
    V3 cam_drifted = sph_to_vec(cam_ra + sid*dt, cam_dec);
    V3 pole_drifted = sph_to_vec(pole_ra + sid*dt, pole_dec);
    V3 cam_after  = rotate(cam_drifted);
    V3 pole_after = rotate(pole_drifted);

    double new_cam_ra, new_cam_dec;
    vec_to_sph(cam_after, new_cam_ra, new_cam_dec);

    auto est = pa.live_pole(static_cast<float>(new_cam_ra),
                            static_cast<float>(new_cam_dec), dt);
    CHECK(est.valid, "live estimate valid");

    // The recovered live pole has a small irreducible residual: any
    // component of the real mount rotation along the camera direction is a
    // twist that doesn't move the camera, so we can't measure it. For a
    // camera 5° from the pole and a 30' rotation, that residual is up to
    // sin(5°)·(axis·cam_unit)·30' ≈ a few arcseconds. Real alt/az knob
    // rotations during PA are perpendicular to the (near-pole) camera, so
    // the residual stays in this range.
    double err_arcsec = angle_between(sph_to_vec(est.ra, est.dec), pole_after)
                      * (180.0 * 3600.0 / M_PI);
    char buf[160];
    std::snprintf(buf, sizeof(buf),
                  "live pole tracks operator rotation (err %.2f″)", err_arcsec);
    CHECK(err_arcsec < 5.0, buf);
}

void test_pafix_tracker_inactive_returns_invalid() {
    std::puts("test: PAFixTracker not started -> invalid live pole");
    PAFixTracker pa;
    auto est = pa.live_pole(0.0f, 0.0f, 0.0);
    CHECK(!est.valid, "inactive tracker returns invalid pole");
}

}  // namespace

int main() {
    test_insufficient_samples_invalid();
    test_perfect_pole_recovery_southern();
    test_perfect_pole_recovery_northern();
    test_well_aligned_reports_near_zero();
    test_misaligned_reports_offset();
    test_small_arc_recovers_pole();
    test_noisy_samples_still_close();
    test_pafix_clean_cycle_reduces_error();
    test_pafix_stationary_camera_is_harmless();
    test_pafix_single_outlier_solve_corrupts_fit();
    test_pafix_marginal_initial_samples_are_unstable();
    test_pafix_adjust_mount_does_not_reduce_numbers();
    test_pafix_mixed_arcs_corrupts_fit();

    test_pafix_tracker_inactive_returns_invalid();
    test_pafix_tracker_idle_returns_drifted_pole();
    test_pafix_tracker_detects_mount_adjustment();

    if (failures) {
        std::printf("\n%d failure(s)\n", failures);
        return 1;
    }
    std::printf("\nAll polar-alignment tests passed.\n");
    return 0;
}
