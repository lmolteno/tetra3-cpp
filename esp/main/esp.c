#include <stdio.h>
#include "star_solver.h"
#include "star_catalog.h"
#include "pattern_catalog.h"

void app_main(void)
{
    printf("ESP32 Star Solver Test\n");

    // Create some test centroids (example star positions from real solve)
    centroid_t test_centroids[] = {
        {692.4872, 491.22507},
        {946.2311, 332.36243},
        {12.089273, 442.55933},
        {943.6792, 298.08932},
        {734.17804, 297.76077},
        {1335.9685, 518.3289},
        {1374.8608, 493.86456},
        {905.30676, 585.91943}
    };
    size_t num_centroids = sizeof(test_centroids) / sizeof(test_centroids[0]);

    printf("Loading catalogs:\n");
    printf("  Stars: %d entries (%d KB)\n", STAR_CATALOG_SIZE, (STAR_CATALOG_SIZE * sizeof(star_entry_t)) / 1024);
    printf("  Patterns: %d entries (%d KB)\n", PATTERN_CATALOG_SIZE, (PATTERN_CATALOG_SIZE * sizeof(pattern_entry_t)) / 1024);

    // Create star solver with real catalog data
    star_solver_handle_t solver = star_solver_create(
        star_catalog, STAR_CATALOG_SIZE,
        pattern_catalog, PATTERN_CATALOG_SIZE
    );

    if (!solver) {
        printf("Failed to create star solver\n");
        return;
    }

    printf("Star solver created successfully\n");
    star_solver_print_memory_usage(solver);

    // Test solve with centroids (using parameters from example)
    solve_result_t result = star_solver_solve_from_centroids(
        solver,
        test_centroids, num_centroids,
        1200, 1600, // height, width (from example)
        7.0,        // fov_estimate_deg (from example)
        8,          // pattern_checking_stars
        0.01f,      // match_radius
        0.001f,     // match_threshold
        false,      // use_distortion
        0.0f        // distortion_coeff
    );

    printf("Solve completed in %.2f ms\n", result.solve_time_ms);
    if (result.solved) {
        printf("SOLVED! RA: %.2f deg, Dec: %.2f deg, Roll: %.2f deg\n",
               result.ra * 180.0f / 3.14159f,
               result.dec * 180.0f / 3.14159f,
               result.roll * 180.0f / 3.14159f);
        printf("FOV: %.2f deg, RMSE: %.2f arcsec, Matches: %d\n",
               result.fov * 180.0f / 3.14159f, result.rmse, result.num_matches);
    } else {
        printf("No solution found\n");
    }

    // Clean up
    star_solver_destroy(solver);
    printf("Test completed\n");
}
