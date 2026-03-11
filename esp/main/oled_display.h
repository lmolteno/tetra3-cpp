#pragma once

#include <esp_err.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// OLED configuration for ESP32-S3-XIAO
#define OLED_SCL_GPIO 2   // D1 on ESP32-S3-XIAO
#define OLED_SDA_GPIO 1   // D0 on ESP32-S3-XIAO
#define OLED_WIDTH    128
#define OLED_HEIGHT   64
#define OLED_I2C_ADDR 0x3C

/**
 * Initialize the OLED display
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_init(void);

/**
 * Show "Hello world!" message on the display
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_startup_screen(void);

/**
 * Clear the display
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_clear(void);

/**
 * Display a text message
 * @param text The text to display
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_text(const char* text);

/**
 * Update display with current buffer content
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_update(void);

/**
 * Draw WiFi status indicator
 * @param connected True if WiFi is connected, false otherwise
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_wifi_status(bool connected);

/**
 * Draw star centroids on display
 * @param centroids Array of centroid coordinates
 * @param num_centroids Number of centroids
 * @param img_width Image width for scaling
 * @param img_height Image height for scaling
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_centroids(const float* centroids_x, const float* centroids_y, int num_centroids, int img_width, int img_height);

/**
 * Display solve results (RA/DEC)
 * @param solved True if solve was successful
 * @param ra_deg Right Ascension in degrees
 * @param dec_deg Declination in degrees
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_solve_results(bool solved, float ra_deg, float dec_deg);

/**
 * Start the continuous UI loop (capture, detect, solve, display)
 * @return ESP_OK on success, error code on failure
 */
esp_err_t oled_display_start_ui_loop(void);

/**
 * Stop the UI capture loop
 */
void oled_display_stop_capture(void);

/**
 * Start the UI capture loop
 */
void oled_display_start_capture(void);

/**
 * Check if UI capture loop is running
 * @return true if running, false otherwise
 */
bool oled_display_is_capture_running(void);

/**
 * Deinitialize the OLED display
 */
void oled_display_deinit(void);

#ifdef __cplusplus
}
#endif