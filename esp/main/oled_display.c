#include "oled_display.h"
#include "wifi_manager.h"
#include "camera_ops.h"
#include "star_detection.h"
#include "star_solver.h"
#include "example_catalogs.h"
#include <esp_log.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <u8g2.h>
#include <u8g2_esp32_hal.h>

static const char* TAG = "oled_display";

// u8g2 handle and initialization state
static u8g2_t u8g2;
static bool u8g2_initialized = false;

// u8g2 HAL structure for ESP32
static u8g2_esp32_hal_t u8g2_esp32_hal = U8G2_ESP32_HAL_DEFAULT;

// Helper function to clear buffer
static void clear_buffer(void) {
    if (u8g2_initialized) {
        u8g2_ClearBuffer(&u8g2);
    }
}

// Helper function to draw text with u8g2 (better fonts)
static void draw_text_u8g2(int x, int y, const char* text) {
    if (u8g2_initialized) {
        // Use a nice 8-pixel font for small text
        u8g2_SetFont(&u8g2, u8g2_font_6x10_tf);
        u8g2_DrawStr(&u8g2, x, y + 8, text);  // y+8 because u8g2 uses baseline
    }
}

// Helper function to draw text with larger font
static void draw_text_large_u8g2(int x, int y, const char* text) {
    if (u8g2_initialized) {
        // Use a larger font for headers/important text
        u8g2_SetFont(&u8g2, u8g2_font_8x13_tf);
        u8g2_DrawStr(&u8g2, x, y + 10, text);  // y+10 because u8g2 uses baseline
    }
}

// Helper function to update display
static void update_display(void) {
    if (u8g2_initialized) {
        u8g2_SendBuffer(&u8g2);
    }
}

// Helper function to draw a dot for centroids
static void draw_dot(int x, int y) {
    if (u8g2_initialized) {
        u8g2_DrawPixel(&u8g2, x, y);
        // Draw a 3x3 dot for better visibility
        u8g2_DrawPixel(&u8g2, x-1, y);
        u8g2_DrawPixel(&u8g2, x+1, y);
        u8g2_DrawPixel(&u8g2, x, y-1);
        u8g2_DrawPixel(&u8g2, x, y+1);
    }
}

esp_err_t oled_display_init(void)
{
    ESP_LOGI(TAG, "Initializing SH1106 OLED display with u8g2");

    // Initialize u8g2 HAL for ESP32 I2C communication
    ESP_LOGI(TAG, "Configuring u8g2 ESP32 HAL for I2C on SCL=%d, SDA=%d", OLED_SCL_GPIO, OLED_SDA_GPIO);
    u8g2_esp32_hal.bus.i2c.sda = OLED_SDA_GPIO;
    u8g2_esp32_hal.bus.i2c.scl = OLED_SCL_GPIO;
    u8g2_esp32_hal_init(u8g2_esp32_hal);

    // Initialize u8g2 for SH1106 display
    ESP_LOGI(TAG, "Setting up u8g2 for SH1106 128x64 display");
    u8g2_Setup_sh1106_i2c_128x64_noname_f(&u8g2, U8G2_R2, u8g2_esp32_i2c_byte_cb, u8g2_esp32_gpio_and_delay_cb);
    u8x8_SetI2CAddress(&u8g2.u8x8, OLED_I2C_ADDR * 2);  // u8g2 expects 8-bit address (0x3C -> 0x78)

    // Initialize display
    ESP_LOGI(TAG, "Initializing display hardware...");
    u8g2_InitDisplay(&u8g2);
    u8g2_SetPowerSave(&u8g2, 0);
    u8g2_ClearDisplay(&u8g2);

    u8g2_initialized = true;
    ESP_LOGI(TAG, "SH1106 OLED display initialized successfully with u8g2");
    return ESP_OK;
}

esp_err_t oled_display_startup_screen(void)
{
    if (!u8g2_initialized) {
        ESP_LOGE(TAG, "OLED not initialized");
        return ESP_ERR_INVALID_STATE;
    }

    ESP_LOGI(TAG, "Displaying startup screen on OLED");

    // Clear the screen first
    clear_buffer();

    // Display startup information with large fonts
    draw_text_large_u8g2(8, 0, "Plate Solving");
    draw_text_large_u8g2(32, 16, "Camera");
    draw_text_u8g2(16, 36, "Linus Molteno");
    draw_text_u8g2(48, 52, "2025-09");

    // Update display
    update_display();

    ESP_LOGI(TAG, "Startup screen displayed successfully");
    return ESP_OK;
}

esp_err_t oled_display_clear(void)
{
    if (!u8g2_initialized) {
        ESP_LOGE(TAG, "OLED not initialized");
        return ESP_ERR_INVALID_STATE;
    }

    clear_buffer();
    update_display();
    return ESP_OK;
}

esp_err_t oled_display_text(const char* text)
{
    if (!u8g2_initialized) {
        ESP_LOGE(TAG, "OLED not initialized");
        return ESP_ERR_INVALID_STATE;
    }

    if (text == NULL) {
        ESP_LOGE(TAG, "Text parameter is NULL");
        return ESP_ERR_INVALID_ARG;
    }

    // Clear screen first
    clear_buffer();

    // Display the text with better fonts
    draw_text_u8g2(2, 16, text);

    // Update display
    update_display();
    return ESP_OK;
}

esp_err_t oled_display_update(void)
{
    if (!u8g2_initialized) {
        ESP_LOGE(TAG, "OLED not initialized");
        return ESP_ERR_INVALID_STATE;
    }

    update_display();
    return ESP_OK;
}

void oled_display_deinit(void)
{
    ESP_LOGI(TAG, "Deinitializing OLED display");

    // Clear display
    if (u8g2_initialized) {
        u8g2_ClearDisplay(&u8g2);
        u8g2_initialized = false;
    }

    ESP_LOGI(TAG, "OLED display deinitialized");
}

esp_err_t oled_display_wifi_status(bool connected)
{
    if (!u8g2_initialized) {
        return ESP_ERR_INVALID_STATE;
    }

    // Draw WiFi indicator in top-right corner with better font
    const char* wifi_symbol = connected ? "W" : "X";
    draw_text_u8g2(118, 0, wifi_symbol);
    return ESP_OK;
}

esp_err_t oled_display_centroids(const float* centroids_x, const float* centroids_y, int num_centroids, int img_width, int img_height)
{
    if (!u8g2_initialized) {
        return ESP_ERR_INVALID_STATE;
    }

    // Scale centroids to display size
    // Use area from y=8 to y=47 (40 pixels height) for centroids
    const int display_area_width = 118;  // Leave space for WiFi indicator and account for 2px offset
    const int display_area_height = 40;
    const int display_area_y_start = 8;
    const int x_offset = 2;  // Account for 2-pixel horizontal shift

    // Draw up to 20 centroids to avoid overcrowding
    int max_centroids = (num_centroids > 20) ? 20 : num_centroids;

    for (int i = 0; i < max_centroids; i++) {
        // Scale centroid position to display area
        int x = x_offset + (int)((centroids_x[i] / img_width) * display_area_width);
        int y = display_area_y_start + (int)((centroids_y[i] / img_height) * display_area_height);

        // Ensure coordinates are within bounds
        if (x >= x_offset && x < (display_area_width + x_offset) && y >= display_area_y_start && y < (display_area_y_start + display_area_height)) {
            // Draw a dot for the centroid with better visibility
            draw_dot(x, y);
        }
    }

    return ESP_OK;
}

esp_err_t oled_display_solve_results(bool solved, float ra_deg, float dec_deg)
{
    if (!u8g2_initialized) {
        return ESP_ERR_INVALID_STATE;
    }

    char result_text[32];

    if (solved) {
        // Convert RA to hours and minutes
        float ra_hours = ra_deg / 15.0f;  // 15 degrees per hour
        int ra_h = (int)ra_hours;
        int ra_m = (int)((ra_hours - ra_h) * 60.0f);

        // Convert DEC to degrees and minutes
        int dec_d = (int)dec_deg;
        int dec_m = (int)((fabs(dec_deg) - abs(dec_d)) * 60.0f);

        snprintf(result_text, sizeof(result_text), "RA:%02dh%02dm DEC:%+03d%02d",
                 ra_h % 24, ra_m, dec_d, dec_m);
    } else {
        snprintf(result_text, sizeof(result_text), "No solution");
    }

    // Display results in bottom area (y=48 to y=63) with better fonts
    draw_text_u8g2(2, 48, result_text);
    return ESP_OK;
}

// Static task handle for the UI loop
static TaskHandle_t ui_task_handle = NULL;

// Global capture control flag
static volatile bool capture_enabled = false; // Start running by default, user can control with 's'/'S'

// Persistent display state to maintain across cycles
static struct {
    float* last_centroids_x;
    float* last_centroids_y;
    int last_num_centroids;
    int last_img_width;
    int last_img_height;
    bool last_solved;
    float last_ra_deg;
    float last_dec_deg;
    char last_star_count_text[17];
} display_state = {
    .last_centroids_x = NULL,
    .last_centroids_y = NULL,
    .last_num_centroids = 0,
    .last_img_width = 0,
    .last_img_height = 0,
    .last_solved = false,
    .last_ra_deg = 0.0f,
    .last_dec_deg = 0.0f,
    .last_star_count_text = "No data"
};

// Helper function to update the persistent display state
static void update_display_state(float* centroids_x, float* centroids_y, int num_centroids,
                                int img_width, int img_height, bool solved,
                                float ra_deg, float dec_deg) {
    // Free old centroid data
    if (display_state.last_centroids_x) {
        free(display_state.last_centroids_x);
        display_state.last_centroids_x = NULL;
    }
    if (display_state.last_centroids_y) {
        free(display_state.last_centroids_y);
        display_state.last_centroids_y = NULL;
    }

    // Store new centroid data
    if (centroids_x && centroids_y && num_centroids > 0) {
        display_state.last_centroids_x = malloc(num_centroids * sizeof(float));
        display_state.last_centroids_y = malloc(num_centroids * sizeof(float));
        if (display_state.last_centroids_x && display_state.last_centroids_y) {
            memcpy(display_state.last_centroids_x, centroids_x, num_centroids * sizeof(float));
            memcpy(display_state.last_centroids_y, centroids_y, num_centroids * sizeof(float));
            display_state.last_num_centroids = num_centroids;
            display_state.last_img_width = img_width;
            display_state.last_img_height = img_height;
            snprintf(display_state.last_star_count_text, sizeof(display_state.last_star_count_text), "%d stars", num_centroids);
        }
    }

    // Store solve results
    display_state.last_solved = solved;
    display_state.last_ra_deg = ra_deg;
    display_state.last_dec_deg = dec_deg;
}

// Helper function to show the current persistent display state
static void show_persistent_display(bool wifi_connected, const char* status_text) {
    clear_buffer();
    oled_display_wifi_status(wifi_connected);

    // Show current status
    draw_text_u8g2(2, 0, status_text);

    // Show last known centroids if available
    if (display_state.last_centroids_x && display_state.last_centroids_y && display_state.last_num_centroids > 0) {
        oled_display_centroids(display_state.last_centroids_x, display_state.last_centroids_y,
                              display_state.last_num_centroids, display_state.last_img_width,
                              display_state.last_img_height);
    }

    // Show last solve results
    oled_display_solve_results(display_state.last_solved, display_state.last_ra_deg, display_state.last_dec_deg);

    update_display();
}

static void ui_loop_task(void *pvParameters)
{
    ESP_LOGI(TAG, "UI loop task started");

    // Create star solver once
    star_solver_handle_t solver = star_solver_create(
        star_catalog, STAR_CATALOG_SIZE,
        pattern_catalog, PATTERN_CATALOG_SIZE
    );

    if (!solver) {
        ESP_LOGE(TAG, "Failed to create star solver for UI loop");
        vTaskDelete(NULL);
        return;
    }

    while (1) {
        // Check if capture is enabled
        if (!capture_enabled) {
            // Show WiFi status with idle display
            bool wifi_connected = wifi_manager_is_connected();
            show_persistent_display(wifi_connected, "Idle - press 's'");
            vTaskDelay(pdMS_TO_TICKS(500));
            continue;
        }

        // Show WiFi status
        bool wifi_connected = wifi_manager_is_connected();

        // Show current status with persistent data
        show_persistent_display(wifi_connected, "Capturing...");

        // Log memory status before capture
        size_t free_heap = esp_get_free_heap_size();
        size_t free_psram = heap_caps_get_free_size(MALLOC_CAP_SPIRAM);
        ESP_LOGI(TAG, "Memory before capture: Heap=%zu bytes, PSRAM=%zu bytes", free_heap, free_psram);

        // Yield before long operation
        vTaskDelay(pdMS_TO_TICKS(100));

        // Capture frame (2 second accumulation)
        camera_fb_t* fb = camera_ops_capture_accumulated(2);
        if (!fb) {
            ESP_LOGE(TAG, "Failed to capture frame - checking memory");
            free_heap = esp_get_free_heap_size();
            free_psram = heap_caps_get_free_size(MALLOC_CAP_SPIRAM);
            ESP_LOGE(TAG, "Memory after failed capture: Heap=%zu bytes, PSRAM=%zu bytes", free_heap, free_psram);

            show_persistent_display(wifi_connected, "Capture FAIL");
            vTaskDelay(pdMS_TO_TICKS(2000)); // Wait longer before retry
            continue;
        }

        // Update status
        show_persistent_display(wifi_connected, "Detecting...");

        // Yield before detection
        vTaskDelay(pdMS_TO_TICKS(50));

#ifdef DOWNSAMPLE

        // Downsample image by 2x2 averaging for faster star detection
        int downsampled_width = fb->width / 2;
        int downsampled_height = fb->height / 2;
        size_t downsampled_size = downsampled_width * downsampled_height;

        uint8_t* downsampled_data = malloc(downsampled_size);
        if (!downsampled_data) {
            ESP_LOGE(TAG, "Failed to allocate memory for downsampled image");
            esp_camera_fb_return(fb);
            continue;
        }

        // Perform 2x2 averaging downsampling
        for (int y = 0; y < downsampled_height; y++) {
            for (int x = 0; x < downsampled_width; x++) {
                // Get 2x2 block from original image
                int orig_x = x * 2;
                int orig_y = y * 2;

                uint32_t sum = 0;
                sum += fb->buf[orig_y * fb->width + orig_x];                    // top-left
                sum += fb->buf[orig_y * fb->width + orig_x + 1];                // top-right
                sum += fb->buf[(orig_y + 1) * fb->width + orig_x];              // bottom-left
                sum += fb->buf[(orig_y + 1) * fb->width + orig_x + 1];          // bottom-right

                // Average the 4 pixels
                downsampled_data[y * downsampled_width + x] = (uint8_t)(sum / 4);
            }
        }

        // Create a temporary framebuffer structure for the downsampled data
        camera_fb_t downsampled_fb = {
            .buf = downsampled_data,
            .len = downsampled_size,
            .width = downsampled_width,
            .height = downsampled_height,
            .format = fb->format,
            .timestamp = fb->timestamp
        };

        ESP_LOGI(TAG, "Downsampled image from %dx%d to %dx%d", fb->width, fb->height, downsampled_width, downsampled_height);

        // Detect stars on downsampled image
        star_detection_result_t* detection = detect_stars_simple(&downsampled_fb, 60.0f, 2, 100);
#else
        // Detect stars on original
        star_detection_result_t* detection = detect_stars_simple(fb, 60.0f, 2, 100);
#endif

        // Scale centroids back to original image coordinates
        if (detection && detection->num_centroids > 0) {
            for (int i = 0; i < detection->num_centroids; i++) {
                detection->centroids[i].x *= 2.0f;  // Scale back to original width
                detection->centroids[i].y *= 2.0f;  // Scale back to original height
            }
        }

#ifdef DOWNSAMPLE
        // Free downsampled data
        free(downsampled_data);
#endif

        if (!detection || detection->num_centroids == 0) {
            ESP_LOGW(TAG, "No stars detected");
            esp_camera_fb_return(fb);

            // Update display state to show no stars
            update_display_state(NULL, NULL, 0, 0, 0, false, 0, 0);
            show_persistent_display(wifi_connected, "No stars");

            if (detection) star_detection_free(detection);
            vTaskDelay(pdMS_TO_TICKS(2000));
            continue;
        }

        ESP_LOGI(TAG, "Detected %d stars", detection->num_centroids);

        // Prepare centroids for display
        float* centroids_x = malloc(detection->num_centroids * sizeof(float));
        float* centroids_y = malloc(detection->num_centroids * sizeof(float));

        if (centroids_x && centroids_y) {
            for (int i = 0; i < detection->num_centroids; i++) {
                centroids_x[i] = detection->centroids[i].x;
                centroids_y[i] = detection->centroids[i].y;
            }

            // Update display state with new centroids (but keep old solve results for now)
            update_display_state(centroids_x, centroids_y, detection->num_centroids, fb->width, fb->height,
                                display_state.last_solved, display_state.last_ra_deg, display_state.last_dec_deg);
        }

        // Show solving status with updated centroids
        show_persistent_display(wifi_connected, "Solving...");

        // Yield before solving
        vTaskDelay(pdMS_TO_TICKS(50));

        // Attempt to solve
        solve_result_t solve_result = star_solver_solve_from_centroids(
            solver,
            detection->centroids, detection->num_centroids,
            fb->height, fb->width,
            30.0,      // fov_estimate_deg
            16,       // pattern_checking_stars
            0.01f,    // match_radius
            0.001f,   // match_threshold
            false,    // use_distortion
            0.0f      // distortion_coeff
        );

        // Process solve results and update display state
        if (solve_result.solved) {
            ESP_LOGI(TAG, "SOLVED! RA: %.2f°, Dec: %.2f°",
                     solve_result.ra * 180.0f / M_PI,
                     solve_result.dec * 180.0f / M_PI);

            // Update display state with solve results
            if (centroids_x && centroids_y) {
                update_display_state(centroids_x, centroids_y, detection->num_centroids, fb->width, fb->height,
                                    true, solve_result.ra * 180.0f / M_PI, solve_result.dec * 180.0f / M_PI);
            }
        } else {
            ESP_LOGI(TAG, "No solution found");
            // Update display state without solve results
            if (centroids_x && centroids_y) {
                update_display_state(centroids_x, centroids_y, detection->num_centroids, fb->width, fb->height,
                                    false, 0, 0);
            }
        }

        // Show final results
        char star_text[16];
        snprintf(star_text, sizeof(star_text), "%d stars", detection->num_centroids);
        show_persistent_display(wifi_connected, star_text);

        // Clean up temporary centroid arrays
        if (centroids_x) free(centroids_x);
        if (centroids_y) free(centroids_y);

        // Cleanup
        esp_camera_fb_return(fb);
        star_detection_free(detection);

        // Log memory status after cleanup
        free_heap = esp_get_free_heap_size();
        free_psram = heap_caps_get_free_size(MALLOC_CAP_SPIRAM);
        ESP_LOGI(TAG, "Memory after cleanup: Heap=%zu bytes, PSRAM=%zu bytes", free_heap, free_psram);

        // Wait before next iteration - yield to other tasks
        vTaskDelay(pdMS_TO_TICKS(1000)); // Longer delay to allow memory recovery

        // Additional yield points to prevent watchdog
        taskYIELD();
    }

    // This should never be reached, but cleanup if it is
    star_solver_destroy(solver);

    // Clean up persistent display state
    if (display_state.last_centroids_x) {
        free(display_state.last_centroids_x);
        display_state.last_centroids_x = NULL;
    }
    if (display_state.last_centroids_y) {
        free(display_state.last_centroids_y);
        display_state.last_centroids_y = NULL;
    }

    vTaskDelete(NULL);
}

esp_err_t oled_display_start_ui_loop(void)
{
    if (!u8g2_initialized) {
        ESP_LOGE(TAG, "OLED not initialized");
        return ESP_ERR_INVALID_STATE;
    }

    ESP_LOGI(TAG, "Creating UI loop task");

    // Create the UI loop task with high priority to start immediately
    BaseType_t result = xTaskCreate(
        ui_loop_task,           // Function that implements the task
        "oled_ui_loop",         // Text name for the task
        8192,                   // Stack size in bytes
        NULL,                   // Parameter passed into the task
        10,                     // Higher priority to start before WiFi blocks
        &ui_task_handle         // Task handle
    );

    if (result != pdPASS) {
        ESP_LOGE(TAG, "Failed to create UI loop task");
        return ESP_FAIL;
    }

    ESP_LOGI(TAG, "UI loop task created successfully");
    return ESP_OK;
}

void oled_display_stop_capture(void)
{
    capture_enabled = false;
    ESP_LOGI(TAG, "UI capture loop stopped");
}

void oled_display_start_capture(void)
{
    capture_enabled = true;
    ESP_LOGI(TAG, "UI capture loop started");
}

bool oled_display_is_capture_running(void)
{
    return capture_enabled;
}