#include <stdio.h>
#include <esp_log.h>
#include <esp_system.h>
#include <nvs_flash.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <driver/uart.h>

#include "camera_ops.h"
#include "http_server.h"
#include "wifi_manager.h"
#include "star_solver.h"
#include "star_catalog.h"
#include "pattern_catalog.h"
// #include "example_catalogs.h"
#include "oled_display.h"

// Define your WiFi credentials here - replace with your network
// #define WIFI_SSID "Farm"
#define WIFI_SSID "auroraborealis"
#define WIFI_PASS "ragamuffin"
// #define WIFI_SSID "Geocam"
// #define WIFI_PASS "geocam360wifi"

static const char* TAG = "astro_explorer";

// Global control for UI capture loop
static bool ui_capture_running = false;

// UART configuration for keyboard input (works with ESP32-S3 USB CDC)
static void init_uart_for_keyboard(void)
{
    // ESP32-S3 uses UART0 for USB CDC by default
    ESP_LOGI(TAG, "UART0 already configured for USB CDC - ready for keyboard input");
}

static void check_keyboard_input(void)
{
    // Use getchar which is non-blocking when properly configured
    int c = getchar();

    if (c != EOF) {
        char key = (char)c;

        if (key == 's') {
            ESP_LOGI(TAG, "Key 's' pressed - starting capture");
            oled_display_start_capture();
        } else if (key == 'S') {
            ESP_LOGI(TAG, "Key 'S' pressed - stopping capture");
            oled_display_stop_capture();
        }
    }
}

void app_main(void)
{
    ESP_LOGI(TAG, "OpenAstroExplorer ESP32 Starting...");

    // Initialize NVS (required for WiFi)
    esp_err_t ret = nvs_flash_init();
    if (ret == ESP_ERR_NVS_NO_FREE_PAGES || ret == ESP_ERR_NVS_NEW_VERSION_FOUND) {
        ESP_ERROR_CHECK(nvs_flash_erase());
        ret = nvs_flash_init();
    }
    ESP_ERROR_CHECK(ret);

    // Initialize OLED display early for boot messages
    ESP_LOGI(TAG, "Initializing OLED display...");
    bool oled_available = false;
    if (oled_display_init() == ESP_OK) {
        oled_display_startup_screen();
        ESP_LOGI(TAG, "OLED display initialized");
        oled_available = true;
    } else {
        ESP_LOGW(TAG, "OLED display initialization failed - continuing without display");
    }

    // Initialize camera
    ESP_LOGI(TAG, "Initializing camera...");
    if (camera_ops_init() != ESP_OK) {
        ESP_LOGE(TAG, "Failed to initialize camera");
        return;
    }

    // Test star solver initialization
    ESP_LOGI(TAG, "Initializing star solver...");
    ESP_LOGI(TAG, "Loading catalogs:");
    ESP_LOGI(TAG, "  Stars: %d entries (%d KB)", STAR_CATALOG_SIZE,
             (STAR_CATALOG_SIZE * sizeof(star_entry_t)) / 1024);
    ESP_LOGI(TAG, "  Patterns: %d entries (%d KB)", PATTERN_CATALOG_SIZE,
             (PATTERN_CATALOG_SIZE * sizeof(pattern_entry_t)) / 1024);

    star_solver_handle_t solver = star_solver_create(
        star_catalog, STAR_CATALOG_SIZE,
        pattern_catalog, PATTERN_CATALOG_SIZE
    );

    if (!solver) {
        ESP_LOGE(TAG, "Failed to create star solver");
        return;
    }

    ESP_LOGI(TAG, "Star solver created successfully");
    star_solver_print_memory_usage(solver);

    // Start the OLED UI loop immediately - this will run in a separate task
    if (oled_available) {
        ESP_LOGI(TAG, "Starting OLED UI loop...");
        esp_err_t ui_result = oled_display_start_ui_loop();
        if (ui_result != ESP_OK) {
            ESP_LOGE(TAG, "Failed to start OLED UI loop");
        } else {
            ESP_LOGI(TAG, "OLED UI loop started - star detection will begin immediately");
            // Give the UI task a chance to start
            vTaskDelay(pdMS_TO_TICKS(100));
        }
    }

    // Initialize WiFi
    ESP_LOGI(TAG, "Initializing WiFi...");
    ESP_ERROR_CHECK(wifi_manager_init());

    // Connect to WiFi - replace with your credentials
    ESP_LOGI(TAG, "Connecting to WiFi...");
    char ip_str[16] = "No WiFi";
    bool wifi_connected = false;

    if (wifi_manager_connect(WIFI_SSID, WIFI_PASS) == ESP_OK) {
        if (wifi_manager_get_ip(ip_str, sizeof(ip_str)) == ESP_OK) {
            ESP_LOGI(TAG, "Connected! IP address: %s", ip_str);
            wifi_connected = true;
        }
    } else {
        ESP_LOGW(TAG, "Failed to connect to WiFi. Continuing without network connectivity.");
        ESP_LOGW(TAG, "Star detection and solving will still work, but no web interface.");
    }


    // Initialize and start HTTP server only if WiFi is connected
    httpd_handle_t server = NULL;
    if (wifi_connected) {
        ESP_LOGI(TAG, "Starting HTTP server...");
        ESP_ERROR_CHECK(http_server_init());
        server = http_server_start();
        if (!server) {
            ESP_LOGE(TAG, "Failed to start HTTP server");
            // Don't return here - continue without web interface
            wifi_connected = false;
        }
    }

    // Initialize UART for keyboard input
    init_uart_for_keyboard();

    ESP_LOGI(TAG, "=== OpenAstroExplorer Ready ===");
    ESP_LOGI(TAG, "Keyboard controls (via USB):");
    ESP_LOGI(TAG, "  's' - Start capture loop");
    ESP_LOGI(TAG, "  'S' - Stop capture loop");
    if (wifi_connected && server) {
        ESP_LOGI(TAG, "WiFi connected - Web interface available:");
        ESP_LOGI(TAG, "  http://%s/capture - Single frame", ip_str);
        ESP_LOGI(TAG, "  http://%s/accumulate?seconds=10 - Stacked frames", ip_str);
        ESP_LOGI(TAG, "  http://%s/stream - MJPEG stream", ip_str);
        ESP_LOGI(TAG, "  http://%s/detect?seconds=5&threshold=28 - Star detection", ip_str);
        ESP_LOGI(TAG, "  POST http://%s/solve - Full star solver", ip_str);
        ESP_LOGI(TAG, "  http://%s/health - Health check", ip_str);
    } else {
        ESP_LOGI(TAG, "Operating in standalone mode (no WiFi/web interface)");
        ESP_LOGI(TAG, "Star detection and solving available via OLED display");
    }

    // Main loop - keep the application running and monitor system
    while (1) {
        // Check for keyboard input
        check_keyboard_input();

        // Only check WiFi connection if we originally had WiFi
        if (wifi_connected && !wifi_manager_is_connected()) {
            ESP_LOGW(TAG, "WiFi connection lost - attempting reconnection");
            if (wifi_manager_connect(WIFI_SSID, WIFI_PASS) == ESP_OK) {
                ESP_LOGI(TAG, "WiFi reconnected successfully");
            }
        }

        vTaskDelay(pdMS_TO_TICKS(100)); // Check every 100ms for responsive keyboard input
    }

    // Cleanup (won't reach here in normal operation)
    if (server) {
        http_server_stop(server);
    }
    star_solver_destroy(solver);
    camera_ops_deinit();
    wifi_manager_deinit();
    oled_display_deinit();
}
