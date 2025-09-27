#include <stdio.h>
#include <esp_log.h>
#include <esp_system.h>
#include <nvs_flash.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>

#include "camera_ops.h"
#include "http_server.h"
#include "wifi_manager.h"
#include "star_solver.h"
#include "star_catalog.h"
#include "pattern_catalog.h"
// #include "example_catalogs.h"

// Define your WiFi credentials here - replace with your network
// #define WIFI_SSID "Farm"
// #define WIFI_PASS "ragamuffin"
#define WIFI_SSID "Geocam"
#define WIFI_PASS "geocam360wifi"

static const char* TAG = "astro_explorer";

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

    // Initialize WiFi
    ESP_LOGI(TAG, "Initializing WiFi...");
    ESP_ERROR_CHECK(wifi_manager_init());

    // Connect to WiFi - replace with your credentials
    ESP_LOGI(TAG, "Connecting to WiFi...");
    if (wifi_manager_connect(WIFI_SSID, WIFI_PASS) != ESP_OK) {
        ESP_LOGE(TAG, "Failed to connect to WiFi. Please check credentials.");
        return;
    }

    char ip_str[16];
    if (wifi_manager_get_ip(ip_str, sizeof(ip_str)) == ESP_OK) {
        ESP_LOGI(TAG, "Connected! IP address: %s", ip_str);
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

    // Initialize and start HTTP server
    ESP_LOGI(TAG, "Starting HTTP server...");
    ESP_ERROR_CHECK(http_server_init());
    httpd_handle_t server = http_server_start();
    if (!server) {
        ESP_LOGE(TAG, "Failed to start HTTP server");
        star_solver_destroy(solver);
        return;
    }

    ESP_LOGI(TAG, "=== OpenAstroExplorer Ready ===");
    ESP_LOGI(TAG, "Camera endpoints:");
    ESP_LOGI(TAG, "  http://%s/capture - Single frame", ip_str);
    ESP_LOGI(TAG, "  http://%s/accumulate?seconds=10 - Stacked frames", ip_str);
    ESP_LOGI(TAG, "Star processing:");
    ESP_LOGI(TAG, "  http://%s/detect?seconds=5&threshold=28 - Star detection", ip_str);
    ESP_LOGI(TAG, "  POST http://%s/solve - Full star solver", ip_str);
    ESP_LOGI(TAG, "System:");
    ESP_LOGI(TAG, "  http://%s/health - Health check", ip_str);

    // Main loop - keep the application running
    while (1) {
        // Check WiFi connection status
        if (!wifi_manager_is_connected()) {
            ESP_LOGW(TAG, "WiFi connection lost");
        }

        vTaskDelay(pdMS_TO_TICKS(10000)); // Check every 10 seconds
    }

    // Cleanup (won't reach here in normal operation)
    http_server_stop(server);
    star_solver_destroy(solver);
    camera_ops_deinit();
    wifi_manager_deinit();
}
