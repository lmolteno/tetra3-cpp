#pragma once

#include <esp_err.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// WiFi initialization and management
esp_err_t wifi_manager_init(void);
esp_err_t wifi_manager_connect(const char* ssid, const char* password);
void wifi_manager_deinit(void);

// Status and info
bool wifi_manager_is_connected(void);
esp_err_t wifi_manager_get_ip(char* ip_str, size_t max_len);

#ifdef __cplusplus
}
#endif