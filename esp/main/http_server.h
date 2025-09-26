#pragma once

#include <esp_err.h>
#include <esp_http_server.h>

#ifdef __cplusplus
extern "C" {
#endif

// HTTP server lifecycle
esp_err_t http_server_init(void);
void http_server_deinit(void);
httpd_handle_t http_server_start(void);
void http_server_stop(httpd_handle_t server);

#ifdef __cplusplus
}
#endif