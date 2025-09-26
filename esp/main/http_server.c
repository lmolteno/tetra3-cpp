#include "http_server.h"
#include "camera_ops.h"
#include "star_solver.h"
#include "star_detection.h"
#include "star_catalog.h"
#include "pattern_catalog.h"
#include <esp_log.h>
#include <esp_timer.h>
#include <string.h>
#include <stdlib.h>

static const char* TAG = "http_server";

// Forward declarations of handler functions
static esp_err_t capture_handler(httpd_req_t *req);
static esp_err_t accumulate_handler(httpd_req_t *req);
static esp_err_t detect_handler(httpd_req_t *req);
static esp_err_t solve_handler(httpd_req_t *req);
static esp_err_t health_handler(httpd_req_t *req);

// Utility functions
static size_t jpg_encode_stream(void *arg, size_t index, const void *data, size_t len);

typedef struct {
    httpd_req_t *req;
    size_t len;
} jpg_chunking_t;

esp_err_t http_server_init(void)
{
    ESP_LOGI(TAG, "HTTP server module initialized");
    return ESP_OK;
}

void http_server_deinit(void)
{
    ESP_LOGI(TAG, "HTTP server module deinitialized");
}

httpd_handle_t http_server_start(void)
{
    httpd_handle_t server = NULL;
    httpd_config_t config = HTTPD_DEFAULT_CONFIG();
    config.max_uri_handlers = 8;
    config.server_port = 80;

    ESP_LOGI(TAG, "Starting HTTP server on port %d", config.server_port);

    if (httpd_start(&server, &config) != ESP_OK) {
        ESP_LOGE(TAG, "Error starting HTTP server");
        return NULL;
    }

    // Register URI handlers
    httpd_uri_t health_uri = {
        .uri = "/health",
        .method = HTTP_GET,
        .handler = health_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &health_uri);

    httpd_uri_t capture_uri = {
        .uri = "/capture",
        .method = HTTP_GET,
        .handler = capture_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &capture_uri);

    httpd_uri_t accumulate_uri = {
        .uri = "/accumulate",
        .method = HTTP_GET,
        .handler = accumulate_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &accumulate_uri);

    httpd_uri_t detect_uri = {
        .uri = "/detect",
        .method = HTTP_GET,
        .handler = detect_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &detect_uri);

    httpd_uri_t solve_uri = {
        .uri = "/solve",
        .method = HTTP_POST,
        .handler = solve_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &solve_uri);

    ESP_LOGI(TAG, "Available endpoints:");
    ESP_LOGI(TAG, "  GET /health - Server health check");
    ESP_LOGI(TAG, "  GET /capture - Single frame JPEG");
    ESP_LOGI(TAG, "  GET /accumulate?seconds=N - Accumulated JPEG (1-60s)");
    ESP_LOGI(TAG, "  GET /detect?seconds=N&threshold=X - Star detection only");
    ESP_LOGI(TAG, "  POST /solve - Star solver endpoint");

    return server;
}

void http_server_stop(httpd_handle_t server)
{
    if (server) {
        httpd_stop(server);
        ESP_LOGI(TAG, "HTTP server stopped");
    }
}

static esp_err_t health_handler(httpd_req_t *req)
{
    const char* response = "{\"status\":\"ok\",\"service\":\"astro-solver\"}";

    httpd_resp_set_type(req, "application/json");
    httpd_resp_send(req, response, strlen(response));

    return ESP_OK;
}

static esp_err_t capture_handler(httpd_req_t *req)
{
    ESP_LOGI(TAG, "Single frame capture request");

    int64_t start_time = esp_timer_get_time();

    camera_fb_t* fb = camera_ops_capture_frame();
    if (!fb) {
        ESP_LOGE(TAG, "Camera capture failed");
        httpd_resp_send_500(req);
        return ESP_FAIL;
    }

    esp_err_t res = httpd_resp_set_type(req, "image/jpeg");
    if (res == ESP_OK) {
        res = httpd_resp_set_hdr(req, "Content-Disposition", "inline; filename=capture.jpg");
    }

    if (res == ESP_OK) {
        if (fb->format == PIXFORMAT_JPEG) {
            res = httpd_resp_send(req, (const char*)fb->buf, fb->len);
        } else {
            // Convert grayscale to JPEG
            jpg_chunking_t jchunk = {req, 0};
            res = frame2jpg_cb(fb, 80, jpg_encode_stream, &jchunk) ? ESP_OK : ESP_FAIL;
            httpd_resp_send_chunk(req, NULL, 0);
        }
    }

    esp_camera_fb_return(fb);

    int64_t end_time = esp_timer_get_time();
    ESP_LOGI(TAG, "Capture completed in %lld ms", (end_time - start_time) / 1000);

    return res;
}

static esp_err_t accumulate_handler(httpd_req_t *req)
{
    ESP_LOGI(TAG, "Accumulated capture request");

    char query[64];
    int seconds = 5; // default

    if (httpd_req_get_url_query_str(req, query, sizeof(query)) == ESP_OK) {
        char param[16];
        if (httpd_query_key_value(query, "seconds", param, sizeof(param)) == ESP_OK) {
            seconds = atoi(param);
            if (seconds < 1 || seconds > 60) seconds = 5;
        }
    }

    ESP_LOGI(TAG, "Starting %d-second accumulated capture", seconds);

    int64_t start_time = esp_timer_get_time();

    camera_fb_t* fb = camera_ops_capture_accumulated(seconds);
    if (!fb) {
        ESP_LOGE(TAG, "Accumulated capture failed");
        httpd_resp_send_500(req);
        return ESP_FAIL;
    }

    esp_err_t res = httpd_resp_set_type(req, "image/jpeg");
    if (res == ESP_OK) {
        res = httpd_resp_set_hdr(req, "Content-Disposition", "inline; filename=accumulated.jpg");
    }

    if (res == ESP_OK) {
        jpg_chunking_t jchunk = {req, 0};
        res = frame2jpg_cb(fb, 80, jpg_encode_stream, &jchunk) ? ESP_OK : ESP_FAIL;
        httpd_resp_send_chunk(req, NULL, 0);
    }

    esp_camera_fb_return(fb);

    int64_t end_time = esp_timer_get_time();
    ESP_LOGI(TAG, "Accumulated capture completed in %lld ms", (end_time - start_time) / 1000);

    return res;
}

static esp_err_t detect_handler(httpd_req_t *req)
{
    ESP_LOGI(TAG, "Star detection request");

    char query[128];
    int seconds = 5;
    float threshold = 28.0f;
    int min_area = 3;

    if (httpd_req_get_url_query_str(req, query, sizeof(query)) == ESP_OK) {
        char param[16];
        if (httpd_query_key_value(query, "seconds", param, sizeof(param)) == ESP_OK) {
            seconds = atoi(param);
            if (seconds < 1 || seconds > 60) seconds = 5;
        }
        if (httpd_query_key_value(query, "threshold", param, sizeof(param)) == ESP_OK) {
            threshold = atof(param);
            if (threshold < 1.0f || threshold > 100.0f) threshold = 28.0f;
        }
        if (httpd_query_key_value(query, "min_area", param, sizeof(param)) == ESP_OK) {
            min_area = atoi(param);
            if (min_area < 1 || min_area > 50) min_area = 3;
        }
    }

    ESP_LOGI(TAG, "Starting detection: %ds exposure, threshold=%.1f, min_area=%d",
             seconds, threshold, min_area);

    int64_t start_time = esp_timer_get_time();

    // Capture accumulated frame
    camera_fb_t* fb = camera_ops_capture_accumulated(seconds);
    if (!fb) {
        ESP_LOGE(TAG, "Failed to capture frame for detection");
        httpd_resp_send_500(req);
        return ESP_FAIL;
    }

    // Detect stars using our algorithm
    star_detection_result_t* detection = detect_stars_simple(fb, threshold, min_area, 100);

    char response[1024];
    if (!detection || detection->num_centroids == 0) {
        ESP_LOGW(TAG, "No stars detected");

        snprintf(response, sizeof(response),
            "{"
            "\"status\":\"no_stars\","
            "\"width\":%d,"
            "\"height\":%d,"
            "\"exposure_seconds\":%d,"
            "\"threshold\":%.1f,"
            "\"min_area\":%d,"
            "\"stars_detected\":0,"
            "\"centroids\":[]"
            "}",
            fb->width, fb->height, seconds, threshold, min_area);
    } else {
        ESP_LOGI(TAG, "Detected %d stars", detection->num_centroids);

        // Build centroid JSON array
        char* centroid_json = malloc(detection->num_centroids * 64); // ~64 chars per centroid
        if (centroid_json) {
            strcpy(centroid_json, "[");
            for (int i = 0; i < detection->num_centroids; i++) {
                char centroid_str[64];
                snprintf(centroid_str, sizeof(centroid_str),
                    "%s{\"x\":%.2f,\"y\":%.2f,\"size\":%.0f}",
                    (i > 0) ? "," : "",
                    detection->centroids[i].x,
                    detection->centroids[i].y,
                    detection->sizes[i]);
                strcat(centroid_json, centroid_str);
            }
            strcat(centroid_json, "]");

            snprintf(response, sizeof(response),
                "{"
                "\"status\":\"detected\","
                "\"width\":%d,"
                "\"height\":%d,"
                "\"exposure_seconds\":%d,"
                "\"threshold\":%.1f,"
                "\"min_area\":%d,"
                "\"stars_detected\":%d,"
                "\"centroids\":%s"
                "}",
                fb->width, fb->height, seconds, threshold, min_area,
                detection->num_centroids, centroid_json);

            free(centroid_json);
        } else {
            // Fallback without detailed centroids
            snprintf(response, sizeof(response),
                "{"
                "\"status\":\"detected\","
                "\"width\":%d,"
                "\"height\":%d,"
                "\"exposure_seconds\":%d,"
                "\"threshold\":%.1f,"
                "\"min_area\":%d,"
                "\"stars_detected\":%d"
                "}",
                fb->width, fb->height, seconds, threshold, min_area,
                detection->num_centroids);
        }
    }

    // Cleanup
    esp_camera_fb_return(fb);
    if (detection) star_detection_free(detection);

    esp_err_t res = httpd_resp_set_type(req, "application/json");
    if (res == ESP_OK) {
        res = httpd_resp_send(req, response, strlen(response));
    }

    int64_t end_time = esp_timer_get_time();
    ESP_LOGI(TAG, "Detection completed in %lld ms", (end_time - start_time) / 1000);

    return res;
}

static esp_err_t solve_handler(httpd_req_t *req)
{
    ESP_LOGI(TAG, "Star solver request");

    char query[128];
    int seconds = 5; // default accumulation time
    float threshold = 28.0f;
    int min_area = 3;

    if (httpd_req_get_url_query_str(req, query, sizeof(query)) == ESP_OK) {
        char param[16];
        if (httpd_query_key_value(query, "seconds", param, sizeof(param)) == ESP_OK) {
            seconds = atoi(param);
            if (seconds < 1 || seconds > 60) seconds = 5;
        }
        if (httpd_query_key_value(query, "threshold", param, sizeof(param)) == ESP_OK) {
            threshold = atof(param);
            if (threshold < 1.0f || threshold > 100.0f) threshold = 28.0f;
        }
        if (httpd_query_key_value(query, "min_area", param, sizeof(param)) == ESP_OK) {
            min_area = atoi(param);
            if (min_area < 1 || min_area > 50) min_area = 3;
        }
    }

    ESP_LOGI(TAG, "Starting solve: %ds exposure, threshold=%.1f, min_area=%d",
             seconds, threshold, min_area);

    int64_t start_time = esp_timer_get_time();

    // Capture accumulated frame
    camera_fb_t* fb = camera_ops_capture_accumulated(seconds);
    if (!fb) {
        ESP_LOGE(TAG, "Failed to capture frame for star solving");
        httpd_resp_send_500(req);
        return ESP_FAIL;
    }

    // Detect stars using our algorithm
    star_detection_result_t* detection = detect_stars_simple(fb, threshold, min_area, 100);
    if (!detection || detection->num_centroids == 0) {
        ESP_LOGW(TAG, "No stars detected in image");

        char response[256];
        snprintf(response, sizeof(response),
            "{"
            "\"status\":\"no_stars\","
            "\"width\":%d,"
            "\"height\":%d,"
            "\"exposure_seconds\":%d,"
            "\"stars_detected\":0"
            "}",
            fb->width, fb->height, seconds);

        esp_camera_fb_return(fb);
        if (detection) star_detection_free(detection);

        esp_err_t res = httpd_resp_set_type(req, "application/json");
        if (res == ESP_OK) {
            res = httpd_resp_send(req, response, strlen(response));
        }
        return res;
    }

    ESP_LOGI(TAG, "Detected %d stars, attempting to solve", detection->num_centroids);

    // Create star solver
    star_solver_handle_t solver = star_solver_create(
        star_catalog, STAR_CATALOG_SIZE,
        pattern_catalog, PATTERN_CATALOG_SIZE
    );

    char response[1024];
    if (!solver) {
        ESP_LOGE(TAG, "Failed to create star solver");

        snprintf(response, sizeof(response),
            "{"
            "\"status\":\"solver_error\","
            "\"width\":%d,"
            "\"height\":%d,"
            "\"exposure_seconds\":%d,"
            "\"stars_detected\":%d"
            "}",
            fb->width, fb->height, seconds, detection->num_centroids);
    } else {
        // Attempt to solve
        solve_result_t solve_result = star_solver_solve_from_centroids(
            solver,
            detection->centroids, detection->num_centroids,
            fb->height, fb->width,    // height, width
            20.0,                     // fov_estimate_deg
            8,                        // pattern_checking_stars
            0.01f,                    // match_radius
            0.001f,                   // match_threshold
            false,                    // use_distortion
            0.0f                      // distortion_coeff
        );

        if (solve_result.solved) {
            ESP_LOGI(TAG, "SOLVED! RA: %.2f°, Dec: %.2f°, Roll: %.2f°",
                     solve_result.ra * 180.0f / 3.14159f,
                     solve_result.dec * 180.0f / 3.14159f,
                     solve_result.roll * 180.0f / 3.14159f);

            snprintf(response, sizeof(response),
                "{"
                "\"status\":\"solved\","
                "\"width\":%d,"
                "\"height\":%d,"
                "\"exposure_seconds\":%d,"
                "\"stars_detected\":%d,"
                "\"ra_deg\":%.6f,"
                "\"dec_deg\":%.6f,"
                "\"roll_deg\":%.6f,"
                "\"fov_deg\":%.6f,"
                "\"rmse_arcsec\":%.2f,"
                "\"num_matches\":%d,"
                "\"solve_time_ms\":%.2f"
                "}",
                fb->width, fb->height, seconds, detection->num_centroids,
                solve_result.ra * 180.0f / 3.14159f,
                solve_result.dec * 180.0f / 3.14159f,
                solve_result.roll * 180.0f / 3.14159f,
                solve_result.fov * 180.0f / 3.14159f,
                solve_result.rmse,
                solve_result.num_matches,
                solve_result.solve_time_ms);
        } else {
            ESP_LOGW(TAG, "No solution found");

            snprintf(response, sizeof(response),
                "{"
                "\"status\":\"no_solution\","
                "\"width\":%d,"
                "\"height\":%d,"
                "\"exposure_seconds\":%d,"
                "\"stars_detected\":%d,"
                "\"solve_time_ms\":%.2f"
                "}",
                fb->width, fb->height, seconds, detection->num_centroids,
                solve_result.solve_time_ms);
        }

        star_solver_destroy(solver);
    }

    // Cleanup
    esp_camera_fb_return(fb);
    star_detection_free(detection);

    esp_err_t res = httpd_resp_set_type(req, "application/json");
    if (res == ESP_OK) {
        res = httpd_resp_send(req, response, strlen(response));
    }

    int64_t end_time = esp_timer_get_time();
    ESP_LOGI(TAG, "Solve request completed in %lld ms", (end_time - start_time) / 1000);

    return res;
}

static size_t jpg_encode_stream(void *arg, size_t index, const void *data, size_t len)
{
    jpg_chunking_t *j = (jpg_chunking_t *)arg;
    if (!index) {
        j->len = 0;
    }
    if (httpd_resp_send_chunk(j->req, (const char *)data, len) != ESP_OK) {
        return 0;
    }
    j->len += len;
    return len;
}