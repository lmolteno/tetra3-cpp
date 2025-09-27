#include "http_server.h"
#include "camera_ops.h"
#include "star_solver.h"
#include "star_detection.h"
#include "star_catalog.h"
#include "pattern_catalog.h"
// #include "example_catalogs.h"
#include <esp_log.h>
#include <esp_timer.h>
#include <string.h>
#include <stdlib.h>

static const char* TAG = "http_server";

// MJPEG stream constants
#define PART_BOUNDARY "123456789000000000000987654321"
static const char* _STREAM_CONTENT_TYPE = "multipart/x-mixed-replace;boundary=" PART_BOUNDARY;
static const char* _STREAM_BOUNDARY = "\r\n--" PART_BOUNDARY "\r\n";
static const char* _STREAM_PART = "Content-Type: image/jpeg\r\nContent-Length: %zu\r\n\r\n";

// Forward declarations of handler functions
static esp_err_t capture_handler(httpd_req_t *req);
static esp_err_t accumulate_handler(httpd_req_t *req);
static esp_err_t detect_handler(httpd_req_t *req);
static esp_err_t solve_handler(httpd_req_t *req);
static esp_err_t capture_solve_handler(httpd_req_t *req);
static esp_err_t health_handler(httpd_req_t *req);
static esp_err_t stream_handler(httpd_req_t *req);

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
    config.max_uri_handlers = 10;
    config.server_port = 80;
    config.stack_size = 8192;  // Increase from default 4096 to 8192

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

    httpd_uri_t capture_solve_uri = {
        .uri = "/capture-solve",
        .method = HTTP_GET,
        .handler = capture_solve_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &capture_solve_uri);

    httpd_uri_t stream_uri = {
        .uri = "/stream",
        .method = HTTP_GET,
        .handler = stream_handler,
        .user_ctx = NULL
    };
    httpd_register_uri_handler(server, &stream_uri);

    ESP_LOGI(TAG, "Available endpoints:");
    ESP_LOGI(TAG, "  GET /health - Server health check");
    ESP_LOGI(TAG, "  GET /capture - Single frame JPEG");
    ESP_LOGI(TAG, "  GET /accumulate?seconds=N - Accumulated JPEG (1-60s)");
    ESP_LOGI(TAG, "  GET /detect?seconds=N&threshold=X&test=1 - Star detection only");
    ESP_LOGI(TAG, "  POST /solve - Star solver endpoint");
    ESP_LOGI(TAG, "  GET /capture-solve?seconds=N&threshold=X&test=1 - Capture, detect and solve");
    ESP_LOGI(TAG, "  GET /stream - MJPEG stream");

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
    bool use_test_stars = false;

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
        if (httpd_query_key_value(query, "test", param, sizeof(param)) == ESP_OK) {
            use_test_stars = (strcmp(param, "1") == 0 || strcmp(param, "true") == 0);
        }
    }

    if (use_test_stars) {
        ESP_LOGI(TAG, "Starting TEST detection with hardcoded stars");
    } else {
        ESP_LOGI(TAG, "Starting detection: %ds exposure, threshold=%.1f, min_area=%d",
                 seconds, threshold, min_area);
    }

    int64_t start_time = esp_timer_get_time();

    // Capture accumulated frame
    camera_fb_t* fb = camera_ops_capture_accumulated(seconds);
    if (!fb) {
        ESP_LOGE(TAG, "Failed to capture frame for detection");
        httpd_resp_send_500(req);
        return ESP_FAIL;
    }

    // Detect stars using our algorithm or test data
    star_detection_result_t* detection;
    if (use_test_stars) {
        detection = detect_stars_test(fb);
    } else {
        detection = detect_stars_simple(fb, threshold, min_area, 100);
    }

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

    ESP_LOGI(TAG, "Instantiated solver");

    char response[512];
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
            7.0,                     // fov_estimate_deg
            16,                        // pattern_checking_stars
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

static esp_err_t capture_solve_handler(httpd_req_t *req)
{
    ESP_LOGI(TAG, "Capture and solve request");

    char query[128];
    int seconds = 5;
    float threshold = 28.0f;
    int min_area = 3;
    bool use_test_stars = false;

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
        if (httpd_query_key_value(query, "test", param, sizeof(param)) == ESP_OK) {
            use_test_stars = (strcmp(param, "1") == 0 || strcmp(param, "true") == 0);
        }
    }

    if (use_test_stars) {
        ESP_LOGI(TAG, "Starting capture-solve with TEST stars");
    } else {
        ESP_LOGI(TAG, "Starting capture-solve: %ds exposure, threshold=%.1f, min_area=%d",
                 seconds, threshold, min_area);
    }

    int64_t start_time = esp_timer_get_time();

    // Capture accumulated frame
    camera_fb_t* fb = camera_ops_capture_accumulated(seconds);
    if (!fb) {
        ESP_LOGE(TAG, "Failed to capture frame for capture-solve");
        httpd_resp_send_500(req);
        return ESP_FAIL;
    }

    // Detect stars using our algorithm or test data
    star_detection_result_t* detection;
    if (use_test_stars) {
        detection = detect_stars_test(fb);
    } else {
        detection = detect_stars_simple(fb, threshold, min_area, 100);
    }
    if (!detection || detection->num_centroids == 0) {
        ESP_LOGW(TAG, "No stars detected in image");

        char response[256];
        snprintf(response, sizeof(response),
            "{"
            "\"status\":\"no_stars\","
            "\"width\":%d,"
            "\"height\":%d,"
            "\"exposure_seconds\":%d,"
            "\"threshold\":%.1f,"
            "\"min_area\":%d,"
            "\"stars_detected\":0"
            "}",
            fb->width, fb->height, seconds, threshold, min_area);

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

    esp_err_t res;
    if (!solver) {
        ESP_LOGE(TAG, "Failed to create star solver");

        const char* error_response =
            "{"
            "\"status\":\"solver_error\","
            "\"stars_detected\":%d"
            "}";

        char response[128];
        snprintf(response, sizeof(response), error_response, detection->num_centroids);

        esp_camera_fb_return(fb);
        star_detection_free(detection);

        res = httpd_resp_set_type(req, "application/json");
        if (res == ESP_OK) {
            res = httpd_resp_send(req, response, strlen(response));
        }

        int64_t end_time = esp_timer_get_time();
        ESP_LOGI(TAG, "Capture-solve request completed in %lld ms", (end_time - start_time) / 1000);
        return res;
    }

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

    // Store star count before cleanup
    int num_stars = detection->num_centroids;

    // Cleanup early to free memory
    esp_camera_fb_return(fb);
    star_detection_free(detection);
    star_solver_destroy(solver);

    // Build response with smaller buffer
    char response[384];
    if (solve_result.solved) {
        ESP_LOGI(TAG, "SOLVED! RA: %.2f°, Dec: %.2f°, Roll: %.2f°",
                 solve_result.ra * 180.0f / 3.14159f,
                 solve_result.dec * 180.0f / 3.14159f,
                 solve_result.roll * 180.0f / 3.14159f);

        snprintf(response, sizeof(response),
            "{"
            "\"status\":\"solved\","
            "\"stars_detected\":%d,"
            "\"ra_deg\":%.6f,"
            "\"dec_deg\":%.6f,"
            "\"roll_deg\":%.6f,"
            "\"fov_deg\":%.6f,"
            "\"rmse_arcsec\":%.2f,"
            "\"num_matches\":%d,"
            "\"solve_time_ms\":%.2f"
            "}",
            num_stars,
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
            "\"stars_detected\":%d,"
            "\"solve_time_ms\":%.2f"
            "}",
            num_stars,
            solve_result.solve_time_ms);
    }

    res = httpd_resp_set_type(req, "application/json");
    if (res == ESP_OK) {
        res = httpd_resp_send(req, response, strlen(response));
    }

    int64_t end_time = esp_timer_get_time();
    ESP_LOGI(TAG, "Capture-solve request completed in %lld ms", (end_time - start_time) / 1000);

    return res;
}

static esp_err_t stream_handler(httpd_req_t *req)
{
    ESP_LOGI(TAG, "MJPEG stream started");

    camera_fb_t *fb = NULL;
    esp_err_t res = ESP_OK;
    size_t jpg_buf_len = 0;
    uint8_t *jpg_buf = NULL;
    char part_buf[64];
    static int64_t last_frame = 0;

    if (!last_frame) {
        last_frame = esp_timer_get_time();
    }

    res = httpd_resp_set_type(req, _STREAM_CONTENT_TYPE);
    if (res != ESP_OK) {
        return res;
    }

    while (true) {
        fb = camera_ops_capture_frame();
        if (!fb) {
            ESP_LOGE(TAG, "Camera capture failed in stream");
            res = ESP_FAIL;
            break;
        }

        if (fb->format != PIXFORMAT_JPEG) {
            // Convert to JPEG
            jpg_chunking_t jchunk = {req, 0};
            bool jpeg_converted = frame2jpg_cb(fb, 80, jpg_encode_stream, &jchunk);
            esp_camera_fb_return(fb);

            if (!jpeg_converted) {
                ESP_LOGE(TAG, "JPEG compression failed in stream");
                res = ESP_FAIL;
                break;
            }
            jpg_buf_len = jchunk.len;
        } else {
            jpg_buf_len = fb->len;
            jpg_buf = fb->buf;
        }

        // Send boundary
        if (res == ESP_OK) {
            res = httpd_resp_send_chunk(req, _STREAM_BOUNDARY, strlen(_STREAM_BOUNDARY));
        }

        // Send part header
        if (res == ESP_OK) {
            int hlen = snprintf(part_buf, sizeof(part_buf), _STREAM_PART, jpg_buf_len);
            if (hlen < 0 || hlen >= sizeof(part_buf)) {
                ESP_LOGE(TAG, "Header truncated (%d bytes needed >= %zu buffer)", hlen, sizeof(part_buf));
                res = ESP_FAIL;
            } else {
                res = httpd_resp_send_chunk(req, part_buf, (size_t)hlen);
            }
        }

        // Send JPEG data
        if (res == ESP_OK) {
            if (fb->format == PIXFORMAT_JPEG) {
                res = httpd_resp_send_chunk(req, (const char *)jpg_buf, jpg_buf_len);
                esp_camera_fb_return(fb);
            } else {
                // For non-JPEG formats, data was already sent via jpg_encode_stream callback
                httpd_resp_send_chunk(req, NULL, 0); // End chunk
            }
        }

        if (res != ESP_OK) {
            if (fb && fb->format == PIXFORMAT_JPEG) {
                esp_camera_fb_return(fb);
            }
            break;
        }

        // Calculate FPS
        int64_t fr_end = esp_timer_get_time();
        int64_t frame_time = fr_end - last_frame;
        last_frame = fr_end;
        frame_time /= 1000;
        float fps = frame_time > 0 ? 1000.0f / (float)frame_time : 0.0f;

        ESP_LOGI(TAG, "MJPEG: %luKB %lums (%.1ffps)",
                 (uint32_t)(jpg_buf_len/1024),
                 (uint32_t)frame_time, fps);
    }

    last_frame = 0;
    ESP_LOGI(TAG, "MJPEG stream ended");
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