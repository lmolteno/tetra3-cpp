#pragma once

#include <esp_err.h>
#include <esp_camera.h>
#include <esp_http_server.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint16_t* accumulator;
    uint32_t frame_count;
    size_t pixel_count;
    size_t width;
    size_t height;
} mono_accumulator_t;

// Camera initialization and configuration
esp_err_t camera_ops_init(void);
void camera_ops_deinit(void);

// Single frame capture
camera_fb_t* camera_ops_capture_frame(void);

// Accumulated frame capture for astrophotography
camera_fb_t* camera_ops_capture_accumulated(int target_seconds);

// Frame stacking utilities
mono_accumulator_t* mono_accumulator_init(size_t width, size_t height);
void mono_accumulator_add_frame(mono_accumulator_t* accum, const uint8_t* frame_data);
void mono_accumulator_get_average(mono_accumulator_t* accum, uint8_t* output);
void mono_accumulator_free(mono_accumulator_t* accum);

// Frame buffer utilities
camera_fb_t* create_framebuffer(size_t width, size_t height, pixformat_t format);
void free_framebuffer(camera_fb_t* fb);

#ifdef __cplusplus
}
#endif