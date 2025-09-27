#include "camera_ops.h"
#include <esp_log.h>
#include <esp_heap_caps.h>
#include <esp_timer.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

static const char* TAG = "camera_ops";

// Board-specific pin definitions (ESP32S3_XIAO as default)
#define CAM_PIN_PWDN -1
#define CAM_PIN_RESET -1
#define CAM_PIN_VSYNC 38
#define CAM_PIN_HREF 47
#define CAM_PIN_PCLK 13
#define CAM_PIN_XCLK 10
#define CAM_PIN_SIOD 40
#define CAM_PIN_SIOC 39
#define CAM_PIN_D0 15
#define CAM_PIN_D1 17
#define CAM_PIN_D2 18
#define CAM_PIN_D3 16
#define CAM_PIN_D4 14
#define CAM_PIN_D5 12
#define CAM_PIN_D6 11
#define CAM_PIN_D7 48

static camera_config_t camera_config = {
    .pin_pwdn = CAM_PIN_PWDN,
    .pin_reset = CAM_PIN_RESET,
    .pin_xclk = CAM_PIN_XCLK,
    .pin_sccb_sda = CAM_PIN_SIOD,
    .pin_sccb_scl = CAM_PIN_SIOC,
    .pin_d7 = CAM_PIN_D7,
    .pin_d6 = CAM_PIN_D6,
    .pin_d5 = CAM_PIN_D5,
    .pin_d4 = CAM_PIN_D4,
    .pin_d3 = CAM_PIN_D3,
    .pin_d2 = CAM_PIN_D2,
    .pin_d1 = CAM_PIN_D1,
    .pin_d0 = CAM_PIN_D0,
    .pin_vsync = CAM_PIN_VSYNC,
    .pin_href = CAM_PIN_HREF,
    .pin_pclk = CAM_PIN_PCLK,

    // Astrophotography optimizations
    .xclk_freq_hz = 40000000,
    .pixel_format = PIXFORMAT_GRAYSCALE, // Monochrome for maximum light sensitivity
    .frame_size = FRAMESIZE_UXGA,         // 640x480 for good balance of detail and memory
    .jpeg_quality = 12,
    .fb_count = 1,
    .fb_location = CAMERA_FB_IN_PSRAM,
    .grab_mode = CAMERA_GRAB_WHEN_EMPTY,
};

esp_err_t camera_ops_init(void)
{
    ESP_LOGI(TAG, "Initializing camera with astrophotography settings");

    esp_err_t err = esp_camera_init(&camera_config);
    if (err != ESP_OK) {
        ESP_LOGE(TAG, "Camera init failed: %s", esp_err_to_name(err));
        return err;
    }

    return ESP_OK;

    sensor_t* s = esp_camera_sensor_get();
    if (!s) {
        ESP_LOGE(TAG, "Failed to get camera sensor");
        return ESP_FAIL;
    }

    // Apply astrophotography settings
    s->set_gain_ctrl(s, 0);      // Manual gain control
    s->set_exposure_ctrl(s, 0);  // Manual exposure control
    s->set_agc_gain(s, 2048);    // High gain for low-light sensitivity
    s->set_aec_value(s, 1048575); // Maximum exposure time
    s->set_whitebal(s, 0);       // Disable auto white balance
    s->set_lenc(s, 0);           // Disable lens correction
    s->set_denoise(s, 0);        // Minimize noise reduction

    ESP_LOGI(TAG, "Camera initialized with astrophotography settings");
    return ESP_OK;
}

void camera_ops_deinit(void)
{
    esp_camera_deinit();
    ESP_LOGI(TAG, "Camera deinitialized");
}

camera_fb_t* camera_ops_capture_frame(void)
{
    camera_fb_t* fb = esp_camera_fb_get();
    if (!fb) {
        ESP_LOGE(TAG, "Failed to capture frame");
        return NULL;
    }

    if (fb->format != PIXFORMAT_GRAYSCALE) {
        ESP_LOGW(TAG, "Frame format is not grayscale");
    }

    ESP_LOGD(TAG, "Captured frame: %dx%d, %zu bytes", fb->width, fb->height, fb->len);
    return fb;
}

mono_accumulator_t* mono_accumulator_init(size_t width, size_t height)
{
    mono_accumulator_t* accum = (mono_accumulator_t*)malloc(sizeof(mono_accumulator_t));
    if (!accum) {
        ESP_LOGE(TAG, "Failed to allocate accumulator structure");
        return NULL;
    }

    accum->width = width;
    accum->height = height;
    accum->pixel_count = width * height;
    accum->frame_count = 0;

    // Try SPIRAM first, fallback to regular RAM
    accum->accumulator = (uint16_t*)heap_caps_malloc(accum->pixel_count * sizeof(uint16_t), MALLOC_CAP_SPIRAM);
    if (!accum->accumulator) {
        ESP_LOGW(TAG, "SPIRAM allocation failed, trying regular RAM");
        accum->accumulator = (uint16_t*)malloc(accum->pixel_count * sizeof(uint16_t));
    }

    if (!accum->accumulator) {
        ESP_LOGE(TAG, "Failed to allocate accumulator buffer (%zu bytes)",
                 accum->pixel_count * sizeof(uint16_t));
        free(accum);
        return NULL;
    }

    memset(accum->accumulator, 0, accum->pixel_count * sizeof(uint16_t));

    ESP_LOGI(TAG, "Accumulator initialized: %zux%zu (%zu pixels, %zu KB)",
             width, height, accum->pixel_count,
             (accum->pixel_count * sizeof(uint16_t)) / 1024);

    return accum;
}

void mono_accumulator_add_frame(mono_accumulator_t* accum, const uint8_t* frame_data)
{
    if (!accum || !frame_data) {
        ESP_LOGE(TAG, "Invalid parameters for mono_accumulator_add_frame");
        return;
    }

    accum->frame_count++;

    for (size_t i = 0; i < accum->pixel_count; i++) {
        if (accum->accumulator[i] < UINT16_MAX) {
            accum->accumulator[i] += frame_data[i];
        }
    }

    ESP_LOGD(TAG, "Frame %" PRIu32 " added to accumulator", accum->frame_count);
}

void mono_accumulator_get_average(mono_accumulator_t* accum, uint8_t* output)
{
    if (!accum || !output || accum->frame_count == 0) {
        ESP_LOGE(TAG, "Invalid parameters or no frames accumulated");
        return;
    }

    for (size_t i = 0; i < accum->pixel_count; i++) {
        uint32_t averaged = accum->accumulator[i] / accum->frame_count;
        output[i] = (uint8_t)(averaged > 255 ? 255 : averaged);
    }

    ESP_LOGI(TAG, "Generated average image from %" PRIu32 " frames", accum->frame_count);
}

void mono_accumulator_free(mono_accumulator_t* accum)
{
    if (!accum) return;

    if (accum->accumulator) {
        free(accum->accumulator);
    }
    free(accum);
}

camera_fb_t* create_framebuffer(size_t width, size_t height, pixformat_t format)
{
    camera_fb_t* fb = (camera_fb_t*)malloc(sizeof(camera_fb_t));
    if (!fb) {
        return NULL;
    }

    fb->width = width;
    fb->height = height;
    fb->format = format;

    size_t pixel_size = 1;
    if (format == PIXFORMAT_GRAYSCALE) {
        pixel_size = 1;
    }

    fb->len = fb->width * fb->height * pixel_size;
    fb->buf = (uint8_t*)malloc(fb->len);
    if (!fb->buf) {
        free(fb);
        ESP_LOGW(TAG, "Could not allocate frame buffer");
        return NULL;
    }

    gettimeofday(&fb->timestamp, NULL);
    return fb;
}

void free_framebuffer(camera_fb_t* fb)
{
    if (!fb) return;
    if (fb->buf) free(fb->buf);
    free(fb);
}

camera_fb_t* camera_ops_capture_accumulated(int target_seconds)
{
    camera_fb_t* fb = camera_ops_capture_frame();
    if (!fb) {
        ESP_LOGE(TAG, "Failed to get initial frame for accumulation");
        return NULL;
    }

    if (fb->format != PIXFORMAT_GRAYSCALE) {
        ESP_LOGE(TAG, "Camera must be in grayscale mode for accumulation");
        esp_camera_fb_return(fb);
        return NULL;
    }

    size_t width = fb->width;
    size_t height = fb->height;
    size_t pixel_count = fb->len;

    mono_accumulator_t* accum = mono_accumulator_init(width, height);
    if (!accum) {
        ESP_LOGE(TAG, "Failed to initialize accumulator");
        esp_camera_fb_return(fb);
        return NULL;
    }

    int target_frames = target_seconds * 2; // ~2 fps target
    if (target_frames < 1) target_frames = 1;

    ESP_LOGI(TAG, "Starting accumulation: %d frames over %d seconds (%zux%zu)",
             target_frames, target_seconds, width, height);

    int64_t start_time = esp_timer_get_time();
    camera_fb_t* result_fb = NULL;

    // Capture and accumulate frames
    for (int i = 0; i < target_frames; i++) {
        camera_fb_t* current_fb = (i == 0) ? fb : camera_ops_capture_frame();

        if (current_fb && current_fb->format == PIXFORMAT_GRAYSCALE && current_fb->len == pixel_count) {
            mono_accumulator_add_frame(accum, current_fb->buf);
            ESP_LOGI(TAG, "Frame %d/%d captured and accumulated", i + 1, target_frames);

            if (i != target_frames - 1) {
                esp_camera_fb_return(current_fb);
            } else {
                result_fb = current_fb; // Keep the last frame
            }
        } else {
            ESP_LOGE(TAG, "Failed to capture frame %d", i + 1);
            if (current_fb) esp_camera_fb_return(current_fb);
            continue;
        }
    }

    int64_t end_time = esp_timer_get_time();
    float total_time_s = (end_time - start_time) / 1000000.0f;

    ESP_LOGI(TAG, "Accumulation complete: %" PRIu32 " frames in %.1f seconds (%.1f fps)",
             accum->frame_count, total_time_s, accum->frame_count / total_time_s);

    if (result_fb && accum->frame_count > 0) {
        // Write accumulated result into the last framebuffer
        mono_accumulator_get_average(accum, result_fb->buf);
        ESP_LOGI(TAG, "Generated accumulated image: %zux%zu (%zu KB)",
                 width, height, result_fb->len / 1024);
    } else {
        ESP_LOGE(TAG, "No valid frames captured for accumulation");
    }

    mono_accumulator_free(accum);
    return result_fb;
}