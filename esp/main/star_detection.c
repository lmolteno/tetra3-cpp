#include "star_detection.h"
#include <esp_log.h>
#include <esp_heap_caps.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static const char* TAG = "star_detection";

typedef struct {
    int y, x;
} point_t;

typedef struct {
    point_t* points;
    int count;
    int capacity;
} cluster_t;

// Stack for flood fill algorithm
typedef struct {
    point_t* data;
    int size;
    int capacity;
} flood_stack_t;

static flood_stack_t* stack_create(int initial_capacity)
{
    flood_stack_t* stack = malloc(sizeof(flood_stack_t));
    if (!stack) return NULL;

    stack->data = malloc(initial_capacity * sizeof(point_t));
    if (!stack->data) {
        free(stack);
        return NULL;
    }

    stack->size = 0;
    stack->capacity = initial_capacity;
    return stack;
}

static void stack_destroy(flood_stack_t* stack)
{
    if (stack) {
        if (stack->data) free(stack->data);
        free(stack);
    }
}

static bool stack_push(flood_stack_t* stack, point_t p)
{
    if (stack->size >= stack->capacity) {
        // Resize if needed
        int new_capacity = stack->capacity * 2;
        point_t* new_data = realloc(stack->data, new_capacity * sizeof(point_t));
        if (!new_data) return false;

        stack->data = new_data;
        stack->capacity = new_capacity;
    }

    stack->data[stack->size++] = p;
    return true;
}

static bool stack_pop(flood_stack_t* stack, point_t* p)
{
    if (stack->size == 0) return false;
    *p = stack->data[--stack->size];
    return true;
}

// Calculate median of uint8_t array
static uint8_t calculate_median(const uint8_t* data, size_t len)
{
    // For large arrays, we'll use a histogram approach for efficiency
    uint32_t histogram[256] = {0};

    // Build histogram
    for (size_t i = 0; i < len; i++) {
        histogram[data[i]]++;
    }

    // Find median
    size_t target = len / 2;
    size_t count = 0;

    for (int i = 0; i < 256; i++) {
        count += histogram[i];
        if (count >= target) {
            return (uint8_t)i;
        }
    }

    return 128; // fallback
}

// Flood fill to find connected component
static int flood_fill(const uint8_t* binary_mask, int width, int height,
                     int start_y, int start_x, uint8_t* visited,
                     point_t* cluster_points, int max_cluster_size)
{
    flood_stack_t* stack = stack_create(256);
    if (!stack) return 0;

    point_t start = {start_y, start_x};
    if (!stack_push(stack, start)) {
        stack_destroy(stack);
        return 0;
    }

    int cluster_size = 0;

    while (stack->size > 0) {
        point_t p;
        if (!stack_pop(stack, &p)) break;

        int y = p.y;
        int x = p.x;

        // Check bounds and visited status
        if (y < 0 || y >= height || x < 0 || x >= width) continue;
        if (visited[y * width + x]) continue;
        if (!binary_mask[y * width + x]) continue;

        // Mark as visited
        visited[y * width + x] = 1;

        // Add to cluster
        if (cluster_size < max_cluster_size) {
            cluster_points[cluster_size].y = y;
            cluster_points[cluster_size].x = x;
        }
        cluster_size++;

        // Add 8-connected neighbors to stack
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dy == 0 && dx == 0) continue;

                point_t neighbor = {y + dy, x + dx};
                stack_push(stack, neighbor);
            }
        }
    }

    stack_destroy(stack);
    return cluster_size;
}

star_detection_result_t* detect_stars_simple(const camera_fb_t* fb,
                                            float threshold,
                                            int min_area,
                                            int max_stars)
{
    if (!fb || fb->format != PIXFORMAT_GRAYSCALE) {
        ESP_LOGE(TAG, "Invalid frame buffer or not grayscale");
        return NULL;
    }

    int width = fb->width;
    int height = fb->height;
    size_t pixel_count = width * height;

    ESP_LOGI(TAG, "Starting star detection on %dx%d image", width, height);

    // Step 1: Calculate median
    uint8_t median = calculate_median(fb->buf, pixel_count);
    ESP_LOGD(TAG, "Image median: %d", median);

    // Step 2: Median subtraction and rescaling
    uint8_t* processed = malloc(pixel_count);
    if (!processed) {
        ESP_LOGE(TAG, "Failed to allocate processed image buffer");
        return NULL;
    }

    // Find max value after median subtraction
    int max_val = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        int val = (int)fb->buf[i] - median;
        if (val > max_val) max_val = val;
    }

    if (max_val == 0) {
        ESP_LOGW(TAG, "No pixels above median, no stars detected");
        free(processed);
        return NULL;
    }

    // Rescale to 0-255
    float scale = 255.0f / max_val;
    for (size_t i = 0; i < pixel_count; i++) {
        int val = (int)fb->buf[i] - median;
        if (val <= 0) {
            processed[i] = 0;
        } else {
            processed[i] = (uint8_t)(val * scale);
        }
    }

    ESP_LOGD(TAG, "Rescaled image, scale factor: %.2f", scale);

    // Step 3: Create binary mask
    uint8_t* binary_mask = malloc(pixel_count);
    if (!binary_mask) {
        ESP_LOGE(TAG, "Failed to allocate binary mask");
        free(processed);
        return NULL;
    }

    int threshold_int = (int)threshold;
    for (size_t i = 0; i < pixel_count; i++) {
        binary_mask[i] = (processed[i] > threshold_int) ? 1 : 0;
    }

    free(processed); // No longer needed

    // Step 4: Find connected components using flood fill
    uint8_t* visited = calloc(pixel_count, 1);
    if (!visited) {
        ESP_LOGE(TAG, "Failed to allocate visited array");
        free(binary_mask);
        return NULL;
    }

    // Temporary storage for cluster points
    point_t* temp_cluster = malloc(min_area * 4 * sizeof(point_t)); // reasonable max cluster size
    if (!temp_cluster) {
        ESP_LOGE(TAG, "Failed to allocate temporary cluster buffer");
        free(binary_mask);
        free(visited);
        return NULL;
    }

    // Storage for results
    centroid_t* centroids = malloc(max_stars * sizeof(centroid_t));
    float* sizes = malloc(max_stars * sizeof(float));
    if (!centroids || !sizes) {
        ESP_LOGE(TAG, "Failed to allocate result arrays");
        free(binary_mask);
        free(visited);
        free(temp_cluster);
        if (centroids) free(centroids);
        if (sizes) free(sizes);
        return NULL;
    }

    int num_centroids = 0;

    // Scan image for star clusters
    for (int y = 0; y < height && num_centroids < max_stars; y++) {
        for (int x = 0; x < width && num_centroids < max_stars; x++) {
            if (binary_mask[y * width + x] && !visited[y * width + x]) {
                int cluster_size = flood_fill(binary_mask, width, height, y, x, visited,
                                            temp_cluster, min_area * 4);

                if (cluster_size >= min_area) {
                    // Calculate centroid
                    float sum_y = 0, sum_x = 0;
                    int points_to_use = (cluster_size < min_area * 4) ? cluster_size : min_area * 4;

                    for (int i = 0; i < points_to_use; i++) {
                        sum_y += temp_cluster[i].y;
                        sum_x += temp_cluster[i].x;
                    }

                    centroids[num_centroids].y = sum_y / points_to_use;
                    centroids[num_centroids].x = sum_x / points_to_use;
                    sizes[num_centroids] = (float)cluster_size;

                    ESP_LOGD(TAG, "Star %d: center=(%.1f, %.1f), size=%d",
                            num_centroids, centroids[num_centroids].x,
                            centroids[num_centroids].y, cluster_size);

                    num_centroids++;
                }
            }
        }
    }

    // Sort stars by size (biggest first) using simple bubble sort
    for (int i = 0; i < num_centroids - 1; i++) {
        for (int j = 0; j < num_centroids - 1 - i; j++) {
            if (sizes[j] < sizes[j + 1]) {
                // Swap sizes
                float temp_size = sizes[j];
                sizes[j] = sizes[j + 1];
                sizes[j + 1] = temp_size;

                // Swap corresponding centroids
                centroid_t temp_centroid = centroids[j];
                centroids[j] = centroids[j + 1];
                centroids[j + 1] = temp_centroid;
            }
        }
    }

    // Cleanup temporary arrays
    free(binary_mask);
    free(visited);
    free(temp_cluster);

    // Create result structure
    star_detection_result_t* result = malloc(sizeof(star_detection_result_t));
    if (!result) {
        ESP_LOGE(TAG, "Failed to allocate result structure");
        free(centroids);
        free(sizes);
        return NULL;
    }

    result->num_centroids = num_centroids;
    result->centroids = centroids;
    result->sizes = sizes;

    ESP_LOGI(TAG, "Star detection complete: found %d stars (sorted by size)", num_centroids);

    // Log the top few stars for debugging
    for (int i = 0; i < num_centroids && i < 5; i++) {
        ESP_LOGI(TAG, "Star %d: center=(%.1f, %.1f), size=%.0f",
                i + 1, centroids[i].x, centroids[i].y, sizes[i]);
    }

    return result;
}

star_detection_result_t* detect_stars_test(const camera_fb_t* fb)
{
    if (!fb) {
        ESP_LOGE(TAG, "Invalid frame buffer");
        return NULL;
    }

    // Hardcoded star positions from 1600x1200 image
    const float reference_stars[][2] = {
        {692.4872f, 491.22507f},
        {946.2311f, 332.36243f},
        {12.089273f, 442.55933f},
        {943.6792f, 298.08932f},
        {734.17804f, 297.76077f},
        {1335.9685f, 518.3289f},
        {1374.8608f, 493.86456f},
        {905.30676f, 585.91943f},
        {1098.4895f, 316.72595f},
        {419.94678f, 1107.7972f}
    };
    const int num_reference_stars = sizeof(reference_stars) / sizeof(reference_stars[0]);

    // Calculate scaling factors
    const float reference_width = 1600.0f;
    const float reference_height = 1200.0f;
    const float scale_x = (float)fb->width / reference_width;
    const float scale_y = (float)fb->height / reference_height;

    ESP_LOGI(TAG, "Test star detection: scaling from %dx%d to %dx%d (scale: %.3f, %.3f)",
             (int)reference_width, (int)reference_height, fb->width, fb->height, scale_x, scale_y);

    // Allocate result structure
    star_detection_result_t* result = malloc(sizeof(star_detection_result_t));
    if (!result) {
        ESP_LOGE(TAG, "Failed to allocate result structure");
        return NULL;
    }

    result->num_centroids = num_reference_stars;
    result->centroids = malloc(num_reference_stars * sizeof(centroid_t));
    result->sizes = malloc(num_reference_stars * sizeof(float));

    if (!result->centroids || !result->sizes) {
        ESP_LOGE(TAG, "Failed to allocate centroid arrays");
        star_detection_free(result);
        return NULL;
    }

    // Scale and copy the reference stars
    for (int i = 0; i < num_reference_stars; i++) {
        result->centroids[i].x = reference_stars[i][1] * scale_x;
        result->centroids[i].y = reference_stars[i][0] * scale_y;
        result->sizes[i] = 100.0f;  // Fixed size for test stars

        ESP_LOGI(TAG, "Test star %d: center=(%.1f, %.1f), size=%.0f",
                i + 1, result->centroids[i].x, result->centroids[i].y, result->sizes[i]);
    }

    ESP_LOGI(TAG, "Generated %d test stars", num_reference_stars);
    return result;
}

void star_detection_free(star_detection_result_t* result)
{
    if (!result) return;

    if (result->centroids) free(result->centroids);
    if (result->sizes) free(result->sizes);
    free(result);
}