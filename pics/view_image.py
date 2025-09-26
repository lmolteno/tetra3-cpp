import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

# Load and display the image
img = mpimg.imread('accumulated.jpg')
if len(img.shape) == 3:
    img = np.mean(img, axis=2)  # Convert to grayscale if needed

# Perform median subtraction
median_val = np.median(img)
img_median_sub = img - median_val

# Clip to above zero and scale to 0-255
img_clipped = np.clip(img_median_sub, 0, None)
img_scaled = (img_clipped / img_clipped.max()) * 255

# Create binary mask for pixels above 28
binary_mask = img_scaled > 28

# Simple flood fill to find connected components
def flood_fill(mask, start_y, start_x, visited):
    """Simple flood fill using a stack (no recursion for ESP32)"""
    h, w = mask.shape
    stack = [(start_y, start_x)]
    component = []

    while stack:
        y, x = stack.pop()
        if (y < 0 or y >= h or x < 0 or x >= w or
            visited[y, x] or not mask[y, x]):
            continue

        visited[y, x] = True
        component.append((y, x))

        # Check 8 neighbors
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dy == 0 and dx == 0:
                    continue
                stack.append((y + dy, x + dx))

    return component

# Find all star clusters
visited = np.zeros_like(binary_mask, dtype=bool)
clusters = []

h, w = binary_mask.shape
for y in range(h):
    for x in range(w):
        if binary_mask[y, x] and not visited[y, x]:
            cluster = flood_fill(binary_mask, y, x, visited)
            if len(cluster) >= 3:  # Minimum cluster size
                clusters.append(cluster)

# Calculate centroids (simple average)
centroids = []
for cluster in clusters:
    cy = sum(p[0] for p in cluster) / len(cluster)
    cx = sum(p[1] for p in cluster) / len(cluster)
    centroids.append((cy, cx, len(cluster)))

print(f"Found {len(clusters)} star clusters")
for i, (cy, cx, size) in enumerate(centroids):
    print(f"Star {i+1}: center=({cx:.1f}, {cy:.1f}), size={size} pixels")

# Visualize
plt.figure(figsize=(10, 6))
plt.imshow(binary_mask, cmap='gray')
plt.title(f'Binary Mask with {len(centroids)} Star Centers')

# Plot centroids
for cy, cx, size in centroids:
    plt.plot(cx, cy, 'r+', markersize=10, markeredgewidth=2)
    plt.text(cx+5, cy, f'{size}', color='red', fontsize=8)

plt.axis('off')
plt.tight_layout()
plt.show()
