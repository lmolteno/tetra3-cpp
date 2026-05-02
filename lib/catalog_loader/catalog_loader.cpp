#include "catalog_loader.h"

// Global instance
CatalogLoader catalog;

bool CatalogLoader::load_from_spiffs(const char* star_file, const char* pattern_file) {
    Serial.println("Loading catalogs from SPIFFS...");
    
    if (!SPIFFS.begin()) {
        Serial.println("SPIFFS mount failed");
        return false;
    }
    
    // Load stars
    File star_fs = SPIFFS.open(star_file, "r");
    if (!star_fs) {
        Serial.printf("Failed to open %s\n", star_file);
        return false;
    }
    
    // Read star count
    if (star_fs.read((uint8_t*)&num_stars, sizeof(uint32_t)) != sizeof(uint32_t)) {
        Serial.println("Failed to read star count");
        star_fs.close();
        return false;
    }
    
    Serial.printf("Loading %d stars from SPIFFS...\n", num_stars);
    
    // Allocate memory for stars
    size_t star_size = num_stars * sizeof(StarCatalogEntry);
    stars = (StarCatalogEntry*)malloc(star_size);
    if (!stars) {
        Serial.printf("Failed to allocate %d bytes for stars\n", star_size);
        star_fs.close();
        return false;
    }
    
    // Read star data
    size_t read_size = star_fs.read((uint8_t*)stars, star_size);
    star_fs.close();
    
    if (read_size != star_size) {
        Serial.printf("Star read size mismatch: %d != %d\n", read_size, star_size);
        free(stars);
        stars = nullptr;
        return false;
    }
    
    // Load patterns
    File pattern_fs = SPIFFS.open(pattern_file, "r");
    if (!pattern_fs) {
        Serial.printf("Failed to open %s\n", pattern_file);
        free(stars);
        stars = nullptr;
        return false;
    }
    
    // Read pattern count
    if (pattern_fs.read((uint8_t*)&num_patterns, sizeof(uint32_t)) != sizeof(uint32_t)) {
        Serial.println("Failed to read pattern count");
        pattern_fs.close();
        free(stars);
        stars = nullptr;
        return false;
    }
    
    Serial.printf("Loading %d patterns from SPIFFS...\n", num_patterns);
    
    // Allocate memory for patterns
    size_t pattern_size = num_patterns * sizeof(PatternCatalogEntry);
    patterns = (PatternCatalogEntry*)malloc(pattern_size);
    if (!patterns) {
        Serial.printf("Failed to allocate %d bytes for patterns\n", pattern_size);
        pattern_fs.close();
        free(stars);
        stars = nullptr;
        return false;
    }
    
    // Read pattern data
    read_size = pattern_fs.read((uint8_t*)patterns, pattern_size);
    pattern_fs.close();
    
    if (read_size != pattern_size) {
        Serial.printf("Pattern read size mismatch: %d != %d\n", read_size, pattern_size);
        cleanup();
        return false;
    }
    
    loaded = true;
    print_info();
    return true;
}

bool CatalogLoader::load_embedded() {
    Serial.println("Loading embedded test catalog...");
    
    // For now, create a small test catalog
    num_stars = 50;
    num_patterns = 25;
    
    // Allocate memory
    stars = (StarCatalogEntry*)malloc(num_stars * sizeof(StarCatalogEntry));
    patterns = (PatternCatalogEntry*)malloc(num_patterns * sizeof(PatternCatalogEntry));
    
    if (!stars || !patterns) {
        cleanup();
        return false;
    }
    
    // Generate test stars
    for (uint32_t i = 0; i < num_stars; i++) {
        stars[i].ra = (float)i * 0.2f;
        stars[i].dec = sin((float)i * 0.3f) * 1.0f;
        stars[i].x = cos(stars[i].dec) * cos(stars[i].ra);
        stars[i].y = cos(stars[i].dec) * sin(stars[i].ra);
        stars[i].z = sin(stars[i].dec);
        stars[i].magnitude = 4.0f + (float)(i % 6);
    }
    
    // Generate test patterns
    for (uint32_t i = 0; i < num_patterns; i++) {
        patterns[i].star_indices[0] = (i + 0) % num_stars;
        patterns[i].star_indices[1] = (i + 1) % num_stars;
        patterns[i].star_indices[2] = (i + 2) % num_stars;
        patterns[i].star_indices[3] = (i + 3) % num_stars;
    }
    
    loaded = true;
    print_info();
    return true;
}

const StarCatalogEntry* CatalogLoader::get_stars(uint32_t* count) const {
    if (count) *count = num_stars;
    return stars;
}

const PatternCatalogEntry* CatalogLoader::get_patterns(uint32_t* count) const {
    if (count) *count = num_patterns;
    return patterns;
}

size_t CatalogLoader::get_memory_usage() const {
    size_t star_mem = num_stars * sizeof(StarCatalogEntry);
    size_t pattern_mem = num_patterns * sizeof(PatternCatalogEntry);
    return star_mem + pattern_mem;
}

void CatalogLoader::cleanup() {
    if (stars) {
        free(stars);
        stars = nullptr;
    }
    if (patterns) {
        free(patterns);
        patterns = nullptr;
    }
    num_stars = 0;
    num_patterns = 0;
    loaded = false;
}

void CatalogLoader::print_info() const {
    size_t total_mem = get_memory_usage();
    Serial.printf("Catalog loaded successfully:\n");
    Serial.printf("  Stars: %d entries (%.1f KB)\n", 
                  num_stars, (num_stars * sizeof(StarCatalogEntry)) / 1024.0f);
    Serial.printf("  Patterns: %d entries (%.1f KB)\n", 
                  num_patterns, (num_patterns * sizeof(PatternCatalogEntry)) / 1024.0f);
    Serial.printf("  Total: %.1f KB\n", total_mem / 1024.0f);
    Serial.printf("  Free heap: %d KB\n", ESP.getFreeHeap() / 1024);
}