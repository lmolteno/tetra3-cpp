#pragma once

#include <Arduino.h>
#include <FS.h>
#include <SPIFFS.h>

// Star catalog structures
struct StarCatalogEntry {
    float ra, dec;
    float x, y, z;
    float magnitude;
};

struct PatternCatalogEntry {
    uint16_t star_indices[4];
};

class CatalogLoader {
private:
    StarCatalogEntry* stars = nullptr;
    PatternCatalogEntry* patterns = nullptr;
    uint32_t num_stars = 0;
    uint32_t num_patterns = 0;
    bool loaded = false;

public:
    CatalogLoader() {}
    
    ~CatalogLoader() {
        cleanup();
    }
    
    // Load catalogs from SPIFFS
    bool load_from_spiffs(const char* star_file = "/stars.bin", 
                         const char* pattern_file = "/patterns.bin");
    
    // Load embedded catalogs (compiled into flash)
    bool load_embedded();
    
    // Get catalog data
    const StarCatalogEntry* get_stars(uint32_t* count = nullptr) const;
    const PatternCatalogEntry* get_patterns(uint32_t* count = nullptr) const;
    
    // Check if catalogs are loaded
    bool is_loaded() const { return loaded; }
    
    // Get memory usage
    size_t get_memory_usage() const;
    
    // Cleanup
    void cleanup();
    
    // Print catalog info
    void print_info() const;
};

// Global catalog loader instance
extern CatalogLoader catalog;