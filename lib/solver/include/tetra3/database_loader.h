#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "types.h"

class BinaryDatabaseLoader {
private:
    std::vector<StarEntry> star_catalog;
    std::vector<PatternEntry> pattern_catalog;

public:
    bool load_database(const std::string &stars_file, const std::string &patterns_file) {
        return load_stars(stars_file) && load_patterns(patterns_file);
    }

    const std::vector<StarEntry> &get_stars() const { return star_catalog; }
    const std::vector<PatternEntry> &get_patterns() const { return pattern_catalog; }

    void print_info() const {
        size_t star_memory = star_catalog.size() * sizeof(StarEntry);
        size_t pattern_memory = pattern_catalog.size() * sizeof(PatternEntry);
        size_t total_memory = star_memory + pattern_memory;

        std::cout << "Loaded " << star_catalog.size() << " stars, "
                << pattern_catalog.size() << " patterns" << std::endl;
        std::cout << "Memory usage: " << total_memory / 1024.0f / 1024.0f << " MB" << std::endl;
    }

private:
    bool load_stars(const std::string &filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) return false;

        uint32_t num_stars;
        file.read(reinterpret_cast<char *>(&num_stars), sizeof(uint32_t));
        if (!file) return false;

        star_catalog.resize(num_stars);
        file.read(reinterpret_cast<char *>(star_catalog.data()),
                  num_stars * sizeof(StarEntry));

        int invalid_count = 0;
        for (size_t i = 0; i < std::min(star_catalog.size(), (size_t) 10); i++) {
            const auto &star = star_catalog[i];
            if (std::isnan(star.ra) || std::isnan(star.dec) ||
                std::isnan(star.x) || std::isnan(star.y) || std::isnan(star.z)) {
                std::cout << "Invalid star " << i << ": RA=" << star.ra
                        << " Dec=" << star.dec << " xyz=(" << star.x
                        << "," << star.y << "," << star.z << ")" << std::endl;
                invalid_count++;
            }
        }

        if (invalid_count > 0) {
            std::cout << "Warning: Found " << invalid_count << " invalid stars in first 10" << std::endl;
        }

        return file.good() || file.eof();
    }

    bool load_patterns(const std::string &filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) return false;

        uint32_t num_patterns;
        file.read(reinterpret_cast<char *>(&num_patterns), sizeof(uint32_t));
        if (!file) return false;

        pattern_catalog.resize(num_patterns);
        file.read(reinterpret_cast<char *>(pattern_catalog.data()),
                  num_patterns * sizeof(PatternEntry));

        return file.good() || file.eof();
    }
};
