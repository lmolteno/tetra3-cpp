#!/usr/bin/env python3
"""
Binary catalog to C header converter for ESP32 star solver.

Converts tetra3 binary star and pattern catalogs to C header files
that can be embedded in ESP32 flash memory.

Usage:
    python generate_catalog_headers.py

Input files:
    ../tracker/tetra3_db_stars.bin
    ../tracker/tetra3_db_patterns.bin

Output files:
    esp/main/star_catalog.h
    esp/main/pattern_catalog.h
"""

import struct
import os
from pathlib import Path

def read_binary_stars(filename):
    """Read binary star catalog file."""
    print(f"Reading star catalog: {filename}")

    with open(filename, 'rb') as f:
        # Read number of stars
        num_stars_data = f.read(4)
        if len(num_stars_data) != 4:
            raise ValueError("Could not read star count")

        num_stars = struct.unpack('<I', num_stars_data)[0]
        print(f"  Number of stars: {num_stars}")

        # Read star entries - each StarEntry is 24 bytes:
        # float ra, dec (4 bytes each)
        # float x, y, z (4 bytes each)
        # float magnitude (4 bytes)
        stars = []
        expected_size = num_stars * 24
        star_data = f.read(expected_size)

        if len(star_data) != expected_size:
            raise ValueError(f"Expected {expected_size} bytes of star data, got {len(star_data)}")

        for i in range(num_stars):
            offset = i * 24
            star_bytes = star_data[offset:offset + 24]
            ra, dec, x, y, z, magnitude = struct.unpack('<ffffff', star_bytes)
            stars.append({
                'ra': ra,
                'dec': dec,
                'x': x,
                'y': y,
                'z': z,
                'magnitude': magnitude
            })

        print(f"  Successfully read {len(stars)} stars")
        return stars

def read_binary_patterns(filename):
    """Read binary pattern catalog file."""
    print(f"Reading pattern catalog: {filename}")

    with open(filename, 'rb') as f:
        # Read number of patterns
        num_patterns_data = f.read(4)
        if len(num_patterns_data) != 4:
            raise ValueError("Could not read pattern count")

        num_patterns = struct.unpack('<I', num_patterns_data)[0]
        print(f"  Number of patterns: {num_patterns}")

        # Read pattern entries - each PatternEntry is 8 bytes:
        # uint16_t star_indices[4] (2 bytes each)
        patterns = []
        expected_size = num_patterns * 8
        pattern_data = f.read(expected_size)

        if len(pattern_data) != expected_size:
            raise ValueError(f"Expected {expected_size} bytes of pattern data, got {len(pattern_data)}")

        for i in range(num_patterns):
            offset = i * 8
            pattern_bytes = pattern_data[offset:offset + 8]
            indices = struct.unpack('<HHHH', pattern_bytes)
            patterns.append({
                'star_indices': list(indices)
            })

        print(f"  Successfully read {len(patterns)} patterns")
        return patterns

def write_star_header(stars, output_file):
    """Write star catalog as C header file."""
    print(f"Writing star header: {output_file}")

    with open(output_file, 'w') as f:
        f.write('#pragma once\n\n')
        f.write('#include "star_solver.h"\n\n')
        f.write('// Auto-generated star catalog header\n')
        f.write('// DO NOT EDIT MANUALLY\n\n')
        f.write(f'#define STAR_CATALOG_SIZE {len(stars)}\n\n')
        f.write('static const star_entry_t star_catalog[] = {\n')

        # Write stars in chunks to avoid overly long lines
        for i, star in enumerate(stars):
            f.write(f'    {{{star["ra"]:.6f}f, {star["dec"]:.6f}f, '
                   f'{star["x"]:.6f}f, {star["y"]:.6f}f, {star["z"]:.6f}f, '
                   f'{star["magnitude"]:.3f}f}}')

            if i < len(stars) - 1:
                f.write(',')
            f.write('\n')

            # Add progress indication for large arrays
            if (i + 1) % 10000 == 0:
                print(f"  Written {i + 1}/{len(stars)} stars")

        f.write('};\n\n')
        f.write(f'// Total memory usage: {len(stars) * 24} bytes ({len(stars) * 24 / 1024:.1f} KB)\n')

    print(f"  Successfully wrote {len(stars)} stars to header")

def write_pattern_header(patterns, output_file):
    """Write pattern catalog as C header file."""
    print(f"Writing pattern header: {output_file}")

    with open(output_file, 'w') as f:
        f.write('#pragma once\n\n')
        f.write('#include "star_solver.h"\n\n')
        f.write('// Auto-generated pattern catalog header\n')
        f.write('// DO NOT EDIT MANUALLY\n\n')
        f.write(f'#define PATTERN_CATALOG_SIZE {len(patterns)}\n\n')
        f.write('static const pattern_entry_t pattern_catalog[] = {\n')

        # Write patterns in chunks
        for i, pattern in enumerate(patterns):
            indices = pattern['star_indices']
            f.write(f'    {{{{{indices[0]}, {indices[1]}, {indices[2]}, {indices[3]}}}}}')

            if i < len(patterns) - 1:
                f.write(',')
            f.write('\n')

            # Add progress indication for large arrays
            if (i + 1) % 50000 == 0:
                print(f"  Written {i + 1}/{len(patterns)} patterns")

        f.write('};\n\n')
        f.write(f'// Total memory usage: {len(patterns) * 8} bytes ({len(patterns) * 8 / 1024:.1f} KB)\n')

    print(f"  Successfully wrote {len(patterns)} patterns to header")

def main():
    """Main function to convert binary catalogs to C headers."""
    print("Binary catalog to C header converter")
    print("====================================\n")

    # Define paths
    script_dir = Path(__file__).parent
    tracker_dir = script_dir.parent / "tracker"
    esp_main_dir = script_dir / "esp" / "main"

    star_binary = tracker_dir / "tetra3_db_stars.bin"
    pattern_binary = tracker_dir / "tetra3_db_patterns.bin"

    star_header = esp_main_dir / "star_catalog.h"
    pattern_header = esp_main_dir / "pattern_catalog.h"

    # Check input files exist
    if not star_binary.exists():
        print(f"ERROR: Star binary file not found: {star_binary}")
        return 1

    if not pattern_binary.exists():
        print(f"ERROR: Pattern binary file not found: {pattern_binary}")
        return 1

    # Create output directory if needed
    esp_main_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Read binary catalogs
        stars = read_binary_stars(star_binary)
        patterns = read_binary_patterns(pattern_binary)

        # Calculate memory usage
        star_memory = len(stars) * 24
        pattern_memory = len(patterns) * 8
        total_memory = star_memory + pattern_memory

        print(f"\nMemory usage summary:")
        print(f"  Stars:    {star_memory:,} bytes ({star_memory / 1024 / 1024:.2f} MB)")
        print(f"  Patterns: {pattern_memory:,} bytes ({pattern_memory / 1024 / 1024:.2f} MB)")
        print(f"  Total:    {total_memory:,} bytes ({total_memory / 1024 / 1024:.2f} MB)")

        # Write header files
        print(f"\nGenerating header files...")
        write_star_header(stars, star_header)
        write_pattern_header(patterns, pattern_header)

        print(f"\nHeader files generated successfully!")
        print(f"  Star catalog:    {star_header}")
        print(f"  Pattern catalog: {pattern_header}")

        print(f"\nTo use in ESP32 code:")
        print(f'  #include "star_catalog.h"')
        print(f'  #include "pattern_catalog.h"')
        print(f'  star_solver_handle_t solver = star_solver_create(')
        print(f'      star_catalog, STAR_CATALOG_SIZE,')
        print(f'      pattern_catalog, PATTERN_CATALOG_SIZE);')

        return 0

    except Exception as e:
        print(f"ERROR: {e}")
        return 1

if __name__ == "__main__":
    exit(main())