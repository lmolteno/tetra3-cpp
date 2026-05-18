cmake_minimum_required(VERSION 3.28)
project(tetra3-cpp)

set(CMAKE_CXX_STANDARD 20)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Prefer static zlib so the Pi binaries can be statically linked against it.
# Falls through to .so if the static build isn't installed.
set(ZLIB_USE_STATIC_LIBS ON)

# Plate-solver library — used by every target.
add_subdirectory(lib/solver)

# Pi tracker daemon + solve_cli + oled_test + tests.
add_subdirectory(pi)
