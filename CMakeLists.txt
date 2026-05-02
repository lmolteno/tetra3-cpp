cmake_minimum_required(VERSION 3.28)
project(tetra3-cpp)

set(CMAKE_CXX_STANDARD 20)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Prefer the static zlib so binaries can be statically linked against it.
# Falls through to .so if the static build isn't installed.
set(ZLIB_USE_STATIC_LIBS ON)

# Build the shared solver library
add_subdirectory(lib/solver)

# Original test executable
find_package(ZLIB REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(solver main.cpp)
target_link_libraries(solver PRIVATE tetra3_solver)

# Pi tracker app
add_subdirectory(pi)
