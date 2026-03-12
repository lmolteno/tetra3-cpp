cmake_minimum_required(VERSION 3.28)
project(tetra3-cpp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Debug)

# Build the shared solver library
add_subdirectory(lib/solver)

# Original test executable
find_package(ZLIB REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(solver main.cpp)
target_link_libraries(solver PRIVATE tetra3_solver)

# Pi tracker app
add_subdirectory(pi)
