cmake_minimum_required(VERSION 3.28)
project(tetra3_solver)

set(CMAKE_CXX_STANDARD 20)

find_package(ZLIB REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(tetra3_solver STATIC
    src/solver.cpp
)

target_include_directories(tetra3_solver PUBLIC include)

target_link_libraries(tetra3_solver PUBLIC
    ZLIB::ZLIB
    Eigen3::Eigen
)
