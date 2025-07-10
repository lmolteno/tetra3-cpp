cmake_minimum_required(VERSION 3.31)
project(solver)

set(CMAKE_CXX_STANDARD 20)

find_package(ZLIB REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(solver main.cpp)

target_link_libraries(solver PRIVATE ZLIB::ZLIB)
target_link_libraries(solver PRIVATE Eigen3::Eigen)
