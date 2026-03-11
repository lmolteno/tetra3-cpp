cmake_minimum_required(VERSION 3.28)
project(pi_tracker)

set(CMAKE_CXX_STANDARD 20)

# Required dependencies
find_package(OpenCV REQUIRED COMPONENTS core imgcodecs imgproc)

# Sources always compiled
set(PI_SOURCES
    src/main.cpp
    src/image_loader.cpp
    src/star_detector.cpp
    src/star_map.cpp
    src/polar_align.cpp
    src/display_i2c.cpp
    src/font_6x8.cpp
)

# Optional: SDL2
find_package(SDL2 QUIET)
if(SDL2_FOUND)
    message(STATUS "SDL2 found - enabling SDL display")
    list(APPEND PI_SOURCES src/display_sdl.cpp)
else()
    message(STATUS "SDL2 not found - SDL display disabled")
endif()

# Optional: libcamera
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(LIBCAMERA QUIET libcamera)
endif()
if(LIBCAMERA_FOUND)
    message(STATUS "libcamera found - enabling camera capture")
    list(APPEND PI_SOURCES src/camera_capture.cpp)
else()
    message(STATUS "libcamera not found - camera capture disabled")
endif()

add_executable(pi_tracker ${PI_SOURCES})

target_include_directories(pi_tracker PRIVATE src)
target_link_libraries(pi_tracker PRIVATE
    tetra3_solver
    ${OpenCV_LIBS}
)

if(SDL2_FOUND)
    target_compile_definitions(pi_tracker PRIVATE HAS_SDL_DISPLAY)
    target_include_directories(pi_tracker PRIVATE ${SDL2_INCLUDE_DIRS})
    target_link_libraries(pi_tracker PRIVATE ${SDL2_LIBRARIES})
endif()

if(LIBCAMERA_FOUND)
    target_compile_definitions(pi_tracker PRIVATE HAS_LIBCAMERA)
    target_include_directories(pi_tracker PRIVATE ${LIBCAMERA_INCLUDE_DIRS})
    target_link_libraries(pi_tracker PRIVATE ${LIBCAMERA_LIBRARIES})
endif()

# Solve CLI tool (no display dependencies)
add_executable(solve_cli
    src/solve_cli.cpp
    src/image_loader.cpp
    src/star_detector.cpp
)
target_include_directories(solve_cli PRIVATE src)
target_link_libraries(solve_cli PRIVATE tetra3_solver ${OpenCV_LIBS})

# Display test executable (SDL only)
if(SDL2_FOUND)
    add_executable(display_test
        src/display_test.cpp
        src/display_sdl.cpp
        src/star_map.cpp
        src/polar_align.cpp
        src/font_6x8.cpp
    )
    target_include_directories(display_test PRIVATE src)
    target_link_libraries(display_test PRIVATE tetra3_solver ${SDL2_LIBRARIES})
    target_compile_definitions(display_test PRIVATE HAS_SDL_DISPLAY)
    target_include_directories(display_test PRIVATE ${SDL2_INCLUDE_DIRS})
endif()
