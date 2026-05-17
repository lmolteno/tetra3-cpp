cmake_minimum_required(VERSION 3.28)
project(pi_tracker)

set(CMAKE_CXX_STANDARD 20)

# pi_tracker daemon: capture + render + JSON publish.
# Solving is delegated to a solve_cli subprocess.
#
# Source tree is grouped by module — capture/, solver/, render/, io/, app/.
# Includes use module-prefixed paths (e.g. #include "render/canvas.h").
set(PI_SOURCES
    src/main.cpp
    src/capture/image_loader.cpp
    src/capture/frame_writer.cpp
    src/solver/solver_process.cpp
    src/render/star_map.cpp
    src/render/font_6x8.cpp
    src/render/text.cpp
    src/render/tracker_view.cpp
    src/io/json_sink.cpp
    src/io/control_input.cpp
    src/app/polar_align.cpp
    src/app/cam_state.cpp
    src/imu/mpu9250.cpp
    src/imu/madgwick.cpp
    src/imu/sky_math.cpp
    src/imu/imu_tracker.cpp
)

# Optional: libcamera (camera frame source). Stays dynamic — libcamera dlopens
# IPA modules at runtime so it can't be static-linked cleanly.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(LIBCAMERA QUIET libcamera)
endif()
if(LIBCAMERA_FOUND)
    message(STATUS "libcamera found - enabling camera capture")
    list(APPEND PI_SOURCES src/capture/camera_capture.cpp)
else()
    message(STATUS "libcamera not found - camera capture disabled (use --image)")
endif()

add_executable(pi_tracker ${PI_SOURCES})
target_include_directories(pi_tracker PRIVATE src)
target_link_libraries(pi_tracker PRIVATE tetra3_solver)

# Pi 2 (Cortex-A7) and up have NEON. Alpine's armhf default is vfpv3-d16 with
# no NEON; flip it on so __ARM_NEON is defined and the SIMD paths (the 2x2 bin
# in capture/, the bilinear mask in solver/star_detector) actually compile.
set(PI_NEON_ARMHF FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|armv7|armhf)" AND
   NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64")
    set(PI_NEON_ARMHF TRUE)
endif()
if(PI_NEON_ARMHF)
    target_compile_options(pi_tracker PRIVATE -mfpu=neon-vfpv4)
endif()

# Static-link the C/C++ runtimes and zlib so the only runtime dep is libcamera.
target_link_options(pi_tracker PRIVATE
    -static-libstdc++ -static-libgcc -Wl,-Bstatic -lz -Wl,-Bdynamic)

if(LIBCAMERA_FOUND)
    target_compile_definitions(pi_tracker PRIVATE HAS_LIBCAMERA)
    target_include_directories(pi_tracker PRIVATE ${LIBCAMERA_INCLUDE_DIRS})
    target_link_libraries(pi_tracker PRIVATE ${LIBCAMERA_LIBRARIES})
endif()

# oled_test: cycles through every TrackerView state, prints ASCII art to stdout,
# and optionally drives a real SSD1306 via --oled.  No camera or solver needed.
add_executable(oled_test
    src/oled_test.cpp
    src/render/tracker_view.cpp
    src/render/star_map.cpp
    src/render/font_6x8.cpp
    src/render/text.cpp
    src/app/polar_align.cpp
)
target_include_directories(oled_test PRIVATE src)
target_link_libraries(oled_test PRIVATE tetra3_solver)
if(PI_NEON_ARMHF)
    target_compile_options(oled_test PRIVATE -mfpu=neon-vfpv4)
endif()

# solve_cli: fully static. No libcamera, no OpenCV — just stb_image + zlib +
# Eigen (header-only) + the solver.
add_executable(solve_cli
    src/solver/solve_cli.cpp
    src/capture/image_loader.cpp
    src/solver/star_detector.cpp
)
target_include_directories(solve_cli PRIVATE src)
target_link_libraries(solve_cli PRIVATE tetra3_solver)
target_link_options(solve_cli PRIVATE -static)
if(PI_NEON_ARMHF)
    target_compile_options(solve_cli PRIVATE -mfpu=neon-vfpv4)
endif()
