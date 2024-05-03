cmake_minimum_required(VERSION 3.9)
project(FrameFieldTracing)

set(CMAKE_CXX_STANDARD 17)

set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" CACHE PATH "Data path for various input files.")

# Dependencies
if (NOT TARGET glfw)
    add_subdirectory(libs/glfw)
endif()

if (NOT TARGET clean-core)
    add_subdirectory(libs/clean-core)
endif()

if (NOT TARGET typed-geometry)
    add_subdirectory(libs/typed-geometry)
endif()

if (NOT TARGET polymesh)
    add_subdirectory(libs/polymesh)
endif()

set(GLOW_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}) # Viewer fonts will be placed here
if (NOT TARGET glow)
    add_subdirectory(libs/glow)
endif()

if (NOT TARGET imgui)
    add_subdirectory(libs/imgui-lean)
endif()

if (NOT TARGET glow-extras)
    add_subdirectory(libs/glow-extras)
endif()

if (NOT TARGET OpenMeshCore)
    add_subdirectory(libs/OpenMesh)
endif()

# Find Eigen if not already present
if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 REQUIRED)
endif()

# FrameFieldTracing Library (library directory)
file(GLOB_RECURSE LIBRARY_SOURCE_FILES "src/FrameFieldTracing/*.cc"
    "src/FrameFieldTracing/*.hh"
    "src/FrameFieldTracing/*.c"
    "src/FrameFieldTracing/*.h")
add_library(FrameFieldTracing ${LIBRARY_SOURCE_FILES})
target_link_libraries(FrameFieldTracing PUBLIC typed-geometry polymesh glow-extras OpenMeshCore Eigen3::Eigen)
target_include_directories(FrameFieldTracing PUBLIC src)
target_compile_definitions(FrameFieldTracing PUBLIC DATA_PATH="${DATA_PATH}")
target_link_libraries(FrameFieldTracing PUBLIC stdc++fs)

add_executable("FrameFieldTracer" "src/main/main.cc")
target_link_libraries("FrameFieldTracer" PUBLIC FrameFieldTracing)
