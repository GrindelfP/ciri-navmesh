cmake_minimum_required(VERSION 3.16)
project(navmesh_lib LANGUAGES CXX)

# ── Standard ──────────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Warnings ──────────────────────────────────────────────────────────────────
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror
            -Wno-unused-parameter
            -Wno-comment)
endif()

# ── Library: geometry ─────────────────────────────────────────────────────────
add_library(geometry STATIC
        lib/geometry/primitives.cpp
        lib/geometry/predicates.cpp
        lib/geometry/dcel.cc
)
target_include_directories(geometry PUBLIC lib)

# ── Library: navmesh_algo ─────────────────────────────────────────────────────
# Shared algorithm code used by both navmesh_bench and navmesh_vis.
add_library(navmesh_algo STATIC
        lib/navmesh/astar.cc
        lib/navmesh/funnel.cc
        lib/navmesh/navmesh.cc
        lib/triangulation/delaunay.cc
        lib/triangulation/greedy.cc
        lib/triangulation/quasi_greedy.cc
)
target_include_directories(navmesh_algo PUBLIC lib)
target_link_libraries(navmesh_algo PUBLIC geometry)

# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()

add_executable(test_geometry tests/test_geometry.cpp)
target_link_libraries(test_geometry PRIVATE geometry)
add_test(NAME geometry_tests COMMAND test_geometry)

# ── Application: navmesh_bench ────────────────────────────────────────────────
add_executable(navmesh_bench
        app/main.cc
)
target_include_directories(navmesh_bench PRIVATE app lib)
target_link_libraries(navmesh_bench PRIVATE navmesh_algo)

# ── Application: navmesh_vis  (SVG visualiser, no external deps) ──────────────
add_executable(navmesh_vis
        app/vis_main.cc
        app/svg_writer.cc
)
target_include_directories(navmesh_vis PRIVATE app lib)
target_link_libraries(navmesh_vis PRIVATE navmesh_algo)
