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

# ── Application: navmesh_interactive  (SDL2 interactive viewer) ───────────────
#
# Requires SDL2 and SDL2_ttf. Install on macOS M1:
#   brew install sdl2 sdl2_ttf
#
find_package(SDL2 QUIET)

if(SDL2_FOUND OR TARGET SDL2::SDL2)

    # --- SDL2_ttf: try CMake config first, then pkg-config -------------------
    find_package(SDL2_ttf QUIET)
    if(NOT SDL2_ttf_FOUND AND NOT TARGET SDL2_ttf::SDL2_ttf)
        find_package(PkgConfig QUIET)
        if(PkgConfig_FOUND)
            pkg_check_modules(SDL2_TTF IMPORTED_TARGET SDL2_ttf)
        endif()
    endif()

    add_executable(navmesh_interactive
            app/interactive_main.cc
    )
    target_include_directories(navmesh_interactive PRIVATE app lib)
    target_include_directories(navmesh_interactive PRIVATE /opt/homebrew/include)
    target_link_libraries(navmesh_interactive PRIVATE navmesh_algo SDL2::SDL2)

    # Link SDL2_ttf via whichever discovery method worked
    if(TARGET SDL2_ttf::SDL2_ttf)
        target_link_libraries(navmesh_interactive PRIVATE SDL2_ttf::SDL2_ttf)
        message(STATUS "SDL2_ttf found via CMake config")
    elseif(TARGET PkgConfig::SDL2_TTF)
        target_link_libraries(navmesh_interactive PRIVATE PkgConfig::SDL2_TTF)
        message(STATUS "SDL2_ttf found via pkg-config")
    else()
        message(WARNING
            "SDL2_ttf not found — text labels will be absent. "
            "Install with: brew install sdl2_ttf")
        target_compile_definitions(navmesh_interactive PRIVATE SDL_NO_TTF_AVAILABLE)
    endif()

    # macOS needs several system frameworks that SDL2 depends on
    if(APPLE)
        target_link_libraries(navmesh_interactive PRIVATE
            "-framework Cocoa"
            "-framework IOKit"
            "-framework CoreVideo"
            "-framework CoreAudio"
            "-framework AudioToolbox"
        )
    endif()

    message(STATUS "navmesh_interactive target enabled (SDL2 found)")

else()
    message(WARNING
        "SDL2 not found — navmesh_interactive will NOT be built.\n"
        "To enable: brew install sdl2 sdl2_ttf")
endif()
