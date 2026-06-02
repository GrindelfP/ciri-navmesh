cmake_minimum_required(VERSION 3.16)
project(navmesh_lib LANGUAGES CXX)

# ── Стандарт ─────────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Предупреждения ────────────────────────────────────────────────────────────
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror
            -Wno-unused-parameter)
endif()

# ── Библиотека: geometry ──────────────────────────────────────────────────────
add_library(geometry STATIC
        lib/geometry/primitives.cpp
        lib/geometry/predicates.cpp
        lib/geometry/dcel.cc
)
target_include_directories(geometry PUBLIC lib)

# ── Тесты ────────────────────────────────────────────────────────────────────
enable_testing()

add_executable(test_geometry tests/test_geometry.cpp)
target_link_libraries(test_geometry PRIVATE geometry)
add_test(NAME geometry_tests COMMAND test_geometry)

# ── Приложение: navmesh_bench ─────────────────────────────────────────────────
add_executable(navmesh_bench
        app/main.cc
        lib/navmesh/astar.cc
        lib/navmesh/funnel.cc
        lib/navmesh/navmesh.cc
        lib/triangulation/delaunay.cc
        lib/triangulation/greedy.cc
        lib/triangulation/quasi_greedy.cc
)

# Добавляем пути к заголовочным файлам, чтобы CMake их видел
target_include_directories(navmesh_bench PRIVATE app lib)

# Линкуем базовую статическую библиотеку геометрии
target_link_libraries(navmesh_bench PRIVATE geometry)
