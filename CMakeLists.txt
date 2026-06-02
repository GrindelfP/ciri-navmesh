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

# ── Приложение (заглушка) ─────────────────────────────────────────────────────
# add_executable(navmesh_app app/main.cpp)
# target_link_libraries(navmesh_app PRIVATE geometry)