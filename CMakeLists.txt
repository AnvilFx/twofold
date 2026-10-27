cmake_minimum_required(VERSION 3.21)
project(twofold VERSION 0.1.0 LANGUAGES C CXX)

option(TWOFOLD_BUILD_TESTS "Build the test suite" ${PROJECT_IS_TOP_LEVEL})
option(TWOFOLD_BUILD_VIEWER "Build the interactive viewer, which downloads raylib" ${PROJECT_IS_TOP_LEVEL})

add_library(twofold
    src/collide2d.cpp
    src/collide3d.cpp
    src/distance.cpp
    src/hull.cpp
    src/outline.cpp
    src/polyhedron.cpp
    src/query.cpp
    src/shape.cpp
    src/solver.cpp
    src/world.cpp
)
add_library(twofold::twofold ALIAS twofold)

target_include_directories(twofold PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(twofold PUBLIC cxx_std_17)

# Identical results across compilers and platforms need strict float semantics.
if(MSVC)
    target_compile_options(twofold PRIVATE /W4 /fp:precise)
else()
    target_compile_options(twofold PRIVATE -Wall -Wextra -ffp-contract=off)
endif()

if(TWOFOLD_BUILD_TESTS)
    enable_testing()
    add_executable(twofold_tests tests/main.cpp tests/test_hull.cpp tests/test_outline.cpp tests/test_world.cpp)
    target_link_libraries(twofold_tests PRIVATE twofold)
    add_test(NAME twofold_tests COMMAND twofold_tests)
endif()

if(TWOFOLD_BUILD_VIEWER)
    add_subdirectory(viewer)
endif()
