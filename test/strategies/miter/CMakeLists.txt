# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

include(GoogleTest)

cmake_minimum_required(VERSION 3.10)
project(MiterTests)

# Enable testing
enable_testing()

# Google Test setup
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Add main and test files
add_executable(miterTests MiterTests.cpp)

target_link_libraries(miterTests
    naja_nl
    ${GTEST_LIBRARIES}
    naja_opt
    naja_snl_visual
    kepler_clauses          # 👈 place first
    formal_strategies       # 👈 after its dependencies
    formal_structures
    naja_dnl
    pthread
    naja_snl_pyloader
)

GTEST_DISCOVER_TESTS(miterTests)