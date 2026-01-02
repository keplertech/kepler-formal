# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

include(GoogleTest)

cmake_minimum_required(VERSION 3.10)
project(SNLTruthTableTreeTests)

# Enable testing
enable_testing()

# Google Test setup
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Add main and test files
add_executable(SNLTruthTableTreeTests SNLTruthTableTreeTests.cpp)

target_link_libraries(SNLTruthTableTreeTests
    ${GTEST_LIBRARIES}
    kepler_clauses          # 👈 place first
    formal_strategies       # 👈 after its dependencies
    formal_structures
    pthread
)

GTEST_DISCOVER_TESTS(SNLTruthTableTreeTests)