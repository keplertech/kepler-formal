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
    naja_nl
    ${GTEST_LIBRARIES}
    naja_opt
    naja_snl_visual
    kepler_clauses          # 👈 place first
    formal_strategies       # 👈 after its dependencies
    formal_structures
    naja_bne
    naja_dnl
    pthread
    naja_snl_pyloader
)

GTEST_DISCOVER_TESTS(SNLTruthTableTreeTests)