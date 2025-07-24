cmake_minimum_required(VERSION 3.10)
project(HelloTest)

# Enable testing
enable_testing()

# Google Test setup
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Add main and test files
add_executable(runTests MiterTests.cpp)

target_link_libraries(runTests
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
)




# Add test to CTest
add_test(NAME HelloWorldTest COMMAND runTests)
