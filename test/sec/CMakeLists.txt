# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(secStrategyTests SequentialEquivalenceStrategyTests.cpp)
set_target_properties(secStrategyTests PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
# Keep the very large SEC test translation unit out of LTO, but do not disable
# LTO for the final executable link: CI builds GoogleTest with global -flto, and
# the linker must keep the LTO plugin enabled to resolve those archive objects.
target_compile_options(secStrategyTests PRIVATE -fno-lto)

target_link_libraries(secStrategyTests
    kepler_config
    kepler_sec
    formal_strategies
    naja_nl
    naja_snl_liberty
    naja_snl_systemverilog
)

if (TARGET GTest::gmock_main)
    target_link_libraries(secStrategyTests GTest::gmock_main)
elseif (TARGET GTest::gtest_main)
    target_link_libraries(secStrategyTests gmock GTest::gtest_main)
elseif (DEFINED GTEST_MAIN_LIBRARIES)
    target_link_libraries(secStrategyTests gmock ${GTEST_MAIN_LIBRARIES})
else()
    target_link_libraries(secStrategyTests gmock gtest_main ${GTEST_LIBRARIES})
endif()

GTEST_DISCOVER_TESTS(secStrategyTests)

add_test(
    NAME SecRegressHelperTests.PartialProofExpectationModes
    COMMAND bash
            ${CMAKE_CURRENT_SOURCE_DIR}/RunSecStrategiesRegressTests.sh
            ${PROJECT_SOURCE_DIR}/regress/run_sec_strategies_regress.sh
)
