# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(secStrategyTests SequentialEquivalenceStrategyTests.cpp)

target_link_libraries(secStrategyTests
    kepler_config
    kepler_sec
    formal_strategies
    naja_nl
    ${GTEST_LIBRARIES}
    gmock
    gtest_main
)

GTEST_DISCOVER_TESTS(secStrategyTests)
