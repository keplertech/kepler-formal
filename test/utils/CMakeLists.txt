# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(SNLTruthTableTreeTests SNLTruthTableTreeTests.cpp)

target_link_libraries(SNLTruthTableTreeTests
  kepler_clauses
  naja_nl_dump
  naja_opt
  gmock gtest_main
)

GTEST_DISCOVER_TESTS(SNLTruthTableTreeTests)

# BoolExpr CNF writer tests
add_executable(BoolExprCnfWriterTests BoolExprCnfWriterTests.cpp)

target_link_libraries(BoolExprCnfWriterTests
    GTest::gtest_main
    formal_structures
    pthread
)

GTEST_DISCOVER_TESTS(BoolExprCnfWriterTests)
