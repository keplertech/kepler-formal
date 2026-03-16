# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(scopeExtractionTests ScopeExtractionTests.cpp)

target_link_libraries(scopeExtractionTests
  PRIVATE
    scope_extraction
    gtest_main
)

GTEST_DISCOVER_TESTS(scopeExtractionTests)
