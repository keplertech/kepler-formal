# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set(unit_designs_compare_tests
    UnitDesignCompare.cpp
)
add_executable(unitDesignCompare ${unit_designs_compare_tests})
target_compile_definitions(unitDesignCompare PRIVATE
  BENCHMARKS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/designs"
)
target_link_libraries(unitDesignCompare
  formal_strategies
  naja_snl_verilog
  gmock gtest_main
)
GTEST_DISCOVER_TESTS(unitDesignCompare DISCOVERY_TIMEOUT 30)