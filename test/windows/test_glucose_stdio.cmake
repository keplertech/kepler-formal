# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
foreach(required_variable IN ITEMS KEPLER_SOURCE_DIR TEST_BINARY_DIR)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "This test requires ${required_variable}")
  endif()
endforeach()

set(configure_options "")
foreach(language IN ITEMS C CXX)
  if(TEST_${language}_COMPILER)
    list(APPEND configure_options
      "-DCMAKE_${language}_COMPILER:FILEPATH=${TEST_${language}_COMPILER}")
  endif()
endforeach()
if(TEST_TOOLCHAIN_FILE)
  list(APPEND configure_options
    "-DCMAKE_TOOLCHAIN_FILE:FILEPATH=${TEST_TOOLCHAIN_FILE}")
endif()

function(run_checked)
  execute_process(COMMAND ${ARGV}
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Glucose stdio regression failed:\n${output}\n${error}")
  endif()
endfunction()

run_checked("${CMAKE_COMMAND}"
  -S "${CMAKE_CURRENT_LIST_DIR}/glucose_stdio" -B "${TEST_BINARY_DIR}"
  "-DKEPLER_SOURCE_DIR:PATH=${KEPLER_SOURCE_DIR}"
  -DCMAKE_BUILD_TYPE:STRING=Release ${configure_options})
run_checked("${CMAKE_COMMAND}" --build "${TEST_BINARY_DIR}"
  --config Release --target glucose_stdio_smoke --parallel 2)
run_checked("${CMAKE_CTEST_COMMAND}" --test-dir "${TEST_BINARY_DIR}"
  -C Release --output-on-failure)
message(STATUS "Windows Glucose private stdio replacement and all 256 bytes passed")
