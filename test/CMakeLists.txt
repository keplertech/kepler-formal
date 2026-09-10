# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

include(GoogleTest)
include(${PROJECT_SOURCE_DIR}/cmake/KeplerGoogleTest.cmake)

add_test(NAME GoogleTestHeaderPrecedence
  COMMAND ${CMAKE_COMMAND}
    -DKEPLER_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DTEST_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/cmake/googletest_headers
    -DTEST_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GoogleTestHeadersTest.cmake
)

add_subdirectory(unit_designs)
add_subdirectory(strategies)
add_subdirectory(utils)
add_subdirectory(scope)
add_subdirectory(sec)

if(BUILD_KEPLER_PYTHON)
  add_subdirectory(python)
endif()
