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

add_test(NAME kepler-formal-windows-glucose-stdio
  COMMAND ${CMAKE_COMMAND}
    "-DKEPLER_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}"
    "-DTEST_BINARY_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/windows-glucose-stdio"
    "-DTEST_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
    "-DTEST_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
    "-DTEST_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/windows/test_glucose_stdio.cmake"
)
set_tests_properties(kepler-formal-windows-glucose-stdio PROPERTIES TIMEOUT 180)

if(BUILD_KEPLER_PYTHON)
  add_subdirectory(python)
endif()
