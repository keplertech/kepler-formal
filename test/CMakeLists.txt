# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

include(GoogleTest)

# Prefer headers from the GoogleTest targets we link over system installations
# exposed by dependencies such as Boost and TBB (e.g. Homebrew's include root).
# Mixing system headers with the bundled libraries can change the GoogleTest ABI.
include_directories(BEFORE SYSTEM
  "$<TARGET_PROPERTY:gmock,INTERFACE_INCLUDE_DIRECTORIES>"
  "$<TARGET_PROPERTY:gtest,INTERFACE_INCLUDE_DIRECTORIES>"
)

add_subdirectory(unit_designs)
add_subdirectory(strategies)
add_subdirectory(utils)
add_subdirectory(scope)
add_subdirectory(sec)

# Configure-only coverage does not require building the Python extension.
add_test(NAME kepler-formal-nested-najaeda-cmake
  COMMAND ${CMAKE_COMMAND}
    "-DKEPLER_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}"
    "-DTEST_BINARY_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/nested-najaeda-cmake"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/python/test_nested_najaeda_cmake.cmake"
)

if(BUILD_KEPLER_PYTHON)
  add_subdirectory(python)
endif()
