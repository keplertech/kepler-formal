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
