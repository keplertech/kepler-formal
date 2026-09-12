# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# Use headers from the GoogleTest targets we link, even when a dependency adds
# an SDK prefix containing a different GoogleTest installation. Keep these
# directories ahead of other system includes for all Kepler test targets.
include_directories(BEFORE SYSTEM
  "$<TARGET_PROPERTY:GTest::gtest,INTERFACE_INCLUDE_DIRECTORIES>"
  "$<TARGET_PROPERTY:GTest::gmock,INTERFACE_INCLUDE_DIRECTORIES>"
)
