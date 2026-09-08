# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_subdirectory(naja)
if(NOT WIN32)
  # Windows builds the serial Glucose library through KeplerWindowsSolvers;
  # upstream's CMake also builds POSIX-only command-line/parallel targets.
  add_subdirectory(glucose)
endif()
add_subdirectory(yaml-cpp)
#add_subdirectory(kissat)
# cadical is built through the root ExternalProject because upstream does not
# provide a CMake subdirectory target.
#add_subdirectory(cpptrace EXCLUDE_FROM_ALL)
