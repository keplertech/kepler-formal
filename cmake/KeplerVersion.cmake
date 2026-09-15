# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# Run at build time so switching commits also refreshes an existing build.
if(NOT DEFINED KEPLER_GIT_HASH OR KEPLER_GIT_HASH STREQUAL "")
  set(KEPLER_GIT_HASH "unknown")
endif()

find_package(Git QUIET)
# Archives nested inside another checkout must not inherit its Git hash.
if(GIT_FOUND AND EXISTS "${KEPLER_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${KEPLER_SOURCE_DIR}"
    RESULT_VARIABLE git_result
    OUTPUT_VARIABLE git_hash
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(git_result EQUAL 0 AND NOT git_hash STREQUAL "")
    set(KEPLER_GIT_HASH "${git_hash}")
  endif()
endif()

configure_file("${KEPLER_SOURCE_DIR}/src/bin/KeplerVersion.h.in"
               "${KEPLER_VERSION_HEADER}" @ONLY)
