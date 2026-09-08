# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
project(naja LANGUAGES NONE)

if(NOT KEPLER_NESTED_NAJAEDA_BUILD
    OR NOT CMAKE_SHARED_LIBRARY_PREFIX STREQUAL "libkepler_najaeda_")
  message(FATAL_ERROR "The nested build lost its isolated library configuration")
endif()

include(FetchContent)
foreach(dependency IN ITEMS fmt tomlplusplus)
  FetchContent_Declare(${dependency}
    GIT_REPOSITORY "${PROJECT_SOURCE_DIR}/missing-${dependency}-repository"
    FIND_PACKAGE_ARGS CONFIG)
  FetchContent_MakeAvailable(${dependency})
  set(dependency_source "")
  FetchContent_GetProperties(${dependency} SOURCE_DIR dependency_source)
  if(NOT dependency_source AND NOT TARGET ${dependency}::${dependency})
    message(FATAL_ERROR "Neither sources nor an installed ${dependency} were found")
  endif()
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${dependency}-source.txt"
    "${dependency_source}")
endforeach()
