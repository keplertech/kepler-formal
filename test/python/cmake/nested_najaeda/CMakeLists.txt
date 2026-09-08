# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
project(KeplerNestedNajaedaFixture LANGUAGES NONE)

include(FetchContent)
foreach(dependency IN ITEMS fmt tomlplusplus)
  string(TOUPPER "${dependency}" dependency_upper)
  set(dependency_mode "${${dependency_upper}_DEPENDENCY_MODE}")
  unset(FETCHCONTENT_SOURCE_DIR_${dependency_upper} CACHE)
  unset(${dependency}_DIR CACHE)
  set(dependency_source "${PROJECT_SOURCE_DIR}/dependency sources/${dependency}")
  if(dependency_mode STREQUAL "installed")
    set(${dependency}_DIR
      "${PROJECT_SOURCE_DIR}/installed packages/${dependency}/cmake"
      CACHE PATH "" FORCE)
  elseif(dependency_mode STREQUAL "override")
    set(FETCHCONTENT_SOURCE_DIR_${dependency_upper} "${dependency_source}"
      CACHE PATH "" FORCE)
  elseif(dependency_mode STREQUAL "declared")
    FetchContent_Declare(${dependency} SOURCE_DIR "${dependency_source}")
  elseif(dependency_mode STREQUAL "default")
    file(COPY "${dependency_source}/"
      DESTINATION "${FETCHCONTENT_BASE_DIR}/${dependency}-src")
    FetchContent_Declare(${dependency}
      SOURCE_DIR "${FETCHCONTENT_BASE_DIR}/${dependency}-src")
  else()
    message(FATAL_ERROR "Unknown ${dependency} dependency mode: ${dependency_mode}")
  endif()
  # Any accidental fallback to fetching must fail locally, without network.
  FetchContent_Declare(${dependency}
    GIT_REPOSITORY "${PROJECT_SOURCE_DIR}/missing-${dependency}-repository"
    FIND_PACKAGE_ARGS CONFIG)
  FetchContent_MakeAvailable(${dependency})
endforeach()

include("${KEPLER_SOURCE_DIR}/src/python/KeplerNestedNajaeda.cmake")
add_custom_target(kepler_formal_native)
kepler_add_nested_najaeda(kepler_formal_native)
ExternalProject_Add_StepTargets(kepler_nested_najaeda configure)
