# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE INTERNAL
    "Allow parsing of subprojects that require CMake < 3.5")

set(GLUCOSE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/glucose)
set(NAJA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja)

message(STATUS "CMake version: ${CMAKE_VERSION}")

project(kepler-formal
  VERSION 1.0.0
  HOMEPAGE_URL https://github.com/keplertech/kepler-formal
)

option(PYTHON_INTERFACE "Enable Python interface" OFF)

if(NOT PYTHON_INTERFACE AND NOT TARGET PYTHON_INTERFACE)
  add_library(PYTHON_INTERFACE INTERFACE)
endif()

set(ARGPARSE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja/thirdparty/argparse-3.1/include)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE)

# ---------------------------------------------------------------------------
# Try system packages first (quiet). We'll fall back to vendored copies below.
# ---------------------------------------------------------------------------
find_package(TBB QUIET)
find_package(spdlog QUIET)

# ---------------------------------------------------------------------------
# If system packages provided legacy variables but no imported targets, create
# lightweight imported targets so subprojects can link consistently.
# ---------------------------------------------------------------------------

if(NOT TARGET spdlog::spdlog AND spdlog_FOUND)
  add_library(spdlog::spdlog UNKNOWN IMPORTED)
  if(DEFINED spdlog_LIBRARIES)
    set_target_properties(spdlog::spdlog PROPERTIES IMPORTED_LOCATION "${spdlog_LIBRARIES}")
  endif()
  if(DEFINED spdlog_INCLUDE_DIRS)
    set_property(TARGET spdlog::spdlog PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${spdlog_INCLUDE_DIRS}")
  endif()
endif()

if(NOT TARGET TBB::tbb AND TBB_FOUND AND DEFINED TBB_LIBRARY)
  add_library(TBB::tbb UNKNOWN IMPORTED)
  set_target_properties(TBB::tbb PROPERTIES IMPORTED_LOCATION "${TBB_LIBRARY}")
  if(DEFINED TBB_INCLUDE_DIRS)
    set_property(TARGET TBB::tbb PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIRS}")
  endif()
endif()

# ---------------------------------------------------------------------------
# Vendored fallback: add_subdirectory for thirdparty copies if system not found.
# These must be added before subdirs that link to them.
# ---------------------------------------------------------------------------

if(NOT TARGET spdlog::spdlog)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/spdlog/CMakeLists.txt")
    message(STATUS "Using vendored spdlog from thirdparty/spdlog")
    add_subdirectory(thirdparty/spdlog EXCLUDE_FROM_ALL)
  endif()
endif()

# Require TBB from system (common on CI). If you vendor TBB, add similar fallback.
if(NOT TARGET TBB::tbb)
  if(NOT TBB_FOUND)
    message(FATAL_ERROR "TBB not found. Install system TBB or vendor it under thirdparty/")
  endif()
endif()

# ---------------------------------------------------------------------------
# Canonical target variables (subprojects should use these names)
# ---------------------------------------------------------------------------

if (TARGET spdlog::spdlog)
  set(SPDLOG_TARGET spdlog::spdlog)
elseif (TARGET spdlog)
  set(SPDLOG_TARGET spdlog)
else()
  message(FATAL_ERROR "spdlog target not available; install libspdlog-dev or add thirdparty/spdlog")
endif()

if (TARGET TBB::tbb)
  set(TBB_TARGET TBB::tbb)
elseif (TARGET TBB)
  set(TBB_TARGET TBB)
else()
  message(FATAL_ERROR "TBB target not available; ensure system TBB is installed")
endif()

message(STATUS "Using YAML_CPP_TARGET = ${YAML_CPP_TARGET}")
message(STATUS "Using SPDLOG_TARGET  = ${SPDLOG_TARGET}")
message(STATUS "Using TBB_TARGET     = ${TBB_TARGET}")

# ---------------------------------------------------------------------------
# Add project subdirectories (these define formal_structures and others).
# ---------------------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(thirdparty)
include(CTest)
add_subdirectory(test)

# ---------------------------------------------------------------------------
# Top-level link/include setup for formal_structures (guarded)
# ---------------------------------------------------------------------------
if(TARGET formal_structures)
  target_include_directories(formal_structures
    PRIVATE
      $<$<TARGET_EXISTS:${TBB_TARGET}>:$<TARGET_PROPERTY:${TBB_TARGET},INTERFACE_INCLUDE_DIRECTORIES>>
      /opt/homebrew/include
  )

  target_link_libraries(formal_structures
    PRIVATE
      $<$<BOOL:${TBB_TARGET}>:${TBB_TARGET}>
      $<$<BOOL:${YAML_CPP_TARGET}>:${YAML_CPP_TARGET}>
      $<$<BOOL:${SPDLOG_TARGET}>:${SPDLOG_TARGET}>
  )
else()
  message(WARNING "formal_structures target was not defined by subdirectories; skipping top-level link/include setup.")
endif()
