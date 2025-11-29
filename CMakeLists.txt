# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
# tell CMake to pretend any cmake_minimum_required(VERSION X) with X<3.5
# is OK.  Mirrors the -DCMAKE_POLICY_VERSION_MINIMUM=3.5 hack.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE INTERNAL
    "Allow parsing of subprojects that require CMake < 3.5")

set(GLUCOSE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/glucose)
set(NAJA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja)

message(STATUS "CMake version: ${CMAKE_VERSION}")

project(kepler-formal
  VERSION 1.0.0
  HOMEPAGE_URL https://github.com/keplertech/kepler-formal
)

# Option to enable/disable the Python interface. CI disables it by default.
option(PYTHON_INTERFACE "Enable Python interface" OFF)

# If the python interface is disabled, provide a no-op INTERFACE target so
# other targets that unconditionally link PYTHON_INTERFACE still configure.
# This must be created before add_subdirectory(...) so subprojects can link to it.
if(NOT PYTHON_INTERFACE AND NOT TARGET PYTHON_INTERFACE)
  add_library(PYTHON_INTERFACE INTERFACE)
endif()

set(ARGPARSE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja/thirdparty/argparse-3.1/include)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE)

# ---------------------------------------------------------------------------
# Find dependencies (prefer system packages; fall back to vendored copies)
# ---------------------------------------------------------------------------

# Try system packages quietly first
find_package(TBB QUIET)
find_package(yaml-cpp QUIET)
find_package(spdlog QUIET)

# If system packages didn't provide modern imported targets, try to create
# lightweight imported targets from legacy variables (keeps compatibility).
# yaml-cpp compatibility (legacy variable names vary by distro/CMake mode)
if(NOT TARGET yaml-cpp::yaml-cpp AND yaml-cpp_FOUND)
  add_library(yaml-cpp::yaml-cpp UNKNOWN IMPORTED)
  if(DEFINED YAML_CPP_LIBRARIES)
    set_target_properties(yaml-cpp::yaml-cpp PROPERTIES IMPORTED_LOCATION "${YAML_CPP_LIBRARIES}")
  elseif(DEFINED yaml-cpp_LIBRARIES)
    set_target_properties(yaml-cpp::yaml-cpp PROPERTIES IMPORTED_LOCATION "${yaml-cpp_LIBRARIES}")
  endif()
  if(DEFINED YAML_CPP_INCLUDE_DIR)
    set_property(TARGET yaml-cpp::yaml-cpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${YAML_CPP_INCLUDE_DIR}")
  elseif(DEFINED yaml-cpp_INCLUDE_DIRS)
    set_property(TARGET yaml-cpp::yaml-cpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${yaml-cpp_INCLUDE_DIRS}")
  endif()
endif()

# spdlog compatibility
if(NOT TARGET spdlog::spdlog AND spdlog_FOUND)
  add_library(spdlog::spdlog UNKNOWN IMPORTED)
  if(DEFINED spdlog_LIBRARIES)
    set_target_properties(spdlog::spdlog PROPERTIES IMPORTED_LOCATION "${spdlog_LIBRARIES}")
  endif()
  if(DEFINED spdlog_INCLUDE_DIRS)
    set_property(TARGET spdlog::spdlog PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${spdlog_INCLUDE_DIRS}")
  endif()
endif()

# TBB compatibility
if(NOT TARGET TBB::tbb AND TBB_FOUND AND DEFINED TBB_LIBRARY)
  add_library(TBB::tbb UNKNOWN IMPORTED)
  set_target_properties(TBB::tbb PROPERTIES IMPORTED_LOCATION "${TBB_LIBRARY}")
  if(DEFINED TBB_INCLUDE_DIRS)
    set_property(TARGET TBB::tbb PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIRS}")
  endif()
endif()

# ---------------------------------------------------------------------------
# Fallback: if system packages not found, try vendored thirdparty directories
# (these add_subdirectory calls must happen before any subdir that links them)
# ---------------------------------------------------------------------------

# yaml-cpp: prefer system, else vendored
if(NOT TARGET yaml-cpp::yaml-cpp)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/yaml-cpp/CMakeLists.txt")
    message(STATUS "Using vendored yaml-cpp from thirdparty/yaml-cpp")
    add_subdirectory(thirdparty/yaml-cpp EXCLUDE_FROM_ALL)
  endif()
endif()

# spdlog: prefer system, else vendored
if(NOT TARGET spdlog::spdlog)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/spdlog/CMakeLists.txt")
    message(STATUS "Using vendored spdlog from thirdparty/spdlog")
    add_subdirectory(thirdparty/spdlog EXCLUDE_FROM_ALL)
  endif()
endif()

# TBB: prefer system; vendored TBB is uncommon, so require system if not present
if(NOT TARGET TBB::tbb)
  if (NOT TBB_FOUND)
    message(FATAL_ERROR "TBB not found. Install system TBB or vendor it under thirdparty/")
  endif()
endif()

# ---------------------------------------------------------------------------
# Canonical target variables (subprojects should use these names)
# ---------------------------------------------------------------------------
if (TARGET yaml-cpp::yaml-cpp)
  set(YAML_CPP_TARGET yaml-cpp::yaml-cpp)
elseif (TARGET yaml-cpp)
  set(YAML_CPP_TARGET yaml-cpp)
else()
  message(FATAL_ERROR "yaml-cpp target not available; ensure system package installed or vendored copy present")
endif()

if (TARGET spdlog::spdlog)
  set(SPDLOG_TARGET spdlog::spdlog)
elseif (TARGET spdlog)
  set(SPDLOG_TARGET spdlog)
else()
  message(FATAL_ERROR "spdlog target not available; ensure system package installed or vendored copy present")
endif()

if (TARGET TBB::tbb)
  set(TBB_TARGET TBB::tbb)
elseif (TARGET TBB)
  set(TBB_TARGET TBB)
else()
  message(FATAL_ERROR "TBB target not available; ensure system package installed")
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
  # Include TBB headers (TBB::tbb provides interface dirs; keep fallback path for macOS Homebrew)
  target_include_directories(formal_structures
    PRIVATE
      $<$<TARGET_EXISTS:${TBB_TARGET}>:$<TARGET_PROPERTY:${TBB_TARGET},INTERFACE_INCLUDE_DIRECTORIES>>
      /opt/homebrew/include
  )

  # Link to available targets using canonical variables
  target_link_libraries(formal_structures
    PRIVATE
      $<$<BOOL:${TBB_TARGET}>:${TBB_TARGET}>
      $<$<BOOL:${YAML_CPP_TARGET}>:${YAML_CPP_TARGET}>
      $<$<BOOL:${SPDLOG_TARGET}>:${SPDLOG_TARGET}>
  )
else()
  message(WARNING "formal_structures target was not defined by subdirectories; skipping top-level link/include setup.")
endif()
