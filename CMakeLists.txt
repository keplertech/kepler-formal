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

# Find packages early so we can create compatibility imported targets if needed.
find_package(TBB QUIET)
find_package(yaml-cpp QUIET)
find_package(spdlog QUIET)

# If find_package didn't create modern imported targets but did set legacy variables,
# create lightweight imported targets for consistent usage below.
# yaml-cpp compatibility
if(NOT TARGET yaml-cpp::yaml-cpp)
  if(yaml-cpp_FOUND)
    # If module-mode provided YAML_CPP_LIBRARIES / INCLUDE dirs, create an imported target
    add_library(yaml-cpp::yaml-cpp UNKNOWN IMPORTED)
    if(DEFINED YAML_CPP_LIBRARIES)
      set_target_properties(yaml-cpp::yaml-cpp PROPERTIES
        IMPORTED_LOCATION "${YAML_CPP_LIBRARIES}")
    elseif(DEFINED yaml-cpp_LIBRARIES)
      set_target_properties(yaml-cpp::yaml-cpp PROPERTIES
        IMPORTED_LOCATION "${yaml-cpp_LIBRARIES}")
    endif()
    if(DEFINED YAML_CPP_INCLUDE_DIR)
      set_property(TARGET yaml-cpp::yaml-cpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${YAML_CPP_INCLUDE_DIR}")
    elseif(DEFINED yaml-cpp_INCLUDE_DIRS)
      set_property(TARGET yaml-cpp::yaml-cpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${yaml-cpp_INCLUDE_DIRS}")
    endif()
  else()
    message(STATUS "yaml-cpp not found via quiet search; will fallback to REQUIRED search momentarily")
  endif()
endif()

# --- after your find_package(...) and compatibility imported-target creation ---

# Canonical target variables (prefer modern names, fall back to vendored unnamespaced targets)
if (TARGET yaml-cpp::yaml-cpp)
  set(YAML_CPP_TARGET yaml-cpp::yaml-cpp)
elseif (TARGET yaml-cpp)
  set(YAML_CPP_TARGET yaml-cpp)
else()
  # yaml-cpp should already be REQUIRED above; this is defensive
  message(FATAL_ERROR "yaml-cpp target not available")
endif()

if (TARGET spdlog::spdlog)
  set(SPDLOG_TARGET spdlog::spdlog)
elseif (TARGET spdlog)
  set(SPDLOG_TARGET spdlog)
else()
  message(FATAL_ERROR "spdlog target not available")
endif()

if (TARGET TBB::tbb)
  set(TBB_TARGET TBB::tbb)
elseif (TARGET TBB)
  set(TBB_TARGET TBB)
else()
  message(FATAL_ERROR "TBB target not available")
endif()

# If any package is required by your project policy, explicitly require them now.
# Using REQUIRED to cause clear configure failure if not present.
if(NOT TARGET yaml-cpp::yaml-cpp)
  find_package(yaml-cpp REQUIRED)  # will error if missing
endif()
if(NOT TARGET spdlog::spdlog)
  find_package(spdlog REQUIRED)
endif()
if(NOT TARGET TBB::tbb)
  find_package(TBB REQUIRED)
endif()

# Add project subdirectories (these define formal_structures and others).
add_subdirectory(src)
add_subdirectory(thirdparty)
include(CTest)
add_subdirectory(test)

# Now, only set include_directories / link libraries if the formal_structures
# target exists. This avoids ordering issues.
if(TARGET formal_structures)
  # Include TBB headers (TBB::tbb provides interface dirs; keep fallback path for macOS Homebrew)
  target_include_directories(formal_structures
    PRIVATE
      $<$<TARGET_EXISTS:TBB::tbb>:$<TARGET_PROPERTY:TBB::tbb,INTERFACE_INCLUDE_DIRECTORIES>>
      /opt/homebrew/include
  )

  # Link to available targets in a guarded way.
  target_link_libraries(formal_structures
    PRIVATE
      $<$<TARGET_EXISTS:TBB::tbb>:TBB::tbb>
      $<$<TARGET_EXISTS:yaml-cpp::yaml-cpp>:yaml-cpp::yaml-cpp>
      $<$<TARGET_EXISTS:spdlog::spdlog>:spdlog::spdlog>
  )
else()
  message(WARNING "formal_structures target was not defined by subdirectories; skipping top-level link/include setup.")
endif()