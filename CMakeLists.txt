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

set(ARGPARSE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja/thirdparty/argparse-3.1/include)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE)

add_subdirectory(src)
add_subdirectory(thirdparty)

option(ENABLE_UNIT_TESTS ON)
if(ENABLE_UNIT_TESTS)
  include(CTest)
  add_subdirectory(test)
endif()

find_package(TBB REQUIRED)

# if you’re using find_package(TBB), ensure it picked up Homebrew’s install

target_include_directories(formal_structures
  PRIVATE
    ${TBB_INCLUDE_DIRS}          # if defined by find_package
    /opt/homebrew/include        # fallback for tbb headers
)
target_link_libraries(formal_structures
  PRIVATE TBB::tbb
)
