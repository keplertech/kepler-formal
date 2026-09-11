# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
project(kepler_windows_glucose_stdio LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)

foreach(solver IN ITEMS KISSAT CADICAL GLUCOSE)
  string(TOLOWER "${solver}" directory)
  set(${solver}_ROOT "${KEPLER_SOURCE_DIR}/thirdparty/${directory}")
endforeach()

# Exercise the real Windows CMake branch without making a POSIX compiler
# pretend to provide the Windows SDK. Only the Glucose target is built.
set(host_is_windows "${WIN32}")
if(host_is_windows)
  add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
set(WIN32 TRUE)
include("${KEPLER_SOURCE_DIR}/src/windows/KeplerWindowsSolvers.cmake")
set(WIN32 "${host_is_windows}")

get_target_property(definitions glucose COMPILE_DEFINITIONS)
if(NOT "putc_unlocked=putc" IN_LIST definitions)
  message(FATAL_ERROR "Windows Glucose must replace putc_unlocked with putc")
endif()
get_target_property(public_definitions glucose INTERFACE_COMPILE_DEFINITIONS)
if("putc_unlocked=putc" IN_LIST public_definitions)
  message(FATAL_ERROR "The stdio replacement must not leak to Glucose consumers")
endif()

if(NOT host_is_windows)
  # /FI is clang-cl syntax; native POSIX headers already provide gettimeofday.
  set_property(TARGET glucose PROPERTY COMPILE_OPTIONS "")
  set_property(TARGET glucose PROPERTY INTERFACE_COMPILE_OPTIONS "")
endif()

add_executable(glucose_stdio_smoke glucose_stdio_smoke.cpp)
target_link_libraries(glucose_stdio_smoke PRIVATE glucose)
enable_testing()
add_test(NAME all-byte-values COMMAND glucose_stdio_smoke)
