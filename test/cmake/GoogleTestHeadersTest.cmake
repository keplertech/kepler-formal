# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)

foreach(required_variable IN ITEMS KEPLER_SOURCE_DIR TEST_BINARY_DIR TEST_CXX_COMPILER)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "${required_variable} must be provided")
  endif()
endforeach()

set(fixture_source "${TEST_BINARY_DIR}/source")
set(fixture_build "${TEST_BINARY_DIR}/build")
file(MAKE_DIRECTORY
  "${fixture_source}/fake-sdk/include/gtest"
  "${fixture_source}/fake-sdk/include/gmock")

# A dependency's broad SDK include directory must never replace the headers
# belonging to the GoogleTest libraries that the executable links.
file(WRITE "${fixture_source}/fake-sdk/include/gtest/gtest.h"
  "#error Wrong GoogleTest headers selected from the SDK include directory\n")
file(WRITE "${fixture_source}/fake-sdk/include/gmock/gmock.h"
  "#error Wrong GoogleMock headers selected from the SDK include directory\n")
file(WRITE "${fixture_source}/GoogleTestHeaders.cpp" [=[
#include <gtest/gtest.h>
#include <gmock/gmock.h>

TEST(GoogleTestHeaders, UsesLinkedLibraryHeaders) {
  EXPECT_THAT(42, ::testing::Eq(42));
}
]=])

string(CONFIGURE [=[
cmake_minimum_required(VERSION 3.21)
project(KeplerGoogleTestHeaders LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK ON CACHE BOOL "" FORCE)
set(gtest_build_tests OFF CACHE BOOL "" FORCE)
set(gmock_build_tests OFF CACHE BOOL "" FORCE)
add_subdirectory("@KEPLER_SOURCE_DIR@/thirdparty/naja/thirdparty/googletest"
                 vendor EXCLUDE_FROM_ALL)

add_library(fake_sdk INTERFACE)
target_include_directories(fake_sdk SYSTEM INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/fake-sdk/include")

include("@KEPLER_SOURCE_DIR@/cmake/KeplerGoogleTest.cmake")
add_executable(googleTestHeaders GoogleTestHeaders.cpp)
target_link_libraries(googleTestHeaders PRIVATE fake_sdk GTest::gmock_main)
]=] fixture_project @ONLY)
file(WRITE "${fixture_source}/CMakeLists.txt" "${fixture_project}")

execute_process(
  COMMAND "${CMAKE_COMMAND}" -S "${fixture_source}" -B "${fixture_build}"
          "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}"
          -DCMAKE_BUILD_TYPE=Debug
  RESULT_VARIABLE configure_result
  OUTPUT_VARIABLE configure_output
  ERROR_VARIABLE configure_error
  TIMEOUT 60)
if(NOT configure_result STREQUAL "0")
  message(FATAL_ERROR
    "GoogleTest header fixture configuration failed (${configure_result}):\n"
    "${configure_output}${configure_error}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${fixture_build}"
          --config Debug --target googleTestHeaders --parallel 2
  RESULT_VARIABLE build_result
  OUTPUT_VARIABLE build_output
  ERROR_VARIABLE build_error
  TIMEOUT 120)
if(NOT build_result STREQUAL "0")
  message(FATAL_ERROR
    "GoogleTest header fixture build failed (${build_result}):\n"
    "${build_output}${build_error}")
endif()

message(STATUS "GoogleTest and GoogleMock headers match the linked libraries")
