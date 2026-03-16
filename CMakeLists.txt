# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
# tell CMake to pretend any cmake_minimum_required(VERSION X) with X<3.5
# is OK.  Mirrors the -DCMAKE_POLICY_VERSION_MINIMUM=3.5 hack.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE INTERNAL
    "Allow parsing of subprojects that require CMake < 3.5")

set(KISSAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/kissat)
set(GLUCOSE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/glucose)
set(NAJA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja)

message(STATUS "CMake version: ${CMAKE_VERSION}")

project(kepler-formal
  VERSION 1.0.0
  HOMEPAGE_URL https://github.com/keplertech/kepler-formal
)

set(ARGPARSE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja/thirdparty/argparse-3.1/include)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE)

include(ExternalProject)


ExternalProject_Add(kissat_external
    SOURCE_DIR ${KISSAT_ROOT}
    CONFIGURE_COMMAND ./configure --ultimate
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS
        ${KISSAT_ROOT}/libkissat.a
        ${KISSAT_ROOT}/build/libkissat.a
)

add_library(kissat STATIC IMPORTED)

# Let CMake choose whichever exists
set_property(TARGET kissat PROPERTY IMPORTED_LOCATION
    ${KISSAT_ROOT}/libkissat.a
)

# Fallback for local macOS
if(NOT EXISTS "${KISSAT_ROOT}/libkissat.a")
    set_property(TARGET kissat PROPERTY IMPORTED_LOCATION
        ${KISSAT_ROOT}/build/libkissat.a
    )
endif()

add_dependencies(kissat kissat_external)

add_subdirectory(src)
add_subdirectory(thirdparty)

option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
if(ENABLE_UNIT_TESTS)
 include(CTest)
 enable_testing()
 add_subdirectory(test)
endif()

option(KEPLER_STATIC "Statically link libstdc++, libgcc, and Cap'n Proto" OFF)
if(KEPLER_STATIC)
  add_link_options(-static-libstdc++ -static-libgcc)
  # Cap'n Proto: prefer static libs when available
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  target_link_options(coverage_config INTERFACE --coverage)
endif(CODE_COVERAGE)

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
