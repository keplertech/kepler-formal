# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
# tell CMake to pretend any cmake_minimum_required(VERSION X) with X<3.5
# is OK.  Mirrors the -DCMAKE_POLICY_VERSION_MINIMUM=3.5 hack.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE INTERNAL
    "Allow parsing of subprojects that require CMake < 3.5")

set(KISSAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/kissat)
set(GLUCOSE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/glucose)
set(CADICAL_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/cadical)
set(NAJA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja)

message(STATUS "CMake version: ${CMAKE_VERSION}")

project(kepler-formal
  VERSION 1.0.0
  HOMEPAGE_URL https://github.com/keplertech/kepler-formal
)

set(ARGPARSE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/naja/thirdparty/argparse-3.1/include)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard" FORCE)

include(ExternalProject)

set(CADICAL_KITTEN_SYMBOL_PREFIX_DEFINES
    -Dcompletely_backtrack_to_root_level=cadical_completely_backtrack_to_root_level
    -Dcitten_clause_with_id=cadical_citten_clause_with_id
    -Dcitten_clause_with_id_and_equivalence=cadical_citten_clause_with_id_and_equivalence
    -Dcitten_clause_with_id_and_exception=cadical_citten_clause_with_id_and_exception
    -Dkitten_add_prime_implicant=cadical_kitten_add_prime_implicant
    -Dkitten_assume=cadical_kitten_assume
    -Dkitten_assume_signed=cadical_kitten_assume_signed
    -Dkitten_binary=cadical_kitten_binary
    -Dkitten_clause=cadical_kitten_clause
    -Dkitten_clause_with_id_and_exception=cadical_kitten_clause_with_id_and_exception
    -Dkitten_clear=cadical_kitten_clear
    -Dkitten_compute_clausal_core=cadical_kitten_compute_clausal_core
    -Dkitten_compute_prime_implicant=cadical_kitten_compute_prime_implicant
    -Dkitten_current_ticks=cadical_kitten_current_ticks
    -Dkitten_failed=cadical_kitten_failed
    -Dkitten_fixed=cadical_kitten_fixed
    -Dkitten_fixed_signed=cadical_kitten_fixed_signed
    -Dkitten_flip_and_implicant_for_signed_literal=cadical_kitten_flip_and_implicant_for_signed_literal
    -Dkitten_flip_literal=cadical_kitten_flip_literal
    -Dkitten_flip_phases=cadical_kitten_flip_phases
    -Dkitten_flip_signed_literal=cadical_kitten_flip_signed_literal
    -Dkitten_init=cadical_kitten_init
    -Dkitten_no_terminator=cadical_kitten_no_terminator
    -Dkitten_no_ticks_limit=cadical_kitten_no_ticks_limit
    -Dkitten_randomize_phases=cadical_kitten_randomize_phases
    -Dkitten_release=cadical_kitten_release
    -Dkitten_set_logging=cadical_kitten_set_logging
    -Dkitten_set_terminator=cadical_kitten_set_terminator
    -Dkitten_set_ticks_limit=cadical_kitten_set_ticks_limit
    -Dkitten_shrink_to_clausal_core=cadical_kitten_shrink_to_clausal_core
    -Dkitten_shuffle_clauses=cadical_kitten_shuffle_clauses
    -Dkitten_signed_value=cadical_kitten_signed_value
    -Dkitten_solve=cadical_kitten_solve
    -Dkitten_status=cadical_kitten_status
    -Dkitten_trace_core=cadical_kitten_trace_core
    -Dkitten_track_antecedents=cadical_kitten_track_antecedents
    -Dkitten_traverse_core_clauses=cadical_kitten_traverse_core_clauses
    -Dkitten_traverse_core_clauses_with_id=cadical_kitten_traverse_core_clauses_with_id
    -Dkitten_traverse_core_ids=cadical_kitten_traverse_core_ids
    -Dkitten_unit=cadical_kitten_unit
    -Dkitten_value=cadical_kitten_value
    -Dnew_learned_klause=cadical_new_learned_klause
)
string(JOIN " " CADICAL_KITTEN_SYMBOL_PREFIX_FLAGS
    ${CADICAL_KITTEN_SYMBOL_PREFIX_DEFINES}
)

ExternalProject_Add(kissat_external
    SOURCE_DIR ${KISSAT_ROOT}
    # Keep Kissat small and quiet for embedded use, but do not use
    # `--ultimate`: it implies `--no-options`, which compiles out runtime
    # option storage and prevents SEC from selecting a SAT profile for its
    # cone-of-influence proof obligations.
    # Honor the configured compiler (e.g. a CMAKE_TOOLCHAIN_FILE such as
    # deps/toolchain.cmake); kissat's configure defaults to host cc.
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
        "CC=${CMAKE_C_COMPILER}"
        ./configure --compact --quiet --no-proofs
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

ExternalProject_Add(cadical_external
    SOURCE_DIR ${CADICAL_ROOT}
    # CC/CXX: honor the configured compiler (e.g. a CMAKE_TOOLCHAIN_FILE
    # such as deps/toolchain.cmake); cadical's configure defaults to host g++.
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
        "CC=${CMAKE_C_COMPILER}"
        "CXX=${CMAKE_CXX_COMPILER}"
        "CFLAGS=${CADICAL_KITTEN_SYMBOL_PREFIX_FLAGS}"
        "CXXFLAGS=${CADICAL_KITTEN_SYMBOL_PREFIX_FLAGS}"
        ./configure -q --no-tracing --no-contrib
    BUILD_COMMAND make -C build libcadical.a
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS
        ${CADICAL_ROOT}/build/libcadical.a
)

add_library(cadical STATIC IMPORTED)
set_property(TARGET cadical PROPERTY IMPORTED_LOCATION
    ${CADICAL_ROOT}/build/libcadical.a
)
set_property(TARGET cadical PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    ${CADICAL_ROOT}/src
    ${CADICAL_ROOT}/contrib
)
add_dependencies(cadical cadical_external)

set(KEPLER_CPU_TARGET "portable" CACHE STRING
    "CPU instruction-set baseline (portable, native, or x86-64-v3)")
set_property(CACHE KEPLER_CPU_TARGET PROPERTY
             STRINGS portable native x86-64-v3)

set(_kepler_cpu_targets portable native x86-64-v3)
if(NOT KEPLER_CPU_TARGET IN_LIST _kepler_cpu_targets)
  message(FATAL_ERROR
    "Unknown KEPLER_CPU_TARGET='${KEPLER_CPU_TARGET}'. "
    "Expected one of: portable, native, x86-64-v3.")
endif()

function(_kepler_classify_processor processor output_variable)
  string(TOLOWER "${processor}" _kepler_processor)
  if(_kepler_processor MATCHES "^(x86_64|x86_64h|amd64|x64)$")
    set(${output_variable} "x86-64" PARENT_SCOPE)
  elseif(_kepler_processor MATCHES "^(arm64|aarch64)$")
    set(${output_variable} "arm64" PARENT_SCOPE)
  else()
    set(${output_variable} "${_kepler_processor}" PARENT_SCOPE)
  endif()
endfunction()

set(_kepler_target_processor "${CMAKE_SYSTEM_PROCESSOR}")
set(_kepler_osx_architecture_count 0)
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
  list(LENGTH CMAKE_OSX_ARCHITECTURES _kepler_osx_architecture_count)
  if(_kepler_osx_architecture_count EQUAL 1)
    list(GET CMAKE_OSX_ARCHITECTURES 0 _kepler_target_processor)
  endif()
endif()
_kepler_classify_processor("${_kepler_target_processor}"
                           _kepler_target_architecture)

set(_kepler_cpu_compile_option "")
if(KEPLER_CPU_TARGET STREQUAL "portable")
  get_cmake_property(_kepler_cmake_variables VARIABLES)
  set(_kepler_native_flag_variables "")
  foreach(_kepler_variable IN LISTS _kepler_cmake_variables)
    if(_kepler_variable MATCHES "^CMAKE_CXX_FLAGS($|_.+)" AND
       NOT _kepler_variable MATCHES "_INIT$")
      string(FIND "${${_kepler_variable}}" "-march=native"
             _kepler_march_native_position)
      string(FIND "${${_kepler_variable}}" "-mcpu=native"
             _kepler_mcpu_native_position)
      if(NOT _kepler_march_native_position EQUAL -1 OR
         NOT _kepler_mcpu_native_position EQUAL -1)
        list(APPEND _kepler_native_flag_variables "${_kepler_variable}")
      endif()
    endif()
  endforeach()

  if(_kepler_native_flag_variables)
    list(REMOVE_DUPLICATES _kepler_native_flag_variables)
    list(SORT _kepler_native_flag_variables)
    list(JOIN _kepler_native_flag_variables ", "
         _kepler_native_flag_variable_list)
    message(WARNING
      "KEPLER_CPU_TARGET=portable adds no CPU-specific compiler flags, but "
      "the following C++ flag variables contain -march=native or "
      "-mcpu=native: ${_kepler_native_flag_variable_list}. This contradicts "
      "the portable baseline and can produce binaries that fail on other "
      "CPUs. Remove the global native flag instead of relying on "
      "KEPLER_CPU_TARGET=portable to override it.")
  endif()
elseif(KEPLER_CPU_TARGET STREQUAL "native")
  if(CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=native is not supported while cross-compiling. "
      "Use KEPLER_CPU_TARGET=portable or an explicit target baseline.")
  endif()
  if(APPLE AND _kepler_osx_architecture_count GREATER 1)
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=native is not valid for a macOS universal build "
      "(${CMAKE_OSX_ARCHITECTURES}). Use KEPLER_CPU_TARGET=portable.")
  endif()

  _kepler_classify_processor("${CMAKE_HOST_SYSTEM_PROCESSOR}"
                             _kepler_host_architecture)
  if(APPLE AND CMAKE_OSX_ARCHITECTURES AND
     NOT _kepler_target_architecture STREQUAL _kepler_host_architecture)
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=native cannot target '${_kepler_target_processor}' "
      "from host processor '${CMAKE_HOST_SYSTEM_PROCESSOR}'. Native tuning "
      "is only valid when building for the current machine.")
  endif()

  if(_kepler_target_architecture STREQUAL "x86-64")
    set(_kepler_cpu_compile_option "-march=native")
  elseif(_kepler_target_architecture STREQUAL "arm64")
    set(_kepler_cpu_compile_option "-mcpu=native")
  else()
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=native supports x86-64 and arm64/aarch64 targets; "
      "detected target processor '${_kepler_target_processor}'.")
  endif()
elseif(KEPLER_CPU_TARGET STREQUAL "x86-64-v3")
  if(APPLE AND _kepler_osx_architecture_count GREATER 1)
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=x86-64-v3 is not valid for a macOS universal build "
      "(${CMAKE_OSX_ARCHITECTURES}). Use KEPLER_CPU_TARGET=portable.")
  endif()
  if(NOT _kepler_target_architecture STREQUAL "x86-64")
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=x86-64-v3 requires an x86-64 target, but the target "
      "processor is '${_kepler_target_processor}'. It is not valid on Apple "
      "Silicon or other arm64/aarch64 targets.")
  endif()
  set(_kepler_cpu_compile_option "-march=x86-64-v3")
endif()

if(_kepler_cpu_compile_option)
  include(CheckCXXCompilerFlag)
  string(MAKE_C_IDENTIFIER "${_kepler_cpu_compile_option}"
         _kepler_cpu_compile_option_identifier)
  string(TOUPPER "${_kepler_cpu_compile_option_identifier}"
         _kepler_cpu_compile_option_identifier)
  set(_kepler_cpu_flag_check
      "KEPLER_CXX_SUPPORTS_${_kepler_cpu_compile_option_identifier}")
  check_cxx_compiler_flag("${_kepler_cpu_compile_option}"
                          ${_kepler_cpu_flag_check})
  if(NOT "${${_kepler_cpu_flag_check}}")
    message(FATAL_ERROR
      "KEPLER_CPU_TARGET=${KEPLER_CPU_TARGET} requires "
      "'${_kepler_cpu_compile_option}', but ${CMAKE_CXX_COMPILER_ID} "
      "${CMAKE_CXX_COMPILER_VERSION} does not support that option.")
  endif()
  # Intentionally inherit the selected baseline in src/ and thirdparty/naja/.
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:${_kepler_cpu_compile_option}>")
endif()

message(STATUS "Kepler CPU target: ${KEPLER_CPU_TARGET}")

add_subdirectory(src)
add_subdirectory(thirdparty)

if(TARGET slang_slang)
  # fmt 12 no longer exposes fmt::format through <fmt/core.h>. Keep vendored
  # slang sources untouched and provide the full formatting declarations here.
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(slang_slang PRIVATE -include fmt/format.h)
  elseif(MSVC)
    target_compile_options(slang_slang PRIVATE /FIfmt/format.h)
  endif()
endif()

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
