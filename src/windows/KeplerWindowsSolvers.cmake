# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# Build only the embedded libraries: the upstream shell configure scripts and
# stand-alone solver applications require POSIX. Keep the submodules pristine.
set(kepler_solver_compat "${CMAKE_CURRENT_LIST_DIR}/compat")
set(kepler_solver_generated "${CMAKE_CURRENT_BINARY_DIR}/windows-solvers")

file(STRINGS "${KISSAT_ROOT}/VERSION" KEPLER_KISSAT_VERSION LIMIT_COUNT 1)
configure_file("${CMAKE_CURRENT_LIST_DIR}/kissat-build.h.in"
  "${kepler_solver_generated}/build.h" @ONLY)

file(GLOB kissat_sources CONFIGURE_DEPENDS "${KISSAT_ROOT}/src/*.c")
foreach(application IN ITEMS main application handle parse witness)
  list(REMOVE_ITEM kissat_sources "${KISSAT_ROOT}/src/${application}.c")
endforeach()
add_library(kissat STATIC ${kissat_sources})
set_target_properties(kissat PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON)
target_include_directories(kissat PRIVATE
  "${kepler_solver_compat}" "${kepler_solver_generated}")
# Match the existing --compact --quiet --no-proofs build; do not use NOPTIONS,
# since Kepler's SEC strategies must still select solver profiles at runtime.
target_compile_definitions(kissat PRIVATE COMPACT QUIET NPROOFS NDEBUG)

file(GLOB cadical_sources CONFIGURE_DEPENDS
  "${CADICAL_ROOT}/src/*.cpp" "${CADICAL_ROOT}/src/*.c")
list(REMOVE_ITEM cadical_sources
  "${CADICAL_ROOT}/src/cadical.cpp" "${CADICAL_ROOT}/src/mobical.cpp")
add_library(cadical STATIC ${cadical_sources})
set_target_properties(cadical PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON)
target_include_directories(cadical
  PRIVATE "${kepler_solver_compat}"
  PUBLIC "${CADICAL_ROOT}/src" "${CADICAL_ROOT}/contrib")
# NBUILD is upstream's supported alternative to its generated build.hpp.
# Windows CRT has no getc_unlocked/putc_unlocked; use the supported fallback.
target_compile_definitions(cadical PRIVATE NBUILD QUIET NTRACING NUNLOCKED NDEBUG)
foreach(prefix_definition IN LISTS CADICAL_KITTEN_SYMBOL_PREFIX_DEFINES)
  string(REGEX REPLACE "^-D" "" prefix_definition "${prefix_definition}")
  target_compile_definitions(cadical PRIVATE "${prefix_definition}")
endforeach()
if(WIN32)
  target_link_libraries(cadical PRIVATE psapi)
endif()

find_package(ZLIB REQUIRED)
add_library(glucose STATIC
  "${GLUCOSE_ROOT}/core/Solver.cc"
  "${GLUCOSE_ROOT}/core/lcm.cc"
  "${GLUCOSE_ROOT}/simp/SimpSolver.cc"
  "${GLUCOSE_ROOT}/utils/Options.cc"
  "${GLUCOSE_ROOT}/utils/System.cc")
target_include_directories(glucose PUBLIC
  "${kepler_solver_compat}" "${GLUCOSE_ROOT}")
target_link_libraries(glucose PUBLIC ZLIB::ZLIB)
if(WIN32)
  # Glucose's public System.h uses timeval/gettimeofday even in its Windows
  # branch, without including sys/time.h. Consumers instantiate those inlines.
  target_compile_options(glucose PUBLIC "/FI${kepler_solver_compat}/sys/time.h")
endif()
