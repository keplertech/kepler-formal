# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.21)
foreach(required_variable IN ITEMS KEPLER_SOURCE_DIR TEST_BINARY_DIR)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "This test requires ${required_variable}")
  endif()
endforeach()

set(fixture_source "${TEST_BINARY_DIR}/source with spaces")
file(MAKE_DIRECTORY "${fixture_source}")
file(COPY "${CMAKE_CURRENT_LIST_DIR}/cmake/nested_najaeda/"
  DESTINATION "${fixture_source}")
foreach(dependency IN ITEMS fmt tomlplusplus)
  set(dependency_source "${fixture_source}/dependency sources/${dependency}")
  set(dependency_package "${fixture_source}/installed packages/${dependency}/cmake")
  file(MAKE_DIRECTORY "${dependency_source}" "${dependency_package}")
  file(WRITE "${dependency_source}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 3.21)\nproject(${dependency} LANGUAGES NONE)\n")
  file(WRITE "${dependency_package}/${dependency}-config.cmake"
    "add_library(${dependency}::${dependency} INTERFACE IMPORTED)\n")
endforeach()

function(check_dependency_modes build_name fmt_mode tomlplusplus_mode)
  set(build_dir "${TEST_BINARY_DIR}/${build_name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${fixture_source}" -B "${build_dir}"
      "-DKEPLER_SOURCE_DIR:PATH=${KEPLER_SOURCE_DIR}"
      "-DFMT_DEPENDENCY_MODE:STRING=${fmt_mode}"
      "-DTOMLPLUSPLUS_DEPENDENCY_MODE:STRING=${tomlplusplus_mode}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Parent ${build_name} configure failed:\n${output}\n${error}")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${build_dir}"
      --target kepler_nested_najaeda-configure
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Nested ${build_name} configure failed:\n${output}\n${error}")
  endif()
  foreach(dependency IN ITEMS fmt tomlplusplus)
    set(mode "${${dependency}_mode}")
    file(READ "${build_dir}/kepler_nested_najaeda-build/${dependency}-source.txt"
      actual_source)
    if(mode STREQUAL "installed")
      set(expected_source "")
    elseif(mode STREQUAL "default")
      set(expected_source "${build_dir}/_deps/${dependency}-src")
    else()
      set(expected_source "${fixture_source}/dependency sources/${dependency}")
    endif()
    if(NOT actual_source STREQUAL expected_source)
      message(FATAL_ERROR
        "${mode}: ${dependency} used '${actual_source}', expected '${expected_source}'")
    endif()
  endforeach()
  message(STATUS
    "Nested NajaEDA dependency configure passed: ${build_name} (${fmt_mode}, ${tomlplusplus_mode})")
endfunction()

check_dependency_modes("default build" default default)
check_dependency_modes("declared build" declared declared)
check_dependency_modes("override build" override override)
check_dependency_modes("installed build" installed installed)
check_dependency_modes("mixed source first build" override installed)
check_dependency_modes("mixed package first build" installed override)
# Reuse the external build cache while changing how dependencies are resolved.
check_dependency_modes("override build" installed installed)
check_dependency_modes("override build" override override)
