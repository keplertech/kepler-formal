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
      ${ARGN}
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

check_dependency_modes("platform build" override override
  "-DTEST_PLATFORM_FORWARDING:BOOL=ON")
set(platform_build "${TEST_BINARY_DIR}/platform build")
file(READ "${platform_build}/kepler_nested_najaeda-build/CMakeCache.txt" nested_cache)
foreach(expected IN ITEMS
    "CMAKE_TOOLCHAIN_FILE:FILEPATH=${fixture_source}/test toolchain.cmake"
    "CMAKE_PREFIX_PATH:STRING=${fixture_source}/prefix one\;${fixture_source}/prefix two"
    "CMAKE_MSVC_RUNTIME_LIBRARY:STRING=MultiThreadedDLL"
    "VCPKG_TARGET_TRIPLET:STRING=x64-windows"
    "VCPKG_HOST_TRIPLET:STRING=x64-windows"
    "Python3_FIND_ABI:STRING=ANY\;ANY\;ANY\;ON"
    "PREGENERATED_PARSER_SOURCES:BOOL=ON")
  string(REPLACE "\\;" ";" expected "${expected}")
  string(FIND "${nested_cache}" "${expected}\n" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Nested platform configuration lost: ${expected}")
  endif()
endforeach()
file(GLOB_RECURSE platform_build_files
  "${platform_build}/*install*.cmake" "${platform_build}/*build.make"
  "${platform_build}/*.ninja")
set(found_pyd FALSE)
foreach(build_file IN LISTS platform_build_files)
  file(READ "${build_file}" contents)
  string(FIND "${contents}" "-DMODULE_SUFFIX:STRING=.pyd" position)
  if(NOT position EQUAL -1)
    set(found_pyd TRUE)
  endif()
endforeach()
if(NOT found_pyd)
  message(FATAL_ERROR "Windows nested staging must look for .pyd, not .dll modules")
endif()
message(STATUS "Nested NajaEDA platform/ABI configuration forwarding passed")

# Exercise package staging with Windows file names on every host.
set(stage_source "${TEST_BINARY_DIR}/windows stage source")
set(stage_destination "${TEST_BINARY_DIR}/windows staged package")
file(MAKE_DIRECTORY "${stage_source}")
file(WRITE "${stage_source}/__init__.py" "from najaeda import netlist\n")
file(WRITE "${stage_source}/netlist.py" "from najaeda.naja import SNLUniverse\n")
file(WRITE "${stage_source}/naja.cp314-win_amd64.pyd" "fixture module\n")
file(WRITE "${stage_source}/libkepler_najaeda_naja_nl.dll" "fixture runtime\n")
execute_process(COMMAND "${CMAKE_COMMAND}"
  "-DSOURCE_PACKAGE_DIR:PATH=${stage_source}"
  "-DDESTINATION_PACKAGE_DIR:PATH=${stage_destination}"
  "-DSHARED_LIBRARY_PREFIX:STRING=libkepler_najaeda_"
  "-DSHARED_LIBRARY_SUFFIX:STRING=.dll"
  "-DMODULE_SUFFIX:STRING=.pyd"
  "-DSTAMP_FILE:FILEPATH=${TEST_BINARY_DIR}/windows stage.stamp"
  -P "${KEPLER_SOURCE_DIR}/src/python/StageKeplerNestedNajaeda.cmake"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Windows package staging failed:\n${output}\n${error}")
endif()
file(READ "${stage_destination}/__init__.py" staged_init)
file(READ "${stage_destination}/netlist.py" staged_netlist)
if(NOT staged_init STREQUAL "from kepler_formal.najaeda import netlist\n"
    OR NOT staged_netlist STREQUAL "from kepler_formal.najaeda.naja import SNLUniverse\n"
    OR NOT EXISTS "${stage_destination}/naja.cp314-win_amd64.pyd"
    OR NOT EXISTS "${stage_destination}/libkepler_najaeda_naja_nl.dll")
  message(FATAL_ERROR "Windows package staging lost modules or nested imports")
endif()
message(STATUS "Nested NajaEDA Windows package staging passed")
check_dependency_modes("installed build" installed installed)
check_dependency_modes("mixed source first build" override installed)
check_dependency_modes("mixed package first build" installed override)
# Reuse the external build cache while changing how dependencies are resolved.
check_dependency_modes("override build" installed installed)
check_dependency_modes("override build" override override)
