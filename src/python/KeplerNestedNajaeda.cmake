# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

include_guard(GLOBAL)
include(ExternalProject)

function(kepler_add_nested_najaeda consumer_target)
  if(NOT TARGET "${consumer_target}")
    message(FATAL_ERROR
      "kepler_add_nested_najaeda expected target ${consumer_target}")
  endif()

  set(external_target kepler_nested_najaeda)
  set(external_source_dir "${PROJECT_SOURCE_DIR}/thirdparty/naja")
  set(external_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/kepler_nested_najaeda-build")
  set(external_install_dir "${CMAKE_CURRENT_BINARY_DIR}/kepler_nested_najaeda-install")
  set(package_dir "${CMAKE_CURRENT_BINARY_DIR}/kepler_nested_najaeda-package")
  set(package_stamp "${CMAKE_CURRENT_BINARY_DIR}/kepler_nested_najaeda-package.stamp")
  set(project_include "${CMAKE_CURRENT_LIST_DIR}/KeplerNestedNajaedaProject.cmake")
  set(stage_script "${CMAKE_CURRENT_LIST_DIR}/StageKeplerNestedNajaeda.cmake")
  set(isolated_library_prefix "libkepler_najaeda_")

  set(fmt_source_dir "${FETCHCONTENT_BASE_DIR}/fmt-src")
  set(tomlplusplus_source_dir "${FETCHCONTENT_BASE_DIR}/tomlplusplus-src")
  foreach(fetched_source IN ITEMS "${fmt_source_dir}" "${tomlplusplus_source_dir}")
    if(NOT IS_DIRECTORY "${fetched_source}")
      message(FATAL_ERROR
        "The nested NajaEDA build requires the parent-fetched source directory "
        "${fetched_source}")
    endif()
  endforeach()

  set(external_cmake_args
    "-DCMAKE_BUILD_TYPE:STRING=Release"
    "-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON"
    "-DCMAKE_POLICY_VERSION_MINIMUM:STRING=3.5"
    "-DCMAKE_INSTALL_PREFIX:PATH=${external_install_dir}"
    "-DCMAKE_PROJECT_INCLUDE:FILEPATH=${project_include}"
    "-DKEPLER_NESTED_NAJAEDA_BUILD:BOOL=ON"
    "-DBUILD_NAJA_PYTHON:BOOL=ON"
    "-DPREGENERATED_PARSER_SOURCES:BOOL=OFF"
    "-DBUILD_BENCHMARKS:BOOL=OFF"
    "-DCODE_COVERAGE:BOOL=OFF"
    "-DENABLE_SANITIZERS:BOOL=OFF"
    "-DENABLE_THREAD_SANITIZER:BOOL=OFF"
    "-DPython3_EXECUTABLE:FILEPATH=${Python3_EXECUTABLE}"
    "-DFETCHCONTENT_SOURCE_DIR_FMT:PATH=${fmt_source_dir}"
    "-DFETCHCONTENT_SOURCE_DIR_TOMLPLUSPLUS:PATH=${tomlplusplus_source_dir}"
  )

  foreach(compiler_variable IN ITEMS CMAKE_C_COMPILER CMAKE_CXX_COMPILER)
    if(DEFINED ${compiler_variable} AND NOT "${${compiler_variable}}" STREQUAL "")
      list(APPEND external_cmake_args
        "-D${compiler_variable}:FILEPATH=${${compiler_variable}}")
    endif()
  endforeach()

  foreach(parser_tool_variable IN ITEMS BISON_EXECUTABLE FLEX_EXECUTABLE)
    if(DEFINED ${parser_tool_variable}
        AND NOT "${${parser_tool_variable}}" STREQUAL "")
      list(APPEND external_cmake_args
        "-D${parser_tool_variable}:FILEPATH=${${parser_tool_variable}}")
    endif()
  endforeach()

  foreach(osx_variable IN ITEMS
      CMAKE_OSX_ARCHITECTURES
      CMAKE_OSX_DEPLOYMENT_TARGET
      CMAKE_OSX_SYSROOT)
    if(DEFINED ${osx_variable} AND NOT "${${osx_variable}}" STREQUAL "")
      string(REPLACE ";" "|" osx_value "${${osx_variable}}")
      list(APPEND external_cmake_args "-D${osx_variable}:STRING=${osx_value}")
    endif()
  endforeach()

  set(install_byproducts)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.26)
    list(APPEND install_byproducts INSTALL_BYPRODUCTS "${package_stamp}")
  endif()

  ExternalProject_Add(${external_target}
    SOURCE_DIR "${external_source_dir}"
    BINARY_DIR "${external_binary_dir}"
    INSTALL_DIR "${external_install_dir}"
    DOWNLOAD_COMMAND ""
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    LIST_SEPARATOR "|"
    CMAKE_ARGS ${external_cmake_args}
    BUILD_COMMAND
      "${CMAKE_COMMAND}" --build <BINARY_DIR> --config Release
    INSTALL_COMMAND
      "${CMAKE_COMMAND}" --install <BINARY_DIR> --config Release
        --prefix <INSTALL_DIR>
      COMMAND
      "${CMAKE_COMMAND}"
        "-DSOURCE_PACKAGE_DIR:PATH=<INSTALL_DIR>/najaeda"
        "-DDESTINATION_PACKAGE_DIR:PATH=${package_dir}"
        "-DSHARED_LIBRARY_PREFIX:STRING=${isolated_library_prefix}"
        "-DSHARED_LIBRARY_SUFFIX:STRING=${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "-DMODULE_SUFFIX:STRING=${CMAKE_SHARED_MODULE_SUFFIX}"
        "-DSTAMP_FILE:FILEPATH=${package_stamp}"
        -P "${stage_script}"
    BUILD_ALWAYS TRUE
    EXCLUDE_FROM_ALL TRUE
    USES_TERMINAL_BUILD TRUE
    USES_TERMINAL_INSTALL TRUE
    ${install_byproducts}
  )

  add_dependencies("${consumer_target}" ${external_target})

  set(KEPLER_NESTED_NAJAEDA_TARGET "${external_target}"
    CACHE INTERNAL "Target that stages the isolated nested NajaEDA runtime" FORCE)
  set(KEPLER_NESTED_NAJAEDA_PACKAGE_DIR "${package_dir}"
    CACHE INTERNAL "Staged kepler_formal.najaeda package directory" FORCE)
  set(KEPLER_NESTED_NAJAEDA_PACKAGE_STAMP "${package_stamp}"
    CACHE INTERNAL "Stamp produced after staging the nested NajaEDA package" FORCE)

  set(KEPLER_NESTED_NAJAEDA_TARGET "${external_target}" PARENT_SCOPE)
  set(KEPLER_NESTED_NAJAEDA_PACKAGE_DIR "${package_dir}" PARENT_SCOPE)
  set(KEPLER_NESTED_NAJAEDA_PACKAGE_STAMP "${package_stamp}" PARENT_SCOPE)
endfunction()
