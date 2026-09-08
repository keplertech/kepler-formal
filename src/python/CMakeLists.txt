# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

execute_process(
  COMMAND git rev-parse --short=12 HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE KEPLER_FORMAL_GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT KEPLER_FORMAL_GIT_HASH)
  set(KEPLER_FORMAL_GIT_HASH "unknown")
endif()

Python3_add_library(kepler_formal_native MODULE WITH_SOABI PyKeplerFormal.cpp)
set_target_properties(kepler_formal_native PROPERTIES
  OUTPUT_NAME _native
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)
target_compile_definitions(kepler_formal_native PRIVATE
  KEPLER_FORMAL_VERSION="${PROJECT_VERSION}"
  KEPLER_FORMAL_GIT_HASH="${KEPLER_FORMAL_GIT_HASH}"
)
target_include_directories(kepler_formal_native PRIVATE ${PROJECT_SOURCE_DIR}/src/bin)
target_link_libraries(kepler_formal_native PRIVATE
  kepler_formal_python_driver
  Python3::Module
)

include(${CMAKE_CURRENT_SOURCE_DIR}/KeplerNestedNajaeda.cmake)
kepler_add_nested_najaeda(kepler_formal_native)

if(APPLE)
  set_target_properties(kepler_formal_native PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(kepler_formal_native PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

install(TARGETS
  kepler_formal_native
  naja_nl
  naja_dnl
  naja_bne
  naja_opt
  LIBRARY DESTINATION kepler_formal COMPONENT python
  RUNTIME DESTINATION kepler_formal COMPONENT python
  PUBLIC_HEADER DESTINATION kepler_formal/include
  COMPONENT python-headers
)

install(
  DIRECTORY "${KEPLER_NESTED_NAJAEDA_PACKAGE_DIR}/"
  DESTINATION kepler_formal/najaeda
  COMPONENT python
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
)

install(
  FILES
    "${PROJECT_SOURCE_DIR}/thirdparty/naja/LICENSE"
    "${PROJECT_SOURCE_DIR}/thirdparty/naja/AUTHORS"
  DESTINATION kepler_formal/najaeda
  COMPONENT python
)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/kepler_formal
  DESTINATION .
  COMPONENT python
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
)
