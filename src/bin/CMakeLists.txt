# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

set(kepler_version_header "${CMAKE_CURRENT_BINARY_DIR}/KeplerVersion.h")
set(KEPLER_GIT_HASH "unknown" CACHE STRING "Kepler Git hash for source archives")
add_custom_target(kepler_version_header
  COMMAND ${CMAKE_COMMAND}
    "-DKEPLER_GIT_HASH=${KEPLER_GIT_HASH}"
    "-DKEPLER_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DKEPLER_VERSION_HEADER=${kepler_version_header}"
    -P "${PROJECT_SOURCE_DIR}/cmake/KeplerVersion.cmake"
  BYPRODUCTS "${kepler_version_header}"
  VERBATIM
)
add_library(kepler_version INTERFACE)
add_dependencies(kepler_version kepler_version_header)
target_include_directories(kepler_version INTERFACE "${CMAKE_CURRENT_BINARY_DIR}")

# Result formatting and reports are shared without pulling in file loaders.
add_library(kepler_verification_support STATIC RunResult.cpp KeplerFormalUtils.cpp)
target_include_directories(kepler_verification_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_verification_support PUBLIC kepler_config kepler_sec)

if(BUILD_KEPLER_PYTHON)
  # Python borrows designs loaded by NajaEDA. Do not compile the file driver
  # or introduce its parser dependencies into the extension.
  return()
endif()

set(KEPLER_FORMAL_DRIVER_LIBRARIES
    kepler_verification_support
    kepler_version
    kepler_config
    kepler_sec
    naja_dnl
    naja_nl_dump
    naja_opt
    naja_snl_liberty
    naja_snl_systemverilog
    naja_snl_verilog
    kissat
    cadical
    formal_strategies
    scope_extraction
    yaml-cpp
    ${SPDLOG_TARGET}
)

function(configure_kepler_formal_driver target)
  target_include_directories(${target} SYSTEM BEFORE PUBLIC ${Boost_INCLUDE_DIR})
  target_include_directories(${target} PUBLIC ${ARGPARSE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${target} PUBLIC ${KEPLER_FORMAL_DRIVER_LIBRARIES})
endfunction()

add_library(kepler_formal_core STATIC KeplerFormal.cpp Btor2ExportConfig.cpp)
configure_kepler_formal_driver(kepler_formal_core)

add_library(kepler_formal_driver STATIC CppDriver.cpp)
target_link_libraries(kepler_formal_driver PUBLIC kepler_formal_core naja_snl_pyloader)

add_executable(kepler-formal KeplerFormalMain.cpp)
target_link_libraries(kepler-formal PRIVATE kepler_formal_driver)

add_dependencies(kepler-formal naja)
add_custom_command(TARGET kepler-formal POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:naja>
    $<TARGET_FILE_DIR:kepler-formal>
)

install(TARGETS kepler-formal DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES $<TARGET_FILE:naja> TYPE BIN)
