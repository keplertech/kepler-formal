# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set(KEPLER_FORMAL_DRIVER_LIBRARIES
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

# Both adapters call the same core. Only the standalone adapter links the
# embedded Python primitive loader; the Python adapter runs in its host runtime.
add_library(kepler_formal_python_driver STATIC ../python/PythonDriver.cpp)
target_link_libraries(kepler_formal_python_driver PUBLIC kepler_formal_core)
target_include_directories(kepler_formal_python_driver PUBLIC
  ${PROJECT_SOURCE_DIR}/src/python)

if(BUILD_KEPLER_PYTHON)
  # Build the standalone executable in a separate, ordinary CMake build. It
  # must not acquire a runtime dependency on an installed NajaEDA wheel.
  return()
endif()

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
