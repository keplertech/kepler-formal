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
  target_compile_definitions(${target} PRIVATE KEPLER_FORMAL_NO_MAIN)
  target_include_directories(${target} SYSTEM BEFORE PUBLIC ${Boost_INCLUDE_DIR})
  target_include_directories(${target} PUBLIC ${ARGPARSE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${target} PUBLIC ${KEPLER_FORMAL_DRIVER_LIBRARIES})
endfunction()

add_library(kepler_formal_driver STATIC KeplerFormal.cpp)
configure_kepler_formal_driver(kepler_formal_driver)
target_link_libraries(kepler_formal_driver PUBLIC naja_snl_pyloader)

if(BUILD_KEPLER_PYTHON)
  # The in-process API rejects Python technology files. Compile that loader
  # path out so a wheel does not embed Naja's separate Python runtime.
  add_library(kepler_formal_python_driver STATIC KeplerFormal.cpp)
  configure_kepler_formal_driver(kepler_formal_python_driver)
  target_compile_definitions(
    kepler_formal_python_driver
    PRIVATE KEPLER_FORMAL_NO_PY_TECH
  )
endif()

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
