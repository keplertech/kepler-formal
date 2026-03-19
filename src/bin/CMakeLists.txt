# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(kepler-formal KeplerFormal.cpp)
target_include_directories(kepler-formal SYSTEM BEFORE PUBLIC ${Boost_INCLUDE_DIR})
target_include_directories(kepler-formal PUBLIC ${ARGPARSE_DIR})

target_link_libraries(kepler-formal
  PRIVATE
    kepler_config
    naja_snl_pyloader
    naja_dnl
    naja_opt
    kissat
    formal_strategies
    scope_extraction
    yaml-cpp
    ${SPDLOG_TARGET}
)

install(TARGETS kepler-formal DESTINATION ${CMAKE_INSTALL_BINDIR})
