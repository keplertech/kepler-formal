# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(kepler_formal KeplerFormal.cpp)
target_include_directories(kepler_formal SYSTEM BEFORE PUBLIC ${Boost_INCLUDE_DIR})
target_include_directories(kepler_formal PUBLIC ${ARGPARSE_DIR})
target_link_libraries(kepler_formal
    naja_snl_pyloader
    naja_dnl naja_opt formal_strategies $<$<TARGET_EXISTS:yaml-cpp::yaml-cpp>:yaml-cpp::yaml-cpp>
    $<$<TARGET_EXISTS:spdlog::spdlog>:spdlog::spdlog> )
install(TARGETS kepler_formal DESTINATION ${CMAKE_INSTALL_BINDIR})
