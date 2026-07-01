# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_executable(miterTests MiterTests.cpp)

target_link_libraries(miterTests
    kepler_config
    naja_nl
    ${GTEST_LIBRARIES}
    naja_opt
    naja_snl_pyloader
    formal_strategies
    gmock gtest_main
)

GTEST_DISCOVER_TESTS(miterTests)

add_executable(keplerFormalCliTests
    KeplerFormalCliTests.cpp
    ${CMAKE_SOURCE_DIR}/src/bin/KeplerFormal.cpp
)

target_compile_definitions(keplerFormalCliTests PRIVATE KEPLER_FORMAL_NO_MAIN)
target_include_directories(keplerFormalCliTests PRIVATE ${CMAKE_SOURCE_DIR}/src/bin)

target_link_libraries(keplerFormalCliTests
    kepler_config
    kepler_c_frontend
    kepler_sec
    naja_snl_pyloader
    naja_dnl
    naja_opt
    kissat
    formal_strategies
    scope_extraction
    yaml-cpp
    ${SPDLOG_TARGET}
    pthread
)

if (TARGET GTest::gtest_main)
    target_link_libraries(keplerFormalCliTests GTest::gtest_main)
elseif (DEFINED GTEST_MAIN_LIBRARIES)
    target_link_libraries(keplerFormalCliTests ${GTEST_MAIN_LIBRARIES})
else()
    target_link_libraries(keplerFormalCliTests ${GTEST_LIBRARIES})
endif()

GTEST_DISCOVER_TESTS(keplerFormalCliTests)

find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
  set(KEPLER_VERILOG_C_BENCHMARK_CASES
    safe/Arbiter
    safe/Barrel
    safe/BufAl
    safe/CRC
    safe/FIFO
    safe/FourbyFour
    safe/RCU
    safe/Rrobin
    safe/Spinner
    safe/TicTacToe
    safe/TreeArb
    safe/Vending
    safe/Vlunc
    safe/Vsa16
    safe/b01
    safe/b02
    safe/b03
    safe/b06
    safe/b07
    safe/b08
    safe/traffic_light_example
    unsafe/DAIO
    unsafe/Dekker
    unsafe/FIFO
    unsafe/Ibuf
    unsafe/NullModem
    unsafe/ReqAck
    unsafe/Rotate
    unsafe/Unidec
    unsafe/Vlunc
    unsafe/b01
    unsafe/non-pipelined-microprocessor
  )

  foreach(_verilog_c_case IN LISTS KEPLER_VERILOG_C_BENCHMARK_CASES)
    string(REPLACE "/" "." _verilog_c_test_suffix "${_verilog_c_case}")
    add_test(
      NAME KeplerFormalCliBenchmarks.VerilogC.${_verilog_c_test_suffix}
      COMMAND ${Python3_EXECUTABLE}
              ${CMAKE_CURRENT_SOURCE_DIR}/RunVerilogCBenchmark.py
              --kepler-formal $<TARGET_FILE:kepler-formal>
              --benchmarks-root ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/verilog-c
              --case ${_verilog_c_case}
              --timeout 120
    )
    set_tests_properties(
      KeplerFormalCliBenchmarks.VerilogC.${_verilog_c_test_suffix}
      PROPERTIES
        LABELS "cli;c2rtl;benchmark;verilog-c"
        SKIP_RETURN_CODE 77
        TIMEOUT 180
    )
  endforeach()

  add_test(
    NAME KeplerFormalCliBenchmarks.VerilogC.NonEquivalent.safe.Spinner.reference_mutated
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/RunVerilogCBenchmark.py
            --kepler-formal $<TARGET_FILE:kepler-formal>
            --benchmarks-root ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/verilog-c
            --case safe/Spinner
            --reference-mutation spinner_reference_identity
            --expect different
            --timeout 120
  )
  set_tests_properties(
    KeplerFormalCliBenchmarks.VerilogC.NonEquivalent.safe.Spinner.reference_mutated
    PROPERTIES
      LABELS "cli;c2rtl;benchmark;verilog-c;verilog-c-negative"
      SKIP_RETURN_CODE 77
      TIMEOUT 180
  )
endif()
