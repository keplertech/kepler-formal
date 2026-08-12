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

# This suite exercises the installed-style CLI path as part of one miter test.
# Give CTest the binary produced by this build tree instead of guessing a
# hard-coded build directory name.
add_dependencies(miterTests kepler-formal)
GTEST_DISCOVER_TESTS(
    miterTests
    PROPERTIES ENVIRONMENT "KEPLER_BIN=$<TARGET_FILE:kepler-formal>"
)

add_executable(keplerFormalCliTests
    KeplerFormalCliTests.cpp
    ${CMAKE_SOURCE_DIR}/src/bin/KeplerFormal.cpp
)

target_compile_definitions(keplerFormalCliTests PRIVATE KEPLER_FORMAL_NO_MAIN)
target_include_directories(keplerFormalCliTests PRIVATE ${CMAKE_SOURCE_DIR}/src/bin)

target_link_libraries(keplerFormalCliTests
    PRIVATE
    kepler_config
    kepler_sec
    naja_snl_pyloader
    naja_dnl
    naja_opt
    kissat
    formal_strategies
    scope_extraction
    yaml-cpp
    ${SPDLOG_TARGET}
    kepler_c2rtl_formal
    kepler_xls_c2rtl
    pthread
)

if (TARGET GTest::gtest_main)
    target_link_libraries(keplerFormalCliTests PRIVATE GTest::gtest_main)
elseif (DEFINED GTEST_MAIN_LIBRARIES)
    target_link_libraries(keplerFormalCliTests PRIVATE ${GTEST_MAIN_LIBRARIES})
else()
    target_link_libraries(keplerFormalCliTests PRIVATE ${GTEST_LIBRARIES})
endif()

GTEST_DISCOVER_TESTS(keplerFormalCliTests)

# TinyRocket is an intentionally large end-to-end SEC test.  Remove CTest's
# global timeout for this test alone so slower runners can reach its verdict.
set(tinyrocket_test_properties
    "${CMAKE_CURRENT_BINARY_DIR}/TinyRocketTestProperties.cmake")
file(GENERATE OUTPUT "${tinyrocket_test_properties}" CONTENT [=[
set_tests_properties(
    KeplerFormalCliTests.ConfigTinyRocketSecVerificationAccepted
    PROPERTIES TIMEOUT 0)
]=])
set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES
             "${tinyrocket_test_properties}")
