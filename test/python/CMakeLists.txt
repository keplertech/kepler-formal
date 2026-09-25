# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

add_test(NAME kepler-formal-python-tests-setup
  COMMAND ${CMAKE_COMMAND}
    -DSRC=${PROJECT_SOURCE_DIR}/src/python/kepler_formal
    -DDST=${CMAKE_CURRENT_BINARY_DIR}/kepler_formal
    -DNATIVE=$<TARGET_FILE:kepler_formal_native>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/setup_test_env.cmake
)
set_tests_properties(kepler-formal-python-tests-setup PROPERTIES
  FIXTURES_SETUP kepler-formal-python-test-env
)

add_test(NAME kepler-formal-python-tests
  COMMAND ${Python3_EXECUTABLE} -m unittest discover -v ${CMAKE_CURRENT_SOURCE_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(kepler-formal-python-tests PROPERTIES
  FIXTURES_REQUIRED kepler-formal-python-test-env
  ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable(kepler-borrowed-native-tests borrowed_native_tests.cpp)
# The shared provider includes its Python bridge. A standalone test executable
# supplies Python symbols normally provided by the host interpreter when the
# borrowed-design verifier is loaded as an extension.
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module Development.Embed)
target_link_libraries(kepler-borrowed-native-tests PRIVATE
  kepler_borrowed_designs Python3::Python)
najaeda_fixup_consumer(kepler-borrowed-native-tests)
add_test(NAME kepler-formal-borrowed-native-tests
  COMMAND kepler-borrowed-native-tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(kepler-formal-borrowed-native-tests PROPERTIES TIMEOUT 120)

add_executable(kepler-borrowed-policy-tests borrowed_policy_tests.cpp)
target_link_libraries(kepler-borrowed-policy-tests PRIVATE kepler_borrowed_designs Python3::Python)
najaeda_fixup_consumer(kepler-borrowed-policy-tests)
add_test(NAME kepler-formal-borrowed-policy-tests COMMAND kepler-borrowed-policy-tests)
set_tests_properties(kepler-formal-borrowed-policy-tests PROPERTIES TIMEOUT 120)

add_executable(kepler-borrowed-latch-options-tests borrowed_latch_options_tests.cpp
  ${PROJECT_SOURCE_DIR}/src/python/BorrowedLatchOptions.cpp
  ${PROJECT_SOURCE_DIR}/src/python/PyLatchOptions.cpp)
target_include_directories(kepler-borrowed-latch-options-tests PRIVATE
  ${PROJECT_SOURCE_DIR}/src/python ${PROJECT_SOURCE_DIR}/src/sec)
target_link_libraries(kepler-borrowed-latch-options-tests PRIVATE Python3::Python)
add_test(NAME kepler-formal-borrowed-latch-options-tests COMMAND kepler-borrowed-latch-options-tests)

add_executable(kepler-borrowed-latch-tests borrowed_latch_tests.cpp)
target_link_libraries(kepler-borrowed-latch-tests PRIVATE kepler_borrowed_designs Python3::Python)
najaeda_fixup_consumer(kepler-borrowed-latch-tests)
add_test(NAME kepler-formal-borrowed-latch-tests COMMAND kepler-borrowed-latch-tests)
set_tests_properties(kepler-formal-borrowed-latch-options-tests
  kepler-formal-borrowed-latch-tests PROPERTIES TIMEOUT 120)
