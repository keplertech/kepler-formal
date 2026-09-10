# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# Use the SDK belonging to this interpreter's installed provider. A CMake
# override alone must not select a different runtime than Python will import.
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c
    "from najaeda import sdk; print(sdk.get_cmake_dir())"
  RESULT_VARIABLE provider_status
  OUTPUT_VARIABLE provider_cmake_dir
  ERROR_VARIABLE provider_error
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT provider_status EQUAL 0 OR NOT IS_DIRECTORY "${provider_cmake_dir}")
  message(FATAL_ERROR
    "Python builds require the NajaEDA shared-runtime SDK in ${Python3_EXECUTABLE}. "
    "Install the matching NajaEDA development wheel first (see docs/python-api.md).\n"
    "${provider_error}")
endif()
# NO_DEFAULT_PATH does not ignore a preexisting <Package>_DIR cache entry.
# Shadow it locally so both fresh and reused builds select Python's provider.
set(NajaEDA_DIR "${provider_cmake_dir}")
find_package(NajaEDA CONFIG REQUIRED PATHS "${provider_cmake_dir}" NO_DEFAULT_PATH)

# These build-policy targets normally come from the vendored Naja project.
# Do not modify the installed provider's compilation or import another runtime.
foreach(policy_target IN ITEMS coverage_config sanitizers_config)
  if(NOT TARGET "${policy_target}")
    add_library("${policy_target}" INTERFACE)
  endif()
endforeach()
find_path(Boost_INCLUDE_DIR NAMES boost/version.hpp REQUIRED)
set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")
find_package(TBB REQUIRED)
