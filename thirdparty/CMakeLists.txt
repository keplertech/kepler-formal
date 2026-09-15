# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: Apache-2.0

add_subdirectory(naja)
add_subdirectory(glucose)
add_subdirectory(yaml-cpp)
#add_subdirectory(kissat)
# cadical is built through the root ExternalProject because upstream does not
# provide a CMake subdirectory target.
#add_subdirectory(cpptrace EXCLUDE_FROM_ALL)
