# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-2.0

cmake_minimum_required(VERSION 3.10)
project(kepler_formal_utils_lib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a static library target
add_library(kepler_formal_utils STATIC
    SNLLogicCone.cpp
)

# Make headers accessible to other targets
target_include_directories(kepler_formal_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_formal_utils PUBLIC formal_structures naja_nl naja_dnl)