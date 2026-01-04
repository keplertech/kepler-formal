# Copyright 2024-2025 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.10)
project(scope_extraction_lib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a static library target
add_library(scope_extraction STATIC
    ScopeExtraction.cpp
)

# Make headers accessible to other targets
target_include_directories(scope_extraction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(scope_extraction PUBLIC naja_nl naja_dnl naja_opt kepler_formal_utils)
