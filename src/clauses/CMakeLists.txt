# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.10)
project(kepler_clauses_lib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a static library target
add_library(kepler_clauses STATIC
    SNLLogicCloud.cpp
    SNLTruthTableTree.cpp
    Tree2BoolExpr.cpp
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    # Keep the ancestry traversal practical in symbolized Debug profiles.
    set_source_files_properties(SNLTruthTableTree.cpp PROPERTIES
        COMPILE_OPTIONS "$<$<CONFIG:Debug>:-O2>"
    )
endif()

# Make headers accessible to other targets
target_include_directories(kepler_clauses PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_clauses PUBLIC formal_structures naja_nl naja_dnl naja_core)
