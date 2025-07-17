cmake_minimum_required(VERSION 3.10)
project(kepler_clauses_lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a static library target
add_library(kepler_clauses STATIC
    SNLTruthTableMerger.cpp
)

# Make headers accessible to other targets
target_include_directories(kepler_clauses PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_clauses formal_structures naja_nl)