cmake_minimum_required(VERSION 3.10)
project(formal_strategies_lib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a static library target
add_library(formal_strategies STATIC
    miter/BuildPrimaryOutputClauses.cpp
    miter/MiterStrategy.cpp
)

# Make headers accessible to other targets
target_include_directories(formal_strategies PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/miter ${GLUCOSE_ROOT})
target_link_libraries(formal_strategies PUBLIC kepler_clauses formal_structures naja_dnl glucose naja_snl_visual kepler_formal_utils naja_nl)
