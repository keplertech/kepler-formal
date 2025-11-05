cmake_minimum_required(VERSION 3.10)
project(formal_structures_lib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a static library target
add_library(formal_structures STATIC
    BoolExpr.cpp
    BoolExprCache.cpp
)

# Make headers accessible to other targets
target_include_directories(formal_structures PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
