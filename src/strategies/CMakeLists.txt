# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_library(formal_strategies STATIC
    miter/BuildPrimaryOutputClauses.cpp
    miter/MiterStrategy.cpp
    ${KISSAT_ROOT}/src
)

# Make headers accessible to other targets
target_include_directories(formal_strategies PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/miter ${GLUCOSE_ROOT} ${KISSAT_ROOT}/src )
target_link_libraries(formal_strategies PUBLIC kepler_config kepler_clauses formal_structures naja_dnl glucose naja_snl_visual kepler_formal_utils naja_nl kissat)
