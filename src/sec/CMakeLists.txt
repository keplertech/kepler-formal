# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_library(kepler_sec STATIC
    common/BoolExprUtils.cpp
    imc/ExactInterpolantSynthesizer.cpp
    imc/IMCEngine.cpp
    kinduction/KInductionEngine.cpp
    model/SequentialDesignModel.cpp
    pdr/PDREngine.cpp
    proof/ProofEngineShared.cpp
    kinduction/BaseCaseSolver.cpp
    kinduction/InductionStepSolver.cpp
    kinduction/SatEncoding.cpp
    strategy/ReachableStateInvariant.cpp
    strategy/StructuralStateInvariant.cpp
    strategy/SequentialEquivalenceStrategy.cpp
)

target_include_directories(kepler_sec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_sec
    PUBLIC
      formal_strategies
      kepler_config
      formal_structures
      naja_dnl
      naja_nl
      kissat
)
