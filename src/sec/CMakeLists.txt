# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_library(kepler_sec STATIC
    common/BoolExprUtils.cpp
    common/ProofProblemDebug.cpp
    imc/CraigInterpolatingModelChecker.cpp
    imc/ExactInterpolantSynthesizer.cpp
    imc/IMCEngine.cpp
    kinduction/KInductionEngine.cpp
    kinduction/OutputBatching.cpp
    model/SequentialDesignModel.cpp
    pdr/PDREngine.cpp
    proof/DualRailEncoding.cpp
    proof/ProofEngineShared.cpp
    proof/TransitionExprResolver.cpp
    kinduction/BaseCaseSolver.cpp
    kinduction/InductionStepSolver.cpp
    kinduction/SatEncoding.cpp
    strategy/SequentialEquivalenceStrategy.cpp
)

target_include_directories(kepler_sec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_sec
    PUBLIC
      formal_strategies
      kepler_clocks
      kepler_config
      formal_structures
      naja_dnl
      naja_nl
      kissat
      cadical
)
