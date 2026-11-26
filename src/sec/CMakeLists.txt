# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

add_library(kepler_sec STATIC
    export/SecBtor2Exporter.cpp
    common/BoolExprUtils.cpp
    common/ProofProblemDebug.cpp
    imc/CraigInterpolatingModelChecker.cpp
    imc/ExactInterpolantSynthesizer.cpp
    imc/IMCEngine.cpp
    kinduction/KInductionEngine.cpp
    kinduction/OutputBatching.cpp
    model/SecNetlistChecks.cpp
    model/SequentialDesignModel.cpp
    pdr/PDREngine.cpp
    proof/DualRailEncoding.cpp
    proof/InternalRelations.cpp
    proof/ProofEngineShared.cpp
    proof/TransitionExprResolver.cpp
    kinduction/BaseCaseSolver.cpp
    kinduction/InductionStepSolver.cpp
    kinduction/SatEncoding.cpp
    strategy/SequentialEquivalenceStrategy.cpp
    property/PropertyCheckStrategy.cpp
    property/PropertyExpression.cpp
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

if(NOT BUILD_KEPLER_PYTHON)
  # Reads the assertions the SystemVerilog frontend captured. It needs the
  # frontend and its parser, which the core proof library must not depend on.
  add_library(kepler_sva STATIC property/SvaPropertyTranslator.cpp)
  target_link_libraries(kepler_sva PUBLIC kepler_sec naja_snl_systemverilog)
endif()
