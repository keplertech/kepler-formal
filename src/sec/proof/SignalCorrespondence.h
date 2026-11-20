// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "../../config/Config.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

// The set F of van Eijk (TCAD 2000, section III): the signals of both designs
// among which correspondences are looked for. Any set is sound; a larger one
// proves more and costs more.
enum class SignalCorrespondenceScope {
  // Every register, every compared output and every register's next-state
  // function.
  RegistersAndOutputs,
  // Every node of the next-state and output logic.
  AllSignals,
};

struct SignalCorrespondenceResult {
  // False when the pass gave up (size or solver budget, or an operator it does
  // not read); nothing is reported then.
  bool completed = false;
  // State symbol pairs that are equal in every frame, the boot state and the
  // reset window included, so they can be handed to the engines like the
  // pairs of the register learner.
  std::vector<std::pair<size_t, size_t>> stateEqualities;
  // One flag per observed output: its two cones have the same value in every
  // frame from the first compared one on.
  std::vector<bool> provedOutputs;
  size_t signals = 0;
  size_t classes = 0;
  size_t rounds = 0;
  size_t solverCalls = 0;
};

// Signal correspondence (C. A. J. van Eijk, "Sequential equivalence checking
// based on structural similarities", IEEE TCAD 19(7), 2000), decided with a
// SAT solver instead of BDDs (Bjesse and Claessen, "SAT-based verification
// without state space traversal", FMCAD 2000).
//
// Signals are grouped by random sequential simulation from the start state,
// each normalised in polarity so that signals holding opposite values share a
// class. A class is kept when its members are equal in the start state for
// every input (eq. 4) and when, with every class assumed in one time frame,
// its members are equal in the next (eq. 5); a counterexample splits the
// classes and the step is repeated until nothing splits. What is left is a
// signal correspondence relation, so its members are equal in every reachable
// state.
//
// The start state is the boot state, or with a reset bootstrap the state the
// reset window leaves, computed by propagating constants through the window
// with every other input unknown. After the window the reset ports are held
// at their inactive value, as in the engines.
//
// Assumed equalities are plain constraints of the current frame. Nothing is
// merged before it is proved.
//
// In the dual-rail encoding a class may relate signals that are unknown, X on
// both sides. A caller that disallows X equality must not run the pass.
SignalCorrespondenceResult proveSignalCorrespondence(
    const KInductionProblem& problem,
    SignalCorrespondenceScope scope,
    Config::SolverType solverType);

struct SignalCorrespondenceLearning {
  // One flag per observed output.
  std::vector<bool> provedOutputs;
  std::vector<std::pair<size_t, size_t>> stateEqualities;
};

// Runs the pass for the outputs no earlier step proved: first on registers
// and outputs, then on all signals if outputs remain. The problem is not
// changed; the caller discharges the proved outputs and records the state
// equalities.
SignalCorrespondenceLearning learnSignalCorrespondence(
    const KInductionProblem& problem,
    Config::SolverType solverType,
    bool diagnostics);

}  // namespace KEPLER_FORMAL::SEC
