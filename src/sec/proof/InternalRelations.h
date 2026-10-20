// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "kinduction/KInductionProblem.h"
#include "../../config/Config.h"

namespace KEPLER_FORMAL::SEC {

struct InternalRelationOptions {
  bool learnInternalRelations = true;
  bool allowXEqualityInInternalRelations = true;
  // Dual-rail candidates the rail check refutes are retried on the binary
  // next-state functions: a pair equal in the common initial state (the same
  // value, or unknown on both sides) and inductive there is equal from every
  // common start (van Eijk, TCAD 2000; Savoj et al., FMCAD 2010). Design 2's
  // register then takes design 1's transition in the proof problem.
  bool binaryRelationRound = true;
};

struct SequentialDesignModel;
void learnInternalStateRelations(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    KInductionProblem& problem,
    const InternalRelationOptions& options,
    Config::SolverType solverType,
    bool diagnostics);

// Van Eijk's merge step (Sequential Equivalence Checking Based on Structural
// Similarities, IEEE TCAD 19(7), 2000): read every register through the
// representative of its proved equality class, hash both designs' observed
// output cones under that merge, and discharge each output whose two cones
// become one formula.  Discharged outputs get a trivially true obligation and
// are flagged in problem.structurallyProvedOutputs; the property, bad and
// induction formulas are rebuilt.  Returns the number of outputs discharged.
size_t dischargeStructurallyEquivalentOutputs(KInductionProblem& problem,
                                              bool diagnostics);

// A candidate relates complete values, never independently matched X rails.
struct InternalRelationCandidate {
  std::vector<std::pair<size_t, size_t>> equalities;
  std::vector<DualRailSymbolPair> values;
  // Display names of the two related registers, for diagnostics only.
  std::string label;
};

// Merges every proved same-frame state equality into one representative symbol
// and rewrites the eager transitions and output formulas under that merge, so
// identical logic over related registers becomes one hash-consed node (van
// Eijk, TCAD 2000; Kuehlmann et al., TCAD 2002). Lazily remapped problems are
// left unchanged.
void mergeProvedStateRelations(KInductionProblem& problem, bool diagnostics);

// Returns only a jointly inductive subset whose base case holds at raw boot.
// No output property or reset-frontier assumption participates in this proof.
std::vector<std::pair<size_t, size_t>> proveInternalRelations(
    const KInductionProblem& problem,
    const std::vector<InternalRelationCandidate>& candidates,
    const InternalRelationOptions& options,
    Config::SolverType solverType,
    bool commonStartBinary = false);

}  // namespace KEPLER_FORMAL::SEC
