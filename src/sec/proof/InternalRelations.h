// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "kinduction/KInductionProblem.h"
#include "../../config/Config.h"

namespace KEPLER_FORMAL::SEC {

struct InternalRelationOptions {
  bool learnInternalRelations = true;
  bool allowXEqualityInInternalRelations = true;
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
};

// Returns only a jointly inductive subset whose base case holds at raw boot.
// No output property or reset-frontier assumption participates in this proof.
std::vector<std::pair<size_t, size_t>> proveInternalRelations(
    const KInductionProblem& problem,
    const std::vector<InternalRelationCandidate>& candidates,
    const InternalRelationOptions& options,
    Config::SolverType solverType);

}  // namespace KEPLER_FORMAL::SEC
