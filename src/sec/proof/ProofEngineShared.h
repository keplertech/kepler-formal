// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <unordered_map>
#include <vector>

#include "../../config/Config.h"
#include "BoolExpr.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

// Shared helpers used by the explicit SEC proof engines. They keep the
// low-level SAT/invariant checks in one place so K_INDUCTION, IMC, and PDR all
// reason about the same extracted transition system.

using FormulaSupportCache = std::unordered_map<BoolExpr*, std::vector<size_t>>;

BoolExpr* buildProofInitFormula(const KInductionProblem& problem);

size_t nextFreshProofSymbol(const KInductionProblem& problem);

std::unordered_map<size_t, size_t> allocateFreshProofSymbols(
    const std::vector<size_t>& originalSymbols,
    size_t& nextSymbol);

BoolExpr* buildOneStepTransitionFormula(
    const KInductionProblem& problem,
    const std::unordered_map<size_t, size_t>& nextStateSymbols,
    BoolExpr* certifiedStateConstraints = nullptr);

BoolExpr* buildCurrentStateLegalityFormula(const KInductionProblem& problem);

BoolExpr* remapProofFormula(
    BoolExpr* formula,
    const std::unordered_map<size_t, size_t>& symbolMap);

bool isProofFormulaSatisfiable(
    BoolExpr* formula,
    KEPLER_FORMAL::Config::SolverType solverType);

bool initialFrontierImplies(
    BoolExpr* initFormula,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType);

BoolExpr* selectValidatedStrengtheningInvariant(
    const KInductionProblem& problem,
    BoolExpr* initFormula,
    KEPLER_FORMAL::Config::SolverType solverType);

bool invariantExcludesBadStates(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType,
    BoolExpr* certifiedStateConstraints = nullptr);

// Optional constraints must already be proved invariants of the original
// transition system. Omission preserves the unconstrained proof query.
bool isInductiveInvariant(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType,
    BoolExpr* certifiedStateConstraints = nullptr);

bool isInductiveInvariant(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType,
    FormulaSupportCache& supportCache,
    BoolExpr* certifiedStateConstraints = nullptr);

}  // namespace KEPLER_FORMAL::SEC
