// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

struct TransitionExprView {
  BoolExpr* expr = nullptr;
  const std::unordered_map<size_t, size_t>* symbolMap = nullptr;
};

class TransitionExprResolver {
 public:
  explicit TransitionExprResolver(const KInductionProblem& problem);

  bool contains(size_t stateSymbol) const;
  BoolExpr* at(size_t stateSymbol) const;
  TransitionExprView expressionView(size_t stateSymbol) const;
  const std::set<size_t>& support(size_t stateSymbol) const;
  const std::vector<BoolExpr*>& encodingPostorder(size_t stateSymbol) const;
  void collectSupportForTargets(
      const std::vector<size_t>& stateSymbols,
      const std::unordered_set<size_t>& knownStateSymbols,
      std::unordered_set<size_t>& stateSupport,
      std::unordered_set<size_t>& allSupport) const;
  size_t nodeCount(size_t stateSymbol) const;
  const std::unordered_set<size_t>& stateSymbols() const;
  const std::unordered_map<size_t, size_t>& primaryByComplement() const;

 private:
  const KInductionProblem& problem_;
  std::unordered_map<size_t, BoolExpr*> eagerByStateSymbol_;
  mutable std::unordered_map<size_t, std::set<size_t>> supportByStateSymbol_;
  mutable std::unordered_map<BoolExpr*, std::vector<BoolExpr*>>
      encodingPostorderByExpr_;
  mutable std::unordered_map<size_t, size_t> nodeCountByStateSymbol_;
  mutable std::unordered_set<size_t> stateSymbols_;
  mutable std::unordered_map<size_t, size_t> primaryByComplement_;
  mutable bool stateSymbolsInitialized_ = false;
  mutable bool primaryByComplementInitialized_ = false;
};

}  // namespace KEPLER_FORMAL::SEC
