// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "BoolExpr.h"
#include "SNLTruthTableTree.h"
#include <tbb/concurrent_unordered_map.h>

namespace KEPLER_FORMAL {

/// Convert a truth-table tree directly into a BoolExpr
class Tree2BoolExpr {
 public:
  static BoolExpr* convert(const SNLTruthTableTree& tree,
                                           const std::vector<size_t>& varNames);
  static tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*> iso2boolExpr_;
};

}  // namespace KEPLER_FORMAL
