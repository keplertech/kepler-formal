// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace KEPLER_FORMAL {

class BoolExpr;

struct CnfFormula {
  int numVars = 0;  // DIMACS variable count (1-based IDs)
  int rootLit = 0;  // literal representing the root (signed DIMACS var)
  std::vector<std::vector<int>> clauses;
  std::unordered_map<std::string, int> varNameToDimacs;
};

// Encode a BoolExpr into CNF using Tseitin transformation.
CnfFormula encodeBoolExprToCnf(BoolExpr* root);

// Write CNF in DIMACS format to a stream.
bool writeDimacsCnf(const CnfFormula& cnf, std::ostream& out, bool assertRoot = true);

// Convenience helpers that write to a file.
bool dumpDimacsCnf(const CnfFormula& cnf, const std::string& path, bool assertRoot = true);
bool dumpBoolExprToDimacs(BoolExpr* root, const std::string& path, bool assertRoot = true);

}  // namespace KEPLER_FORMAL
