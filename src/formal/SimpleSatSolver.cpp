#include "SimpleSatSolver.h"

using namespace KEPLER_FORMAL;

bool SimplSatSolver::solve(
    std::unordered_map<std::string, bool>& outAssignment) {
  extractVarNames();
  uint64_t total = 1ULL << numVars_;  // 2^n assignments
  std::unordered_map<std::string, bool> env;
  for (uint64_t mask = 0; mask < total; ++mask) {
    // build one assignment from bits of mask
    env.clear();
    for (size_t i = 0; i < numVars_; ++i) {
      env[varNames_[i]] = ((mask >> i) & 1) != 0;
    }
    // check all clauses under this assignment
    bool allTrue = true;
    for (auto& clause : clauses_) {
      if (!clause.evaluate(env)) {
        allTrue = false;
        break;
      }
    }
    if (allTrue) {
      outAssignment = env;
      return true;
    }
  }
  return false;  // no assignment satisfied all clauses
}

// Recursively collect variable names from a BoolExpr tree
void SimplSatSolver::collectVars(const BoolExpr& e,
                                 std::unordered_set<std::string>& s) {
  switch (e.getOp()) {  // requires solver to be in same namespace/file
    case Op::VAR:
      s.insert(e.getName());
      break;
    case Op::NOT:
      collectVars(*e.getLeft(), s);
      break;
    case Op::AND:
    case Op::OR:
    case Op::XOR:
      collectVars(*e.getLeft(), s);
      collectVars(*e.getRight(), s);
      break;
  }
}

// Extract and sort the unique var names from all clauses
void SimplSatSolver::extractVarNames() {
  std::unordered_set<std::string> vs;
  for (auto& c : clauses_) {
    collectVars(c, vs);
  }
  // Insert the set content into the varNames_ vector
  varNames_.clear();
  varNames_.insert(varNames_.end(), vs.begin(), vs.end());
  std::sort(varNames_.begin(), varNames_.end());
  numVars_ = varNames_.size();
}