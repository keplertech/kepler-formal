// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "BoolExprCnfWriter.h"

#include "BoolExpr.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stack>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace KEPLER_FORMAL {

namespace {

struct Frame {
  BoolExpr* expr = nullptr;
  bool visited = false;
  int leftLit = 0;
  int rightLit = 0;
};

struct StableKey {
  uint64_t h0 = 0;
  uint64_t h1 = 0;
  size_t size = 0;
  Op op = Op::NONE;
  size_t varID = 0;
};

uint64_t mix64(uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

bool isCommutative(Op op) {
  return op == Op::AND || op == Op::OR || op == Op::XOR;
}

bool stableLess(const StableKey& a, const StableKey& b) {
  return std::tie(a.h0, a.h1, a.size, a.op, a.varID) <
         std::tie(b.h0, b.h1, b.size, b.op, b.varID);
}

StableKey combineKey(Op op, size_t varID, const StableKey* left, const StableKey* right) {
  StableKey key;
  key.op = op;
  key.varID = varID;
  key.size = 1 + (left ? left->size : 0) + (right ? right->size : 0);
  key.h0 = mix64(static_cast<uint64_t>(op) ^ (mix64(varID) << 1));
  key.h1 = mix64((static_cast<uint64_t>(op) << 32) ^ mix64(varID + 0x517cc1b727220a95ULL));
  if (left) {
    key.h0 = mix64(key.h0 ^ left->h0 ^ (left->h1 << 1));
    key.h1 = mix64(key.h1 ^ left->h1 ^ (left->h0 >> 1));
  }
  if (right) {
    key.h0 = mix64(key.h0 ^ right->h0 ^ (right->h1 << 1));
    key.h1 = mix64(key.h1 ^ right->h1 ^ (right->h0 >> 1));
  }
  return key;
}

std::unordered_map<BoolExpr*, StableKey> buildStableKeys(BoolExpr* root) {
  // BoolExpr canonicalization can depend on pointer order; use structural keys
  // here so DIMACS export is stable enough for byte-for-byte regression goldens.
  std::unordered_map<BoolExpr*, StableKey> keys;
  std::stack<Frame> stk;
  stk.push({root, false, 0, 0});

  while (!stk.empty()) {
    Frame fr = stk.top();
    stk.pop();
    BoolExpr* e = fr.expr;
    if (!e || keys.count(e)) {
      continue;
    }
    if (!fr.visited) {
      stk.push({e, true, 0, 0});
      if (e->getRight()) {
        stk.push({e->getRight(), false, 0, 0});
      }
      if (e->getLeft()) {
        stk.push({e->getLeft(), false, 0, 0});
      }
      continue;
    }

    const StableKey* leftKey = nullptr;
    const StableKey* rightKey = nullptr;
    if (e->getLeft()) {
      leftKey = &keys.at(e->getLeft());
    }
    if (e->getRight()) {
      rightKey = &keys.at(e->getRight());
    }
    if (rightKey && isCommutative(e->getOp()) && stableLess(*rightKey, *leftKey)) {
      std::swap(leftKey, rightKey);
    }
    keys.emplace(e, combineKey(e->getOp(), e->getId(), leftKey, rightKey));
  }

  return keys;
}

}  // namespace

CnfFormula encodeBoolExprToCnf(BoolExpr* root) {
  if (!root) {
    throw std::invalid_argument("encodeBoolExprToCnf: root is null");
  }

  CnfFormula cnf;
  std::unordered_map<BoolExpr*, int> node2lit;
  std::unordered_set<std::string> forcedConstants;
  const auto stableKeys = buildStableKeys(root);

  int nextVar = 0; // DIMACS variables are 1-based.

  auto allocVar = [&]() -> int {
    return ++nextVar;
  };

  auto getOrCreateVar = [&](const std::string& key) -> int {
    auto it = cnf.varNameToDimacs.find(key);
    if (it != cnf.varNameToDimacs.end()) {
      return it->second;
    }
    int v = allocVar();
    cnf.varNameToDimacs.emplace(key, v);
    return v;
  };

  auto constVar = [&](bool value) -> int {
    const std::string key = value ? "$__CONST_TRUE__" : "$__CONST_FALSE__";
    int v = getOrCreateVar(key);
    if (forcedConstants.insert(key).second) {
      cnf.clauses.push_back({value ? v : -v});
    }
    return v;
  };

  std::stack<Frame> stk;
  stk.push({root, false, 0, 0});

  while (!stk.empty()) {
    Frame& fr = stk.top();
    BoolExpr* e = fr.expr;

    if (node2lit.count(e)) {
      stk.pop();
      continue;
    }

    if (!fr.visited && e->getOp() == Op::VAR) {
      int lit = 0;
      if (e->getId() == 0) {
        lit = constVar(false);
      } else if (e->getId() == 1) {
        lit = constVar(true);
      } else {
        lit = getOrCreateVar(e->getName());
      }
      node2lit[e] = lit;
      stk.pop();
      continue;
    }

    if (!fr.visited) {
      fr.visited = true;
      BoolExpr* left = e->getLeft();
      BoolExpr* right = e->getRight();
      if (right && isCommutative(e->getOp()) &&
          stableLess(stableKeys.at(right), stableKeys.at(left))) {
        std::swap(left, right);
      }
      if (right) {
        stk.push({right, false, 0, 0});
      }
      if (left) {
        stk.push({left, false, 0, 0});
      }
      continue;
    }

    BoolExpr* left = e->getLeft();
    BoolExpr* right = e->getRight();
    if (right && isCommutative(e->getOp()) &&
        stableLess(stableKeys.at(right), stableKeys.at(left))) {
      std::swap(left, right);
    }
    if (left) {
      fr.leftLit = node2lit[left];
    }
    if (right) {
      fr.rightLit = node2lit[right];
    }

    int v = allocVar();
    int lit = v;
    node2lit[e] = lit;

    switch (e->getOp()) {
      case Op::NOT:
        cnf.clauses.push_back({-lit, -fr.leftLit});
        cnf.clauses.push_back({ lit,  fr.leftLit});
        break;

      case Op::AND:
        cnf.clauses.push_back({-lit, fr.leftLit});
        cnf.clauses.push_back({-lit, fr.rightLit});
        cnf.clauses.push_back({ lit, -fr.leftLit, -fr.rightLit});
        break;

      case Op::OR:
        cnf.clauses.push_back({-fr.leftLit,  lit});
        cnf.clauses.push_back({-fr.rightLit, lit});
        cnf.clauses.push_back({-lit,         fr.leftLit, fr.rightLit});
        break;

      case Op::XOR:
        cnf.clauses.push_back({-lit, -fr.leftLit, -fr.rightLit});
        cnf.clauses.push_back({-lit,  fr.leftLit,  fr.rightLit});
        cnf.clauses.push_back({ lit, -fr.leftLit,  fr.rightLit});
        cnf.clauses.push_back({ lit,  fr.leftLit, -fr.rightLit});
        break;

      case Op::VAR:
      case Op::NONE:
      default:
        throw std::runtime_error("encodeBoolExprToCnf: unsupported op");
    }

    stk.pop();
  }

  cnf.numVars = nextVar;
  cnf.rootLit = node2lit[root];
  return cnf;
}

bool writeDimacsCnf(const CnfFormula& cnf, std::ostream& out, bool assertRoot) {
  if (!out.good()) {
    return false;
  }
  int extraClauses = assertRoot ? 1 : 0;
  out << "p cnf " << cnf.numVars << " " << (cnf.clauses.size() + extraClauses) << "\n";
  for (const auto& clause : cnf.clauses) {
    for (int lit : clause) {
      out << lit << " ";
    }
    out << "0\n";
  }
  if (assertRoot) {
    out << cnf.rootLit << " 0\n";
  }
  return out.good();
}

bool dumpDimacsCnf(const CnfFormula& cnf, const std::string& path, bool assertRoot) {
  std::ofstream out(path);
  if (!out.is_open()) {
    return false;
  }
  return writeDimacsCnf(cnf, out, assertRoot);
}

bool dumpBoolExprToDimacs(BoolExpr* root, const std::string& path, bool assertRoot) {
  CnfFormula cnf = encodeBoolExprToCnf(root);
  return dumpDimacsCnf(cnf, path, assertRoot);
}

}  // namespace KEPLER_FORMAL
