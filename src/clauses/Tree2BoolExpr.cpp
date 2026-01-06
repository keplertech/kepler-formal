// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "Tree2BoolExpr.h"
#include "BoolExpr.h"
#include "DNL.h"
#include "SNLTruthTable.h"
#include "SNLTruthTableTree.h"
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb_allocator.h>
#include <bitset>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

// #define DEBUG_CHECKS
// #define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace naja::NL;
using namespace KEPLER_FORMAL;

typedef std::pair<std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>, size_t> TermsPair;
tbb::enumerable_thread_specific<TermsPair> termsETS;
tbb::concurrent_vector<TermsPair*> termsETSvector;

void initTermsETS() {
  if (termsETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
    for (size_t i = termsETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
      termsETSvector.push_back(nullptr);
    }
  }
  if (termsETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
    termsETSvector[tbb::this_task_arena::current_thread_index()] = &termsETS.local();
  }
}

TermsPair& getTErmsETS() {
  initTermsETS();
  size_t idx = tbb::this_task_arena::current_thread_index();
  if (idx >= termsETSvector.size() || termsETSvector[idx] == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error("getTErmsETS: not initialized for this thread");
    // LCOV_EXCL_STOP
  }
  return *termsETSvector[idx];
}

size_t sizeOfTermsETS() {
  return getTErmsETS().second;
}

void clearTermsETS() {
  // If vector reach size larger than 1024, clear it to save memory
  auto& termsLocal = getTErmsETS();
  if (termsLocal.first.size() > 1024) {
    termsLocal.first.clear();
  }
  termsLocal.second = 0;
}

void pushBackTermsETS(const std::shared_ptr<BoolExpr>& term) {
  auto& termsLocal = getTErmsETS();
  auto& vec = termsLocal.first;
  auto& sz = termsLocal.second;
  if (vec.size() > sz) {
    vec[sz] = term;
    sz++;
    return;
  }
  vec.push_back(term);
  sz++;
}

void reserveTermsETS(size_t n) {
  auto& termsLocal = getTErmsETS();
  if (termsLocal.first.size() >= n) return;
  termsLocal.first.reserve(n);
}

bool emptyTermsETS() {
  return getTErmsETS().second == 0;
}

// same for std::vector<bool, tbb::tbb_allocator<bool>>
typedef std::pair<std::vector<bool, tbb::tbb_allocator<bool>>, size_t> RelevantPair;
tbb::enumerable_thread_specific<RelevantPair> relevantETS;
tbb::concurrent_vector<RelevantPair*> relevantETSvector;

void initRelevantETS() {
  if (relevantETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
    for (size_t i = relevantETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
      relevantETSvector.push_back(nullptr);
    }
  }
  if (relevantETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
    relevantETSvector[tbb::this_task_arena::current_thread_index()] = &relevantETS.local();
  }
}

RelevantPair& getRelevantETS() {
  initRelevantETS();
  size_t idx = tbb::this_task_arena::current_thread_index();
  if (idx >= relevantETSvector.size() || relevantETSvector[idx] == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error("getRelevantETS: not initialized for this thread");
    // LCOV_EXCL_STOP
  }
  return *relevantETSvector[idx];
}

size_t sizeOfRelevantETS() {
  return getRelevantETS().second;
}

void clearRelevantETS() {
  // If vector reach size larger than 1024, clear it to save memory
  auto& relevantLocal = getRelevantETS();
  if (relevantLocal.first.size() > 1024) {
    relevantLocal.first.clear();
  }
  relevantLocal.second = 0;
}

// void pushBackRelevantETS(bool b) {
//   auto& relevantLocal = getRelevantETS();
//   auto& vec = relevantLocal.first;
//   auto& sz = relevantLocal.second;
//   if (vec.size() > sz) {
//     vec[sz] = b;
//     sz++;
//     return;
//   }
//   vec.push_back(b);
//   sz++;
// }

void setRelevantETS(size_t i, bool b) {
  auto& relevantLocal = getRelevantETS();
  if (i >= relevantLocal.second) {
    // LCOV_EXCL_START
    assert(false && "setRelevantETS: index out of range");
    // LCOV_EXCL_STOP
  }
  relevantLocal.first[i] = b;
}

bool getRelevantETS(size_t i) {
  auto& relevantLocal = getRelevantETS();
  if (i >= relevantLocal.second) {
    // LCOV_EXCL_START
    throw std::out_of_range("getRelevantETS: index out of range");
    // LCOV_EXCL_STOP
  }
  return relevantLocal.first[i];
}

void reserveRelevantETSwithFalse(size_t n) {
  auto& relevantLocal = getRelevantETS();
  auto& vec = relevantLocal.first;
  auto& sz = relevantLocal.second;
  if (vec.size() >= n) {
    vec.assign(n, false);
    sz = n;
    return;
  }
  size_t oldSize = vec.size();
  vec.resize(n, false);
  vec.assign(n, false);
  sz = n;
}

// do same for std::vector<std::shared_ptr<BoolExpr>,
// tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> memo;
typedef std::pair<std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>, size_t> MemoPair;
tbb::enumerable_thread_specific<MemoPair> memoETS;
tbb::concurrent_vector<MemoPair*> memoETSvector;

void initMemoETS() {
  if (memoETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
    for (size_t i = memoETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
      memoETSvector.push_back(nullptr);
    }
  }
  if (memoETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
    memoETSvector[tbb::this_task_arena::current_thread_index()] = &memoETS.local();
  }
}

MemoPair& getMemoETS() {
  initMemoETS();
  size_t idx = tbb::this_task_arena::current_thread_index();
  if (idx >= memoETSvector.size() || memoETSvector[idx] == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error("getMemoETS: not initialized for this thread");
    // LCOV_EXCL_STOP
  }
  return *memoETSvector[idx];
}

size_t sizeOfMemoETS() {
  return getMemoETS().second;
}

void clearMemoETS() {
  // If vector reach size larger than 1024, clear it to save memory
  auto& memoLocal = getMemoETS();
  if (memoLocal.first.size() > 1024) {
    memoLocal.first.clear();
  }
  memoLocal.second = 0;
}

// void pushBackMemoETS(const std::shared_ptr<BoolExpr>& expr) {
//   auto& memoLocal = getMemoETS();
//   auto& vec = memoLocal.first;
//   auto& sz = memoLocal.second;
//   if (vec.size() > sz) {
//     vec[sz] = expr;
//     sz++;
//     return;
//   }
//   vec.push_back(expr);
//   sz++;
// }

void reserveMemoETS(size_t n) {
  auto& memoLocal = getMemoETS();
  auto& vec = memoLocal.first;
  auto& sz = memoLocal.second;
  if (vec.size() >= n) {
    sz = n;
    vec.assign(n, nullptr);
    return;
  }
  vec.resize(n);
  sz = n;
  vec.assign(n, nullptr);
}

void setMemoETS(size_t i, const std::shared_ptr<BoolExpr>& expr) {
  auto& memoLocal = getMemoETS();
  if (i >= memoLocal.second) {
    assert(false && "setMemoETS: index out of range");
  }
  memoLocal.first[i] = expr;
}

const std::shared_ptr<BoolExpr>& getMemoETS(size_t i) {
  auto& memoLocal = getMemoETS();
  if (i >= memoLocal.second) {
    assert(false && "getMemoETS: index out of range");
  }
  return memoLocal.first[i];
}

// same for std::vector<std::shared_ptr<BoolExpr>>,
// tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> childF;
typedef std::pair<std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>, size_t> ChildFETSPair;
tbb::enumerable_thread_specific<ChildFETSPair> childFETS;
tbb::concurrent_vector<ChildFETSPair*> childFETSvector;

void initChildFETS() {
  if (childFETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
    for (size_t i = childFETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
      childFETSvector.push_back(nullptr);
    }
  }
  if (childFETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
    childFETSvector[tbb::this_task_arena::current_thread_index()] = &childFETS.local();
  }
}

ChildFETSPair& getChildFETS() {
  initChildFETS();
  size_t idx = tbb::this_task_arena::current_thread_index();
  if (idx >= childFETSvector.size() || childFETSvector[idx] == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error("getChildFETS: not initialized for this thread");
    // LCOV_EXCL_STOP
  }
  return *childFETSvector[idx];
}

size_t sizeOfChildFETS() {
  return getChildFETS().second;
}

void clearChildFETS() {
  // If vector reach size larger than 1024, clear it to save memory
  auto& childLocal = getChildFETS();
  if (childLocal.first.size() > 1024) {
    childLocal.first.clear();
  }
  childLocal.second = 0;
}

// void pushBackChildFETS(const std::shared_ptr<BoolExpr>& expr) {
//   auto& childLocal = getChildFETS();
//   auto& vec = childLocal.first;
//   auto& sz = childLocal.second;
//   if (vec.size() > sz) {
//     vec[sz] = expr;
//     sz++;
//     return;
//   }
//   vec.push_back(expr);
//   sz++;
// }

void reserveChildFETS(size_t n) {
  auto& childLocal = getChildFETS();
  auto& vec = childLocal.first;
  auto& sz = childLocal.second;
  if (vec.size() >= n) {
    sz = n;
    vec.assign(n, nullptr);
    return;
  }
  vec.resize(n);
  sz = n;
  vec.assign(n, nullptr);
}

const std::shared_ptr<BoolExpr>& getChildFETS(size_t i) {
  auto& childLocal = getChildFETS();
  if (i >= childLocal.second) {
    assert(false && "getChildFETS: index out of range");
  }
  return childLocal.first[i];
}

void setChildFETS(size_t i, const std::shared_ptr<BoolExpr>& expr) {
  auto& childLocal = getChildFETS();
  if (i >= childLocal.second) {
    assert(false && "setChildFETS: index out of range");
  }
  childLocal.first[i] = expr;
}

// size_t toSizeT(const std::string& s) {
//   if (s.empty()) {
//     assert(false && "toSizeT: empty string");
//   }
//   size_t result = 0;
//   const size_t max = std::numeric_limits<size_t>::max();
//   for (unsigned char uc : s) {
//     if (!std::isdigit(uc)) {
//       throw std::invalid_argument("toSizeT: invalid character '" + std::string(1, static_cast<char>(uc)) + "' in input");
//     }
//     size_t digit = static_cast<size_t>(uc - '0');
//     // Check for overflow: result * 10 + digit > max
//     if (result > (max - digit) / 10) {
//       throw std::out_of_range("toSizeT: value out of range for size_t");
//     }
//     result = result * 10 + digit;
//   }
//   return result;
// }

// Fold a list of literals into a single AND
// static std::shared_ptr<BoolExpr> mkAnd(
//   const std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>& lits) {
//   if (lits.empty()) return BoolExpr::createTrue();
//   auto cur = lits[0];
//   for (size_t i = 1; i < lits.size(); ++i) cur = BoolExpr::And(cur, lits[i]);
//   return cur;
// }

// // Fold a list of terms into a single OR
// static std::shared_ptr<BoolExpr> mkOr(
//   const std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>& terms) {
//   if (terms.empty()) return BoolExpr::createFalse();
//   auto cur = terms[0];
//   for (size_t i = 1; i < terms.size(); ++i) cur = BoolExpr::Or(cur, terms[i]);
//   return cur;
// }

std::shared_ptr<BoolExpr> Tree2BoolExpr::convert(
  const SNLTruthTableTree& tree, const std::vector<size_t>& varNames) {

  initChildFETS();
  initMemoETS();
  initRelevantETS();
  initTermsETS();

  const auto root = tree.getRoot();
  if (!root) return nullptr;

  // 1) find maxID
  // size_t maxID = 0;
  // {
  // using NodePtr = const SNLTruthTableTree::Node*;
  // std::vector<NodePtr, tbb::tbb_allocator<NodePtr>> dfs;
  // dfs.reserve(128);
  // dfs.push_back(root.get());
  // while (!dfs.empty()) {
  // NodePtr n = dfs.back();
  // dfs.pop_back();
  // maxID = std::max(maxID, (size_t) n->nodeID);
  // if (n->type == SNLTruthTableTree::Node::Type::Table || // n->type == SNLTruthTableTree::Node::Type::P)
  // {
  // for (const auto& c : n->childrenIds) {
  // assert(n->tree->nodeFromId(c).get() != nullptr);
  // dfs.push_back(n->tree->nodeFromId(c).get());
  // }
  // }
  // }
  // }

  size_t maxID = tree.getMaxID();

  // 2) memo table
  clearMemoETS();
  reserveMemoETS(maxID + 1);

  // 3) post-order build
  using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
  std::vector<Frame, tbb::tbb_allocator<Frame>> stack;
  stack.reserve(maxID + 1);
  stack.emplace_back(root.get(), false);

  while (!stack.empty()) {
    Frame f = stack.back();
    stack.pop_back();
    const SNLTruthTableTree::Node* node = f.first;
    bool visited = f.second;
    size_t id = node->nodeID;

    if (!visited) {
      if (getMemoETS(id) != nullptr) continue;
      if (node->type == SNLTruthTableTree::Node::Type::Table || node->type == SNLTruthTableTree::Node::Type::P) {
        stack.emplace_back(node, true);
        for (const auto& c : node->childrenIds) stack.emplace_back(node->tree->nodeFromId(c).get(), false);
      } else {
        assert(node->type == SNLTruthTableTree::Node::Type::Input);
        if (node->parentIds.size() > 1) {
          #ifdef DEBUG_PRINTS
          for (const auto& pid : node->parentIds) {
            DEBUG_LOG("%s\n", naja::DNL::get()->getDNLTerminalFromID(tree.nodeFromId(pid)->data.termid)
                     .getSnlBitTerm()->getString().c_str());
            DEBUG_LOG("of model %s\n", naja::DNL::get()->getDNLTerminalFromID(tree.nodeFromId(pid)->data.termid)
                   .getDNLInstance().getSNLModel()->getString().c_str());
          }
          #endif
        }
        if (node->parentIds.empty()) { 
          // LCOV_EXCL_START
          throw std::runtime_error("Input node has no parent"); 
          // LCOV_EXCL_STOP
        }
        assert(node->parentIds.size() == 1);
        auto parent = node->tree->nodeFromId(node->parentIds[0]);
        assert(parent && parent->type == SNLTruthTableTree::Node::Type::P);
        if (parent->data.termid >= varNames.size()) {
          DEBUG_LOG("varNames size: %zu, parent data.termid: %zu\n", varNames.size(), (size_t)parent->data.termid);
          assert(parent->data.termid < varNames.size());
        }
        if (varNames[parent->data.termid] == (size_t)-1) {
          // LCOV_EXCL_START
          throw std::runtime_error("Input variable index is SIZE_MAX");
          // LCOV_EXCL_STOP
        }
        if (varNames[parent->data.termid] == 0) {
           setMemoETS(id, BoolExpr::createFalse());
        } else if (varNames[parent->data.termid] == 1) {
           setMemoETS(id, BoolExpr::createTrue());
        } else {
          setMemoETS(id, BoolExpr::Var(varNames[parent->data.termid]));
        }
      }
    } else {
      // post-visit for Table / P
      const SNLTruthTable& tbl = node->getTruthTable();
      uint32_t k = tbl.size();
      uint64_t rows = uint64_t{1} << k;

      if (tbl.all0()) {
        setMemoETS(id, BoolExpr::createFalse());
      } else if (tbl.all1()) {
        setMemoETS(id, BoolExpr::createTrue());
      } else {
        // gather children
        clearChildFETS();
        reserveChildFETS(k);
        for (uint32_t i = 0; i < k; ++i) {
          size_t cid = node->tree->nodeFromId(node->childrenIds[i])->nodeID;
          setChildFETS(i, getMemoETS(cid));
        }

        // find which inputs actually matter
        clearRelevantETS();
        reserveRelevantETSwithFalse(k);
        for (uint32_t j = 0; j < k; ++j) {
          for (uint64_t m = 0; m < rows; ++m) {
            bool b0 = tbl.bits().bit(m);
            bool b1 = tbl.bits().bit(m ^ (uint64_t{1} << j));
            if (b0 != b1) { setRelevantETS(j, true); break; }
          }
        }

        // collect the indices of relevant vars
        std::vector<uint32_t, tbb::tbb_allocator<uint32_t>> relIdx;
        for (uint32_t j = 0; j < k; ++j) { if (getRelevantETS(j)) relIdx.push_back(j); }

        // if nothing matters, fall back to constant-false
        if (relIdx.empty()) {
          setMemoETS(id, BoolExpr::createFalse());
        } else {
          // build the DNF terms
          clearTermsETS();
          reserveTermsETS(static_cast<size_t>(rows));
          for (uint64_t m = 0; m < rows; ++m) {
            if (!tbl.bits().bit(m)) continue;
            std::shared_ptr<BoolExpr> term = nullptr;
            bool firstLit = true;
            std::shared_ptr<BoolExpr> lit = nullptr;
            for (uint32_t j : relIdx) {
              bool bit1 = ((m >> j) & 1) != 0;
              lit = bit1 ? getChildFETS(j) : BoolExpr::Not(getChildFETS(j));
              if (firstLit) { term = lit; firstLit = false; }
              else { assert(term != nullptr); assert(lit != nullptr); term = BoolExpr::And(term, lit); }
            }
            // only push if we actually got a literal
            if (term) { pushBackTermsETS(std::move(term)); }
          }

          // guard against an empty terms list
          if (emptyTermsETS()) { setMemoETS(id, BoolExpr::createFalse()); }
          else {
            // fold into OR
            std::shared_ptr<BoolExpr> expr = getTErmsETS().first[0];
            for (size_t t = 1; t < sizeOfTermsETS(); ++t) {
              expr = BoolExpr::Or(expr, getTErmsETS().first[t]);
            }
            setMemoETS(id, expr);
          }
        }
      }
    }
  }

  // 4) return root
  return getMemoETS(root->nodeID);
}
