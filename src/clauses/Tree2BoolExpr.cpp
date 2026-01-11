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

tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*> Tree2BoolExpr::iso2boolExpr_ =
    tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*>();

typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> TermsPair;

thread_local TermsPair termsETS;

TermsPair& getTErmsETS() {
  return termsETS;
}

size_t sizeOfTermsETS() {
  return getTErmsETS().second;
}

void clearTermsETS() {
  auto& termsLocal = getTErmsETS();
  termsLocal.first.clear();
}

void pushBackTermsETS(BoolExpr* term) {
  getTErmsETS().first.emplace_back(term);
}

bool emptyTermsETS() {
  return getTErmsETS().first.empty();
}

typedef std::pair<std::vector<uint8_t, tbb::tbb_allocator<uint8_t>>, size_t> RelevantPair;

thread_local RelevantPair relevantETS;

RelevantPair& getRelevantETS() {
  return relevantETS;
}

size_t sizeOfRelevantETS() {
  return getRelevantETS().second;
}

void clearRelevantETS() {
  auto& relevantLocal = getRelevantETS();
  relevantLocal.first.clear();
}

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

typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> MemoPair;

thread_local MemoPair memoETS;

MemoPair& getMemoETS() {
  return memoETS;
}

void clearMemoETS() {
  auto& memoLocal = getMemoETS();
  memoLocal.first.clear();
}

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

void setMemoETS(size_t i, BoolExpr* expr) {
  auto& memoLocal = getMemoETS();
  assert(i < memoLocal.second && "setMemoETS: index out of range");
  memoLocal.first[i] = expr;
}

BoolExpr* getMemoETS(size_t i) {
  auto& memoLocal = getMemoETS();
  assert(i < memoLocal.second && "getMemoETS: index out of range");
  return memoLocal.first[i];
}

typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> ChildFETSPair;

thread_local ChildFETSPair childFETS;

ChildFETSPair& getChildFETS() {
  return childFETS;
}

void clearChildFETS() {
  auto& childLocal = getChildFETS();
  childLocal.second = 0;
}

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

BoolExpr* getChildFETS(size_t i) {
  auto& childLocal = getChildFETS();
  assert(i < childLocal.second && "getChildFETS: index out of range");
  return childLocal.first[i];
}

void setChildFETS(size_t i, BoolExpr* expr) {
  auto& childLocal = getChildFETS();
  assert(i < childLocal.second && "setChildFETS: index out of range");
  childLocal.first[i] = expr;
}

// same for stuck with frame
using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
thread_local std::vector<Frame, tbb::tbb_allocator<Frame>> stackETS;

std::vector<Frame, tbb::tbb_allocator<Frame>>& getStackETS() {
  return stackETS;
}

BoolExpr* Tree2BoolExpr::convert(
  const SNLTruthTableTree& tree, const std::vector<size_t>& varNames) {

  const auto root = tree.getRoot();
  if (!root) return nullptr;

  size_t maxID = tree.getMaxID();

  // 2) memo table
  clearMemoETS();
  reserveMemoETS(maxID + 1);

  // 3) post-order build
  auto & stack = getStackETS();
  stack.clear();
  //stack.reserve(maxID + 1);
  stack.emplace_back(root.get(), false);

  while (!stack.empty()) {
    Frame f = stack.back();
    stack.pop_back();
    const SNLTruthTableTree::Node* node = f.first;
    naja::DNL::DNLID isoID = naja::DNL::DNLID_MAX;
    if (node->type != SNLTruthTableTree::Node::Type::Input) {
      isoID = naja::DNL::get()->getDNLTerminalFromID(node->data.termid).getIsoID();
    }
    bool visited = f.second;
    size_t id = node->nodeID;

    if (!visited) {
      naja::DNL::DNLID isoID = naja::DNL::DNLID_MAX;
      if (node->type != SNLTruthTableTree::Node::Type::Input) {
        isoID = naja::DNL::get()->getDNLTerminalFromID(node->data.termid).getIsoID();
        if (iso2boolExpr_.find(isoID) != iso2boolExpr_.end() && isoID != naja::DNL::DNLID_MAX) {
          setMemoETS(id, iso2boolExpr_[isoID]);
        }
      }
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
        SNLTruthTableTree::Node* const parent = node->tree->nodeFromId(node->parentIds[0]).get();
        assert(parent && parent->type == SNLTruthTableTree::Node::Type::P);
        assert(parent->data.termid < varNames.size());
        const auto& name = varNames[parent->data.termid];
        if (name == (size_t)-1) {
          // LCOV_EXCL_START
          throw std::runtime_error("Input variable index is SIZE_MAX");
          // LCOV_EXCL_STOP
        }
        if (name == 0) {
           setMemoETS(id, BoolExpr::createFalse());
           iso2boolExpr_[isoID] = BoolExpr::createFalse();
        } else if (name == 1) {
           setMemoETS(id, BoolExpr::createTrue());
           iso2boolExpr_[isoID] = BoolExpr::createTrue();
        } else {
          setMemoETS(id, BoolExpr::Var(name));
          iso2boolExpr_[isoID] = BoolExpr::Var(name);
        }
      }
    } else {
      // post-visit for Table / P
      const SNLTruthTable& tbl = node->getTruthTable();
      uint32_t k = tbl.size();
      uint64_t rows = uint64_t{1} << k;

      assert(!tbl.all0()); // Should be a PI
      assert(!tbl.all1()); // Should be a PI
      // if (tbl.all0()) {
      //   setMemoETS(id, BoolExpr::createFalse());
      //   iso2boolExpr_[isoID] = BoolExpr::createFalse();
      // } else if (tbl.all1()) {
      //   setMemoETS(id, BoolExpr::createTrue());
      //   iso2boolExpr_[isoID] = BoolExpr::createTrue();
      // } else 
      {
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

        size_t numRelIdx = 0;
        for (uint32_t j = 0; j < k; ++j) { if (getRelevantETS(j)) numRelIdx++; }

        // if nothing matters, fall back to constant-false
        // if (numRelIdx == 0) {
        //   setMemoETS(id, BoolExpr::createFalse());
        //   iso2boolExpr_[isoID] = BoolExpr::createFalse();
        // } else 
        assert(numRelIdx > 0 && "No relevant inputs for node");
        {
          // build the DNF terms
          clearTermsETS();
          for (uint64_t m = 0; m < rows; ++m) {
            if (!tbl.bits().bit(m)) continue;
            BoolExpr* term = nullptr;
            bool firstLit = true;
            BoolExpr* lit = nullptr;
            //for (uint32_t j : relIdx) {
            for (uint32_t j = 0; j < k; ++j) { 
              if (!getRelevantETS(j)) {
                continue;
              }
              bool bit1 = ((m >> j) & 1) != 0;
              lit = bit1 ? getChildFETS(j) : BoolExpr::Not(getChildFETS(j));
              if (firstLit) { term = lit; firstLit = false; }
              else { assert(term != nullptr); assert(lit != nullptr); term = BoolExpr::And(term, lit); }
            }
            // only push if we actually got a literal
            if (term) { pushBackTermsETS(std::move(term)); }
          }

          // guard against an empty terms list
          assert(!emptyTermsETS()); // Should be a PI
          // if (emptyTermsETS()) { 
          //   setMemoETS(id, BoolExpr::createFalse());
          //   iso2boolExpr_[isoID] = BoolExpr::createFalse();
          // }
          // else 
          {
            // fold into OR
            BoolExpr* expr = getTErmsETS().first[0];
            for (size_t t = 1; t < sizeOfTermsETS(); ++t) {
              expr = BoolExpr::Or(expr, getTErmsETS().first[t]);
            }
            setMemoETS(id, expr);
            iso2boolExpr_[isoID] = expr;
          }
        }
      }
    }
  }

  // 4) return root
  return getMemoETS(root->nodeID);
}
