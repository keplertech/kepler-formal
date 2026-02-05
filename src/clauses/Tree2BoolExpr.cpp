// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only
//
// Annotated version: comments added to explain the flow, data structures,
// and algorithmic steps. No executable code has been changed; only
// explanatory comments were inserted and commented-out code blocks
// that were previously disabled have been removed.

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

// Global concurrent map from DNL iso IDs to already-created BoolExpr pointers.
// This allows sharing of BoolExpr objects across conversions and threads.
tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*> Tree2BoolExpr::iso2boolExpr_ =
    tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*>();

// Helper typedefs for thread-local containers. Each pair stores a vector
// allocated with TBB allocator and a size_t representing the logical size.
// These thread-local containers avoid repeated allocations during traversal.

typedef std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>> TermsPair;

// Thread-local storage for DNF terms while building expressions for a node.
thread_local TermsPair terms;

TermsPair& getTErms() {
  return terms;
}

size_t sizeOfTerms() {
  return getTErms().size();
}

void clearTerms() {
  auto& termsLocal = getTErms();
  termsLocal.clear();
}

void pushBackTerms(BoolExpr* term) {
  getTErms().emplace_back(term);
}

bool emptyTerms() {
  return getTErms().empty();
}

// Relevant inputs bitset per node: which child inputs actually affect the table.
// Stored as a vector<uint8_t> (bool-like) with a logical size.
typedef std::pair<std::vector<uint8_t, tbb::tbb_allocator<uint8_t>>, size_t> RelevantPair;

thread_local RelevantPair relevant;

RelevantPair& getRelevant() {
  return relevant;
}

size_t sizeOfRelevant() {
  return getRelevant().second;
}

void clearRelevant() {
  auto& relevantLocal = getRelevant();
  relevantLocal.first.clear();
}

void setRelevant(size_t i, bool b) {
  auto& relevantLocal = getRelevant();
  if (i >= relevantLocal.second) {
    // LCOV_EXCL_START
    assert(false && "setRelevant: index out of range");
    // LCOV_EXCL_STOP
  }
  relevantLocal.first[i] = b;
}

bool getRelevant(size_t i) {
  auto& relevantLocal = getRelevant();
  if (i >= relevantLocal.second) {
    // LCOV_EXCL_START
    throw std::out_of_range("getRelevant: index out of range");
    // LCOV_EXCL_STOP
  }
  return relevantLocal.first[i];
}

// Ensure the relevant vector has at least n entries and initialize them to false.
// The logical size is stored in the second element of the pair.
void reserveRelevantwithFalse(size_t n) {
  auto& relevantLocal = getRelevant();
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

// Memoization table (thread-local) mapping node IDs to BoolExpr* results.
// The vector is preallocated and indexed by node ID for O(1) lookup.
typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> MemoPair;

thread_local MemoPair memo;

MemoPair& getMemo() {
  return memo;
}

void clearMemo() {
  auto& memoLocal = getMemo();
  memoLocal.first.clear();
}

// Reserve and initialize the memo table to a given logical size.
// The vector elements are set to nullptr to indicate "not computed".
void reserveMemo(size_t n) {
  auto& memoLocal = getMemo();
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

void setMemo(size_t i, BoolExpr* expr) {
  auto& memoLocal = getMemo();
  assert(i < memoLocal.second && "setMemo: index out of range");
  memoLocal.first[i] = expr;
}

BoolExpr* getMemo(size_t i) {
  auto& memoLocal = getMemo();
  assert(i < memoLocal.second && "getMemo: index out of range");
  return memoLocal.first[i];
}

// Temporary storage for child BoolExpr pointers while processing a table node.
typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> ChildFPair;

thread_local ChildFPair childF;

ChildFPair& getChildF() {
  return childF;
}

void clearChildF() {
  auto& childLocal = getChildF();
  childLocal.second = 0;
}

// Reserve child storage and initialize entries to nullptr.
void reserveChildF(size_t n) {
  auto& childLocal = getChildF();
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

BoolExpr* getChildF(size_t i) {
  auto& childLocal = getChildF();
  assert(i < childLocal.second && "getChildF: index out of range");
  return childLocal.first[i];
}

void setChildF(size_t i, BoolExpr* expr) {
  auto& childLocal = getChildF();
  assert(i < childLocal.second && "setChildF: index out of range");
  childLocal.first[i] = expr;
}

// Frame type used for explicit stack-based post-order traversal.
// Each frame holds a pointer to a node and a boolean indicating whether
// the node has been visited (post-visit) or not (pre-visit).
using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
thread_local std::vector<Frame, tbb::tbb_allocator<Frame>> stack;

std::vector<Frame, tbb::tbb_allocator<Frame>>& getStack() {
  return stack;
}

// Main conversion routine: converts an SNLTruthTableTree into a BoolExpr.
// The varNames vector maps SNL variable indices to desired variable IDs
// (or special markers like 0/1 for constants).
BoolExpr* Tree2BoolExpr::convert(
  const SNLTruthTableTree& tree, const std::vector<size_t>& varNames) {

  const auto root = tree.getRoot();
  if (!root) return nullptr;

  // Determine maximum node ID to size memoization structures.
  size_t maxID = tree.getMaxID();

  // 2) memo table: clear and reserve memoization storage for all node IDs.
  clearMemo();
  reserveMemo(maxID + 1);

  // 3) post-order build using an explicit stack to avoid recursion.
  auto & stack = getStack();
  stack.clear();
  stack.emplace_back(root.get(), false);

  while (!stack.empty()) {
    Frame f = stack.back();
    stack.pop_back();
    const SNLTruthTableTree::Node* node = f.first;

    // isoID is used to map DNL terminals to shared BoolExpr instances.
    naja::DNL::DNLID isoID = naja::DNL::DNLID_MAX;
    if (node->type != SNLTruthTableTree::Node::Type::Input) {
      isoID = naja::DNL::get()->getDNLTerminalFromID(node->data.termid).getIsoID();
    }

    bool visited = f.second;
    size_t id = node->nodeID;

    if (!visited) {
      // Pre-visit: attempt to reuse an existing BoolExpr from iso2boolExpr_
      // if the node corresponds to a DNL terminal that was already converted.
      naja::DNL::DNLID isoID = naja::DNL::DNLID_MAX;
      if (node->type != SNLTruthTableTree::Node::Type::Input) {
        isoID = naja::DNL::get()->getDNLTerminalFromID(node->data.termid).getIsoID();
        if (iso2boolExpr_.find(isoID) != iso2boolExpr_.end() && isoID != naja::DNL::DNLID_MAX) {
          setMemo(id, iso2boolExpr_[isoID]);
        }
      }
      // If memo already contains an expression for this node, skip processing.
      if (getMemo(id) != nullptr) continue;

      // If node is a Table or P node, push it back as visited and push children.
      if (node->type == SNLTruthTableTree::Node::Type::Table || node->type == SNLTruthTableTree::Node::Type::P) {
        stack.emplace_back(node, true);
        for (const auto& c : node->childrenIds) stack.emplace_back(node->tree->nodeFromId(c).get(), false);
      } else {
        // Input node handling: map input nodes to variables or constants.
        assert(node->type == SNLTruthTableTree::Node::Type::Input);
        if (node->parentIds.size() > 1) {
          #ifdef DEBUG_PRINTS
          // Debug logging for inputs with multiple parents (should be rare).
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
        // Special handling for constant mappings: 0 -> false, 1 -> true.
        if (name == 0) {
           setMemo(id, BoolExpr::createFalse());
           iso2boolExpr_[isoID] = BoolExpr::createFalse();
        } else if (name == 1) {
           setMemo(id, BoolExpr::createTrue());
           iso2boolExpr_[isoID] = BoolExpr::createTrue();
        } else {
          // Normal variable mapping.
          setMemo(id, BoolExpr::Var(name));
          iso2boolExpr_[isoID] = BoolExpr::Var(name);
        }
      }
    } else {
      // Post-visit: node is a Table or P node and all children have been processed.
      const SNLTruthTable& tbl = node->getTruthTable();
      DEBUG_LOG("Processing node ID %zu with table:\n%s\n", id, tbl.toString().c_str());
      uint32_t k = tbl.size();
      uint64_t rows = uint64_t{1} << k;
      DEBUG_LOG("Node ID %zu has %u inputs and %llu rows\n", id, k, rows);

      // The code expects the table to represent a prime implicant (not all-0 or all-1).
      assert(!tbl.all0()); // Should be a PI
      assert(!tbl.all1()); // Should be a PI

      {
        // Gather child BoolExpr pointers into a temporary array for quick access.
        clearChildF();
        reserveChildF(k);
        for (uint32_t i = 0; i < k; ++i) {
          size_t cid = node->tree->nodeFromId(node->childrenIds[i])->nodeID;
          setChildF(i, getMemo(cid));
        }

        // Determine which inputs actually matter for this truth table.
        // For each input j, check if flipping bit j changes the table output.
        clearRelevant();
        reserveRelevantwithFalse(k);
        for (uint32_t j = 0; j < k; ++j) {
          for (uint64_t m = 0; m < rows; ++m) {
            bool b0 = tbl.bits().bit(m);
            bool b1 = tbl.bits().bit(m ^ (uint64_t{1} << j));
            if (b0 != b1) { setRelevant(j, true); break; }
          }
        }

        // Count how many inputs are relevant.
        size_t numRelIdx = 0;
        for (uint32_t j = 0; j < k; ++j) { if (getRelevant(j)) numRelIdx++; }

        // The algorithm expects at least one relevant input for a PI node.
        assert(numRelIdx > 0 && "No relevant inputs for node");
        {
          // Build DNF terms by iterating over rows where the table output is 1.
          // For each such row, create a conjunction of literals for relevant inputs.
          clearTerms();
          for (uint64_t m = 0; m < rows; ++m) {
            if (!tbl.bits().bit(m)) continue;
            BoolExpr* term = nullptr;
            bool firstLit = true;
            BoolExpr* lit = nullptr;
            // For each relevant input, pick the literal (child or its negation)
            // according to the bit value in row m.
            for (uint32_t j = 0; j < k; ++j) { 
              if (!getRelevant(j)) {
                continue;
              }
              bool bit1 = ((m >> j) & 1) != 0;
              lit = bit1 ? getChildF(j) : BoolExpr::Not(getChildF(j));
              if (firstLit) { term = lit; firstLit = false; }
              else { assert(term != nullptr); assert(lit != nullptr); term = BoolExpr::And(term, lit); }
            }
            // Only push the term if at least one literal was included.
            if (term) { pushBackTerms(std::move(term)); }
          }

          // There must be at least one term for a PI node.
          assert(!emptyTerms()); // Should be a PI

          {
            // Fold the list of terms into a single expression by OR-ing them.
            BoolExpr* expr = getTErms()[0];
            DEBUG_LOG("number of rows for node ID %zu: %zu\n", id, sizeOfTerms());
            for (size_t t = 1; t < sizeOfTerms(); ++t) {
              expr = BoolExpr::Or(expr, getTErms()[t]);
              DEBUG_LOG("Intermediate OR expr for node ID %zu: %s\n", id, expr->toString().c_str());
            }
            // Store the resulting expression in the memo table and in the iso map.
            setMemo(id, expr);
            iso2boolExpr_[isoID] = expr;
            DEBUG_LOG("Bool expression for node ID %zu: %s\n", id, expr->toString().c_str());
          }
        }
      }
    }
  }

  // 4) return root expression from memo table.
  return getMemo(root->nodeID);
}
