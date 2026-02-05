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

// Temporary storage for DNF terms (as BoolExpr*) while building expressions for a node.
// Used to accumulate DNF terms before folding them into an OR expression.
typedef std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>> TempDNFTermsVec;
// Thread-local storage for DNF terms during node processing.
thread_local TempDNFTermsVec tempDNFTerms;
// Accessor for thread-local DNF terms vector.
TempDNFTermsVec& getTempDNFTerms() {
  return tempDNFTerms;
}
// Returns the number of DNF terms currently stored.
size_t sizeOfTempDNFTerms() {
  return getTempDNFTerms().size();
}
// Clears the DNF terms vector.
void clearTempDNFTerms() {
  auto& local = getTempDNFTerms();
  local.clear();
}
// Appends a DNF term (BoolExpr*) to the vector.
void pushTempDNFTerm(BoolExpr* term) {
  getTempDNFTerms().emplace_back(term);
}
// Returns true if no DNF terms are present.
bool emptyTempDNFTerms() {
  return getTempDNFTerms().empty();
}

// Bitset per node indicating which child inputs actually affect the table (relevant inputs).
// Used to optimize DNF generation by ignoring irrelevant inputs.
typedef std::pair<std::vector<uint8_t, tbb::tbb_allocator<uint8_t>>, size_t> RelevantInputsPair;
// Thread-local storage for relevant input bits and logical size.
thread_local RelevantInputsPair relevantInputs;
// Accessor for thread-local relevant input bits.
RelevantInputsPair& getRelevantInputs() {
  return relevantInputs;
}
// Returns the current logical size of relevant inputs.
size_t sizeOfRelevantInputs() {
  return getRelevantInputs().second;
}
// Clears the relevant inputs vector.
void clearRelevantInputs() {
  auto& local = getRelevantInputs();
  local.first.clear();
}
// Sets the relevant flag for input i.
void setRelevantInput(size_t i, bool b) {
  auto& local = getRelevantInputs();
  if (i >= local.second) {
    // LCOV_EXCL_START
    assert(false && "setRelevantInput: index out of range");
    // LCOV_EXCL_STOP
  }
  local.first[i] = b;
}
// Gets the relevant flag for input i.
bool getRelevantInput(size_t i) {
  auto& local = getRelevantInputs();
  if (i >= local.second) {
    // LCOV_EXCL_START
    throw std::out_of_range("getRelevantInput: index out of range");
    // LCOV_EXCL_STOP
  }
  return local.first[i];
}
// Ensure the relevant inputs vector has at least n entries and initialize them to false.
void reserveRelevantInputsWithFalse(size_t n) {
  auto& local = getRelevantInputs();
  auto& vec = local.first;
  auto& sz = local.second;
  if (vec.size() >= n) {
    vec.assign(n, false);
    sz = n;
    return;
  }
  vec.resize(n, false);
  vec.assign(n, false);
  sz = n;
}

// Memoization table (thread-local) mapping node IDs to BoolExpr* results.
// Used to cache results of subtrees during traversal.
typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> NodeExprCachePair;
thread_local NodeExprCachePair nodeExprCache;
// Accessor for thread-local node expression cache.
NodeExprCachePair& getNodeExprCache() {
  return nodeExprCache;
}
// Clears the node expression cache.
void clearNodeExprCache() {
  auto& local = getNodeExprCache();
  local.first.clear();
}
// Reserve and initialize the node expression cache to a given logical size.
void reserveNodeExprCache(size_t n) {
  auto& local = getNodeExprCache();
  auto& vec = local.first;
  auto& sz = local.second;
  if (vec.size() >= n) {
    sz = n;
    vec.assign(n, nullptr);
    return;
  }
  vec.resize(n);
  sz = n;
  vec.assign(n, nullptr);
}
// Sets the cached BoolExpr* for node i.
void setNodeExprCache(size_t i, BoolExpr* expr) {
  auto& local = getNodeExprCache();
  assert(i < local.second && "setNodeExprCache: index out of range");
  local.first[i] = expr;
}
// Gets the cached BoolExpr* for node i.
BoolExpr* getNodeExprCache(size_t i) {
  auto& local = getNodeExprCache();
  assert(i < local.second && "getNodeExprCache: index out of range");
  return local.first[i];
}

// Temporary storage for child BoolExpr pointers while processing a table node.
// Used to gather child expressions of the current node for DNF construction.
typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> ChildExprsPair;
thread_local ChildExprsPair childExprs;
// Accessor for thread-local child expressions.
ChildExprsPair& getChildExprs() {
  return childExprs;
}
// Clears the child expressions (sets logical size to 0).
void clearChildExprs() {
  auto& local = getChildExprs();
  local.second = 0;
}
// Reserve child expression storage and initialize entries to nullptr.
void reserveChildExprs(size_t n) {
  auto& local = getChildExprs();
  auto& vec = local.first;
  auto& sz = local.second;
  if (vec.size() >= n) {
    sz = n;
    vec.assign(n, nullptr);
    return;
  }
  vec.resize(n);
  sz = n;
  vec.assign(n, nullptr);
}
// Gets the child expression at index i.
BoolExpr* getChildExpr(size_t i) {
  auto& local = getChildExprs();
  assert(i < local.second && "getChildExpr: index out of range");
  return local.first[i];
}
// Sets the child expression at index i.
void setChildExpr(size_t i, BoolExpr* expr) {
  auto& local = getChildExprs();
  assert(i < local.second && "setChildExpr: index out of range");
  local.first[i] = expr;
}

// Frame type used for explicit stack-based post-order traversal.
// Each frame holds a pointer to a node and a boolean indicating whether
// the node has been visited (post-visit) or not (pre-visit).
using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
// Explicit traversal stack used to avoid recursion during tree traversal.
thread_local std::vector<Frame, tbb::tbb_allocator<Frame>> traversalStack;
// Accessor for thread-local traversal stack.
std::vector<Frame, tbb::tbb_allocator<Frame>>& getTraversalStack() {
  return traversalStack;
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

  // 2) node expression cache: clear and reserve memoization storage for all node IDs.
  clearNodeExprCache();
  reserveNodeExprCache(maxID + 1);

  // 3) post-order build using an explicit stack to avoid recursion.
  auto & stack = getTraversalStack();
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
          setNodeExprCache(id, iso2boolExpr_[isoID]);
        }
      }
      // If memo already contains an expression for this node, skip processing.
      if (getNodeExprCache(id) != nullptr) continue;

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
           setNodeExprCache(id, BoolExpr::createFalse());
           iso2boolExpr_[isoID] = BoolExpr::createFalse();
        } else if (name == 1) {
           setNodeExprCache(id, BoolExpr::createTrue());
           iso2boolExpr_[isoID] = BoolExpr::createTrue();
        } else {
          // Normal variable mapping.
          setNodeExprCache(id, BoolExpr::Var(name));
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
        clearChildExprs();
        reserveChildExprs(k);
        for (uint32_t i = 0; i < k; ++i) {
          size_t cid = node->tree->nodeFromId(node->childrenIds[i])->nodeID;
          setChildExpr(i, getNodeExprCache(cid));
        }

        // Determine which inputs actually matter for this truth table.
        // For each input j, check if flipping bit j changes the table output.
        clearRelevantInputs();
        reserveRelevantInputsWithFalse(k);
        for (uint32_t j = 0; j < k; ++j) {
          for (uint64_t m = 0; m < rows; ++m) {
            bool b0 = tbl.bits().bit(m);
            bool b1 = tbl.bits().bit(m ^ (uint64_t{1} << j));
            if (b0 != b1) { setRelevantInput(j, true); break; }
          }
        }

        // Count how many inputs are relevant.
        size_t numRelIdx = 0;
        for (uint32_t j = 0; j < k; ++j) { if (getRelevantInput(j)) numRelIdx++; }

        // The algorithm expects at least one relevant input for a PI node.
        assert(numRelIdx > 0 && "No relevant inputs for node");
        {
          // DNF (Disjunctive Normal Form) construction:
          // For each row in the truth table where the output is 1, create a term that is
          // an AND of literals corresponding to the relevant inputs (child expressions
          // or their negations depending on the row). Then these terms will be OR-ed
          // together to form the full expression for this node.
          // Build DNF terms by iterating over rows where the table output is 1.
          // For each such row, create a conjunction of literals for relevant inputs.
          clearTempDNFTerms();
          for (uint64_t m = 0; m < rows; ++m) {
            if (!tbl.bits().bit(m)) continue;
            BoolExpr* term = nullptr;
            bool firstLit = true;
            BoolExpr* lit = nullptr;
            // For each relevant input, pick the literal (child or its negation)
            // according to the bit value in row m.
            for (uint32_t j = 0; j < k; ++j) { 
              if (!getRelevantInput(j)) {
                continue;
              }
              bool bit1 = ((m >> j) & 1) != 0;
              lit = bit1 ? getChildExpr(j) : BoolExpr::Not(getChildExpr(j));
              if (firstLit) { term = lit; firstLit = false; }
              else { assert(term != nullptr); assert(lit != nullptr); term = BoolExpr::And(term, lit); }
            }
            // Only push the term if at least one literal was included.
            if (term) { pushTempDNFTerm(std::move(term)); }
          }

          // There must be at least one term for a PI node.
          assert(!emptyTempDNFTerms()); // Should be a PI

          {
            // Fold the list of terms into a single expression by OR-ing them.
            BoolExpr* expr = getTempDNFTerms()[0];
            DEBUG_LOG("number of rows for node ID %zu: %zu\n", id, sizeOfTempDNFTerms());
            for (size_t t = 1; t < sizeOfTempDNFTerms(); ++t) {
              expr = BoolExpr::Or(expr, getTempDNFTerms()[t]);
              DEBUG_LOG("Intermediate OR expr for node ID %zu: %s\n", id, expr->toString().c_str());
            }
            // Store the resulting expression in the node expression cache and in the iso map.
            setNodeExprCache(id, expr);
            iso2boolExpr_[isoID] = expr;
            DEBUG_LOG("Bool expression for node ID %zu: %s\n", id, expr->toString().c_str());
          }
        }
      }
    }
  }

  // 4) return root expression from node expression cache.
  return getNodeExprCache(root->nodeID);
}
