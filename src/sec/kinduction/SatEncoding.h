// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <set>
#include <unordered_map>
#include <vector>

#include "BoolExpr.h"
#include "../../sat/SATSolverWrapper.h"

namespace KEPLER_FORMAL::SEC {

// Owns the SAT literals that represent each symbolic SEC variable in each time
// frame of the unrolled problem.
class FrameVariableStore {
 public:
  FrameVariableStore(SATSolverWrapper& solver,
                     const std::vector<size_t>& symbols,
                     size_t numFrames);

  void addSymbols(SATSolverWrapper& solver,
                  const std::vector<size_t>& symbols);
  bool hasSymbol(size_t symbol) const;
  int getLiteral(size_t symbol, size_t frame) const;
  std::unordered_map<size_t, int> makeLeafLits(size_t frame) const;
  std::unordered_map<size_t, int> makeLeafLits(
      size_t frame,
      const std::vector<size_t>& symbols) const;
  std::unordered_map<size_t, int> makeLeafLits(
      size_t frame,
      const std::set<size_t>& symbols) const;

 private:
  std::unordered_map<size_t, std::vector<int>> symbolFrameLits_;
  size_t numFrames_ = 0;
};

// Converts a BoolExpr DAG into SAT clauses over one specific frame using a
// Tseitin-style encoding.
class FrameFormulaEncoder {
 public:
  FrameFormulaEncoder(SATSolverWrapper& solver,
                      std::unordered_map<size_t, int> leafLits);
  FrameFormulaEncoder(SATSolverWrapper& solver,
                      std::unordered_map<size_t, int> leafLits,
                      size_t expectedNodeHint);
  FrameFormulaEncoder(SATSolverWrapper& solver,
                      std::unordered_map<size_t, int> leafLits,
                      bool createMissingLeaves);
  FrameFormulaEncoder(SATSolverWrapper& solver,
                      std::unordered_map<size_t, int> leafLits,
                      bool createMissingLeaves,
                      size_t expectedNodeHint);
  FrameFormulaEncoder(SATSolverWrapper& solver,
                      std::unordered_map<size_t, int> leafLits,
                      const std::unordered_map<size_t, size_t>* symbolMap,
                      bool createMissingLeaves,
                      size_t expectedNodeHint);

  int encode(BoolExpr* expr);
  int encode(BoolExpr* expr, const std::vector<BoolExpr*>& postorder);
  void addLeafLiteral(size_t symbol, int literal);
  const std::unordered_map<size_t, int>& leafLits() const;

 private:
  struct BoolExprPtrHash {
    size_t operator()(const BoolExpr* node) const noexcept;
  };

  struct BoolExprPtrEqual {
    bool operator()(const BoolExpr* lhs, const BoolExpr* rhs) const noexcept;
  };

  struct CachedNodeLit {
    BoolExpr* node = nullptr;
    int lit = 0;
  };

  size_t mappedSymbol(size_t symbol) const;
  void reserveNodeCache();
  static size_t nodeCacheBucketCountFor(size_t desiredEntries);
  static size_t nodeCacheSlotFor(BoolExpr* node, size_t mask);
  void reserveNodeCacheSlots(size_t desiredEntries);
  void insertCachedLiteral(BoolExpr* node, int lit);
  int findCachedLiteral(BoolExpr* node) const;
  int cachedLiteral(BoolExpr* node) const;
  void cacheEncodedLiteral(BoolExpr* node, int lit);
  int getConstLit(bool value);
  bool isConstLit(int lit, bool value);
  void encodeReadyNode(BoolExpr* node);

  SATSolverWrapper& solver_;
  std::unordered_map<size_t, int> leafLits_;
  const std::unordered_map<size_t, size_t>* symbolMap_ = nullptr;
  bool createMissingLeaves_ = false;
  size_t expectedNodeHint_ = 0;
  // SEC creates many short-lived encoders while asking local SAT queries.
  // Keep this cache flat so each BoolExpr DAG node insertion is one table probe
  // instead of a per-node unordered_map allocation.
  std::vector<CachedNodeLit> nodeToLit_;
  // Hot encoder probes run millions of times in dual-rail transition CNF
  // construction. Keep raw accessors beside the owning vector so lookup avoids
  // repeated std::vector size/data calls while preserving vector ownership.
  CachedNodeLit* nodeToLitData_ = nullptr;
  const CachedNodeLit* nodeToLitConstData_ = nullptr;
  size_t nodeToLitMask_ = 0;
  size_t nodeToLitSize_ = 0;
  size_t nodeMapReservedEntries_ = 0;
  std::optional<int> trueLit_;
};

void addLiteralEquivalence(SATSolverWrapper& solver, int lhs, int rhs);
int createXorLiteral(SATSolverWrapper& solver, int lhs, int rhs);
// Enforces the noncyclic-path refinement from k-induction by requiring every
// pair of frames to differ in at least one state bit.
void addSimplePathConstraint(SATSolverWrapper& solver,
                             const FrameVariableStore& variables,
                             const std::vector<size_t>& stateSymbols,
                             size_t numFrames);
void addSimplePathConstraint(SATSolverWrapper& solver,
                             const FrameVariableStore& variables,
                             const std::vector<size_t>& stateSymbols,
                             size_t firstFrame,
                             size_t numFrames);

}  // namespace KEPLER_FORMAL::SEC
