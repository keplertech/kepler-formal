// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include "../../config/Config.h"
#include "kinduction/KInductionProblem.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <list>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace KEPLER_FORMAL::SEC {

enum class PDRStatus {
  Equivalent,
  Different,
  Inconclusive,
};

struct PDRResult {
  PDRStatus status = PDRStatus::Inconclusive;
  size_t bound = 0;
};

struct PDRQueryLimits {
  unsigned predecessorConflictLimit = 0;
  unsigned predecessorDecisionLimit = 0;
  // Blocking is the mandatory relative-induction query in IC3/PDR. A caller
  // may give it a deeper allowance while keeping optional generalization and
  // propagation queries inexpensive. Two-field aggregate initializers retain
  // one uniform limit for every query role.
  unsigned blockingConflictLimit = predecessorConflictLimit;
  unsigned blockingDecisionLimit = predecessorDecisionLimit;
  // Broad probes can bound transition-size estimation and return UNKNOWN so
  // the caller splits the exact property. Zero retains the engine defaults.
  size_t predecessorEncodingNodeLimit = 0;
  size_t predecessorNodeHintTargetLimit = 0;
  // Negative retains the deterministic default. This cumulative limit spans
  // optional invariant certification across every output batch sharing one
  // SEC model; it never limits an IC3/PDR property query.
  int64_t invariantCertificationTotalTickLimit = -1;
};

// Output batches from one SEC problem have the same immutable transition and
// startup model. This serial, scoped cache keeps exact model preparation, F[0]
// work, and certified inductive clauses across PDR runs. Property-specific
// bad-state searches and proof obligations remain local to each engine.
class PDRExactInitCache {
 public:
  struct Impl;

  PDRExactInitCache(
      const KInductionProblem& sourceProblem,
      KEPLER_FORMAL::Config::SolverType solverType);
  ~PDRExactInitCache();

  PDRExactInitCache(const PDRExactInitCache&) = delete;
  PDRExactInitCache& operator=(const PDRExactInitCache&) = delete;

 private:
  std::unique_ptr<Impl> impl_;

  friend class PDREngine;
};

namespace detail {

std::vector<size_t> makeDeterministicPdrWorklist(
    const std::unordered_set<size_t>& symbols);

bool pdrCubeLiteralOrderLess(size_t lhsSymbol,
                             bool lhsValue,
                             size_t rhsSymbol,
                             bool rhsValue);

bool pdrCubeAssignmentOrderLess(
    const std::vector<std::pair<size_t, bool>>& lhs,
    const std::vector<std::pair<size_t, bool>>& rhs);

bool pdrProofObligationPriorityLess(size_t lhsLevel,
                                    size_t lhsSequence,
                                    size_t rhsLevel,
                                    size_t rhsSequence);

inline std::vector<size_t> mergeSortedPdrSymbolVectors( // LCOV_EXCL_LINE
    const std::vector<size_t>& lhs,
    const std::vector<size_t>& rhs) {
  std::vector<size_t> merged; // LCOV_EXCL_LINE
  merged.reserve(lhs.size() + rhs.size()); // LCOV_EXCL_LINE
  std::set_union( // LCOV_EXCL_LINE
      lhs.begin(), // LCOV_EXCL_LINE
      lhs.end(), // LCOV_EXCL_LINE
      rhs.begin(), // LCOV_EXCL_LINE
      rhs.end(), // LCOV_EXCL_LINE
      std::back_inserter(merged)); // LCOV_EXCL_LINE
  return merged; // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

inline bool widenSortedPdrSymbolSurface( // LCOV_EXCL_LINE
    std::vector<size_t>& stableSurface,
    const std::vector<size_t>& requestedSurface) {
  if (std::includes( // LCOV_EXCL_LINE
          stableSurface.begin(), // LCOV_EXCL_LINE
          stableSurface.end(), // LCOV_EXCL_LINE
          requestedSurface.begin(), // LCOV_EXCL_LINE
          requestedSurface.end())) { // LCOV_EXCL_LINE
    return false; // LCOV_EXCL_LINE
  }
  // Keep the widened surface sorted and unique so it can be reused directly as
  // a FrameVariableStore symbol list and as part of the cache key.
  stableSurface = // LCOV_EXCL_LINE
      mergeSortedPdrSymbolVectors(stableSurface, requestedSurface); // LCOV_EXCL_LINE
  return true; // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

class PdrFrameSymbolSurfaceCache {
 public:
  const std::vector<size_t>& widen(
      const void* modelIdentity,
      size_t level,
      const std::vector<size_t>& requestedSurface,
      bool* widened = nullptr) {
    if (modelIdentity_ != modelIdentity) {
      // Symbol numbers belong to one transition model. A new model must not
      // inherit a surface from the previous PDR run.
      surfacesByLevel_.clear();
      modelIdentity_ = modelIdentity;
    }
    auto& stableSurface = surfacesByLevel_[level];
    const bool surfaceWidened =
        widenSortedPdrSymbolSurface(stableSurface, requestedSurface);
    if (widened != nullptr) {
      *widened = surfaceWidened;
    }
    return stableSurface;
  }

 private:
  const void* modelIdentity_ = nullptr;
  // IC3 keeps a distinct incremental SAT context for every frame. Preserve
  // each context's monotonic symbol surface when queries move between frames.
  std::unordered_map<size_t, std::vector<size_t>> surfacesByLevel_;
};

template <typename Key, typename Value, typename Hash>
class PdrWeightedLruCache {
 public:
  struct InsertResult {
    Value* value = nullptr;
    size_t evictedEntries = 0;
  };

  explicit PdrWeightedLruCache(size_t maxWeight) : maxWeight_(maxWeight) {}

  PdrWeightedLruCache(const PdrWeightedLruCache&) = delete;
  PdrWeightedLruCache& operator=(const PdrWeightedLruCache&) = delete;

  Value* find(const Key& key) {
    const auto existing = entries_.find(key);
    if (existing == entries_.end()) {
      return nullptr;
    }
    recency_.splice(recency_.begin(), recency_, existing->second.recency);
    return &existing->second.value;
  }

  InsertResult insert(Key key, Value value, size_t weight) {
    if (weight > maxWeight_) {
      return {};
    }
    if (Value* existing = find(key); existing != nullptr) {
      return {existing, 0};
    }

    size_t evictedEntries = 0;
    while (!entries_.empty() && weight > maxWeight_ - retainedWeight_) {
      evictLeastRecent();
      ++evictedEntries;
    }

    auto [inserted, insertedNew] = entries_.emplace(
        std::move(key), Entry{std::move(value), weight, {}});
    (void)insertedNew;
    recency_.push_front(&inserted->first);
    inserted->second.recency = recency_.begin();
    retainedWeight_ += weight;
    return {&inserted->second.value, evictedEntries};
  }

  size_t size() const { return entries_.size(); }
  size_t retainedWeight() const { return retainedWeight_; }

 private:
  struct Entry {
    Value value;
    size_t weight = 0;
    typename std::list<const Key*>::iterator recency;
  };

  void evictLeastRecent() {
    const Key* key = recency_.back();
    const auto existing = entries_.find(*key);
    retainedWeight_ -= existing->second.weight;
    recency_.pop_back();
    entries_.erase(existing);
  }

  size_t maxWeight_ = 0;
  size_t retainedWeight_ = 0;
  // Unordered-map references survive rehashing, so recency nodes can point at
  // map keys without storing a second ASIC-sized cube.
  std::list<const Key*> recency_;
  std::unordered_map<Key, Entry, Hash> entries_;
};

inline bool shouldResetPdrStableUnsatCache(size_t stableUnsatEntries,
                                           size_t maxEntries) {
  // Exact query results use a separate LRU. Stable UNSAT facts remain valid
  // across frame strengthening and therefore expire only at their own bound.
  return stableUnsatEntries >= maxEntries;
}

inline bool shouldSharePredecessorUnsatCore( // LCOV_EXCL_LINE
    size_t frameFingerprint,
    bool excludeTargetOnCurrentFrame) {
  // A predecessor core is reusable for stronger target cubes only in the base
  // PDR context. Do not share proofs that depended on the Q2 cube-exclusion
  // selector.
  return frameFingerprint == 0 && // LCOV_EXCL_LINE
         !excludeTargetOnCurrentFrame; // LCOV_EXCL_LINE
}

}  // namespace detail

// Top-level clause-based Property Directed Reachability strategy for SEC. It
// follows the classic proof-obligation/blocking loop over the already-built
// SEC transition system.
class PDREngine {
 public:
  PDREngine(const KInductionProblem& problem,
            KEPLER_FORMAL::Config::SolverType solverType,
            size_t maxPredecessorQueries = 0,
            std::shared_ptr<PDRExactInitCache> exactInitCache = nullptr);

  PDRResult run(size_t maxFrames) const;
  // Run the same transition system against an alternate safety property.
  // The target is deliberately separate from the model so it cannot alter
  // exact F[0]; the corresponding bad predicate is derived internally.
  PDRResult run(size_t maxFrames, BoolExpr* property) const;
  PDRResult run(size_t maxFrames,
                BoolExpr* property,
                const PDRQueryLimits& queryLimits) const;

 private:
  PDRResult runWithQueryLimits(size_t maxFrames,
                               BoolExpr* property,
                               const PDRQueryLimits* queryLimits) const;

  const KInductionProblem& problem_;
  KEPLER_FORMAL::Config::SolverType solverType_;
  size_t maxPredecessorQueries_ = 0;
  std::shared_ptr<PDRExactInitCache> exactInitCache_;
};

}  // namespace KEPLER_FORMAL::SEC
