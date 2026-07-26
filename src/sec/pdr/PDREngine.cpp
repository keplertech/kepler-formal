// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only
#include "pdr/PDREngine.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <list>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/BoolExprUtils.h"
#include "common/ProofProblemDebug.h"
#include "common/SecDiag.h"
#include "proof/ProofEngineShared.h"
#include "proof/TransitionExprResolver.h"
#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC {

namespace detail {

std::vector<size_t> makeDeterministicPdrWorklist(
    const std::unordered_set<size_t>& symbols) {
  std::vector<size_t> worklist(symbols.begin(), symbols.end());
  std::sort(worklist.begin(), worklist.end());
  return worklist;
}

bool pdrCubeLiteralOrderLess(size_t lhsSymbol,
                             bool lhsValue,
                             size_t rhsSymbol,
                             bool rhsValue) {
  if (lhsSymbol != rhsSymbol) {
    return lhsSymbol < rhsSymbol;
  }
  return lhsValue < rhsValue;
}

bool pdrCubeAssignmentOrderLess(
    const std::vector<std::pair<size_t, bool>>& lhs,
    const std::vector<std::pair<size_t, bool>>& rhs) {
  if (lhs.size() != rhs.size()) {
    return lhs.size() < rhs.size();
  }
  return std::lexicographical_compare(
      lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), [](const auto& a,
                                                         const auto& b) {
        return pdrCubeLiteralOrderLess(a.first, a.second, b.first, b.second);
      });
}

bool pdrProofObligationPriorityLess(size_t lhsLevel,
                                    size_t lhsSequence,
                                    size_t rhsLevel,
                                    size_t rhsSequence) {
  if (lhsLevel != rhsLevel) {
    return lhsLevel < rhsLevel;
  }
  // Figure 6 of the FMCAD'11 PDR paper uses stack order within one frame.
  return lhsSequence > rhsSequence;
}

}  // namespace detail

// Overall PDR algorithm:
// 1. Build Init from the SEC startup constraints and reuse any already
//    validated strengthening invariant when it is sound to do so.
// 2. Maintain frames F[0], F[1], ... where each frame stores clauses known to
//    hold for all states reachable within that many steps.
// 3. At each level, ask whether a bad state still survives the current frame.
// 4. If so, recursively search for predecessors until either Init is reached
//    (real counterexample) or the bad cube is blocked by a learned clause.
// 5. Generalize learned blocking clauses, add them to all earlier frames, and
//    then propagate them forward when the transition relation preserves them.
// 6. Stop once two adjacent frames converge, when a real bug is found, or when
//    the requested frame budget is exhausted.

namespace {

// The init-intersection fast path runs inside literal-dropping generalization.
// On ASICs the complemented-state table can be enormous while each cube is
// tiny, so scanning the full table per literal costs more than the SAT queries
// it was meant to avoid. Above this limit we skip only the cheap contradiction
// shortcut; the exact Init SAT query still decides intersection.
constexpr size_t kMaxComplementPairsForCheapInitCheck = 1024;
// Per-group node counts are reserve hints only. Skip expensive reserve sizing
// for very wide groups; the exact whole-query resource guard is counted
// independently and must never treat a missing hint as proof exhaustion.
constexpr size_t kMaxExactTransitionNodeCountHintTargets = 512;
// Full-state bad cubes require discovering the complete state support. If the
// formula walk exceeds this resource bound, PDR returns inconclusive.
constexpr size_t kMaxPreciseBadCubeSupportNodes = 262144;
constexpr unsigned kDefaultDualRailBadCubeConflictLimit = 20000;
constexpr unsigned kDefaultDualRailPredecessorConflictLimit = 250 * 1000;
// Incremental assumption solving counts this as a propagation budget, so it
// needs more room than the conflict cap for ordinary exact predecessor queries.
constexpr unsigned kDefaultDualRailPredecessorDecisionLimit =
    10 * 1000 * 1000;
// Blocking a proof obligation is the mandatory relative-induction query in
// Figure 6. Keep its role-specific floor equal to the measured exact-query
// default so an explicit lower user limit can still override both values.
constexpr unsigned kDefaultDualRailBlockingConflictLimit = 250 * 1000;
constexpr unsigned kDefaultDualRailBlockingDecisionLimit = 10 * 1000 * 1000;
// Figure 7 asks a local sequence of status-only Q2 queries while removing
// literals from one blocked cube. Give its narrow exact solver a small probe
// budget; UNKNOWN falls through to the existing persistent solver with the
// original full limits below.
constexpr unsigned kNarrowGeneralizationProbeConflictLimit = 10 * 1000;
constexpr unsigned kNarrowGeneralizationProbeDecisionLimit = 150 * 1000;
// Reusable-invariant certification is optional strengthening, not an IC3
// proof obligation. Bound both one query and the aggregate CaDiCaL work across
// output batches so candidate reuse cannot dominate the exact PDR checks.
constexpr int64_t kDualRailInvariantCertificationPerQueryTickLimit =
    100 * 1000 * 1000;
constexpr uint64_t kDualRailInvariantCertificationTotalTickLimit =
    100 * 1000 * 1000;
// Encoding guards are based only on the exact predecessor cone. Every output
// batch receives the same finite limits; enclosing design or port counts never
// select a different PDR problem or resource policy.
constexpr size_t kDefaultDualRailPredecessorEncodingNodeLimit =
    7500 * 1000;
constexpr size_t kDefaultDualRailPredecessorEncodingSupportLimit = 64 * 1024;
constexpr const char* kDualRailPredecessorConflictLimitEnv =
    "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_CONFLICT_LIMIT";
constexpr const char* kDualRailPredecessorDecisionLimitEnv =
    "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT";
constexpr size_t kMaxInitExcludedCubeGeneralizationAttempts = 2;
constexpr size_t kDefaultPdrStatsInterval = 1000;
constexpr size_t kInitialPdrStatsQueries = 20;
// Query-result caching is an accelerator only.  Keep it bounded so a long SEC
// run cannot trade the predecessor-encoding wall for unbounded retained cubes.
// Measured dual-rail leaves issue 70-80k exact predecessor queries. Retaining
// one complete leaf avoids clearing a still-hot cache near the end of PDR.
constexpr size_t kMaxPredecessorQueryResultCacheEntries = 256 * 1024;
constexpr size_t kMaxPredecessorUnsatCoresPerContext = 4096;
constexpr size_t kMaxPredecessorTargetSurfaceCacheBytes = 64 * 1024 * 1024;
// Clauses learned by prior PDR runs are only candidates until the invariant
// finder proves that a subset is initial and inductive. Bound this optional
// cache by literals so it cannot grow with every output batch.
constexpr size_t kMaxReusableInvariantCandidateLiterals = 256 * 1024;
// Higher-frame predecessor solvers absorb learned PDR frame clauses, so keep
// those caches bounded on giant dual-rail leaves. Exact F[0] and its transition
// relation are immutable across obligations and output batches; retaining that
// single solver avoids re-encoding the full startup frontier for every cube.
constexpr size_t kMaxDualRailPredecessorSolverCacheStateSymbols = 256 * 1024;
// Keep a small set of higher-frame SAT contexts per IC3 level. A single
// monotonically widened context eventually carries unrelated output cones;
// a bounded five-way cache retains nearby broad/strict passes while LRU
// eviction still caps ASIC-sized solver ownership.
constexpr size_t kMaxSharedPredecessorSolversPerLevel = 5;
// Retired property selectors leave satisfied clauses in a shared SAT solver.
// Rebuild the cache after a bounded window so unrelated guarded output leaves
// do not spend their query budget traversing an ever-growing inactive history.
// This discards cached SAT state only; PDR frames and every query stay exact.
constexpr size_t kMaxRetiredGuardedPredecessorContexts = 16;
// Higher-frame bad-cube solvers permanently absorb learned frame clauses.
// Keep those caches bounded on giant dual-rail leaves. Exact F[0] is immutable,
// however, and is shared across output batches without accumulating PDR frame
// clauses, so retaining that one solver avoids rebuilding the same startup
// frontier for every output.
constexpr size_t kMaxDualRailBadCubeSolverCacheStateSymbols =
    kMaxDualRailPredecessorSolverCacheStateSymbols;

enum class PredecessorQueryPurpose {
  BlockObligation,
  GeneralizeBlocker,
  LiftBlocker,
  PropagateClause,
};

bool predecessorQueryNeedsModel(PredecessorQueryPurpose purpose) {
  // Figure 6 requests EXTRACTMODEL only while recursively blocking an
  // obligation. Figure 7 generalization, blocker lifting, and Figure 9
  // propagation need only the SAT status of the same exact query.
  return purpose == PredecessorQueryPurpose::BlockObligation;
}

const char* predecessorQueryPurposeName(PredecessorQueryPurpose purpose) {
  switch (purpose) {
    case PredecessorQueryPurpose::BlockObligation:
      return "block";
    case PredecessorQueryPurpose::GeneralizeBlocker:
      return "generalize";
    case PredecessorQueryPurpose::LiftBlocker:
      return "lift";
    case PredecessorQueryPurpose::PropagateClause:
      return "propagate";
  }
  return "unknown";  // LCOV_EXCL_LINE
}
// FrameFormulaEncoder already makes a small generic Tseitin reservation, but
// sampled dual-rail PDR leaves still spent most time growing CaDiCaL variable
// vectors while streaming known-large transition cones. Reserve a larger,
// bounded chunk from PDR when we have the transition DAG estimate.
constexpr size_t kMinPdrTransitionSolverReserveNodes = 64 * 1024;
constexpr size_t kMaxPdrTransitionSolverReserveHint = 512 * 1024;
// Cubes represent a concrete bad/predecessor state, while clauses are the
// blocked generalization of such a state stored in a PDR frame.
struct CubeLiteral {  // LCOV_EXCL_LINE
  size_t symbol = 0;  // LCOV_EXCL_LINE
  bool value = false;  // LCOV_EXCL_LINE

  bool operator==(const CubeLiteral& other) const {
    return symbol == other.symbol && value == other.value;
  }
};

struct CubeLiteralHash {
  size_t operator()(const CubeLiteral& literal) const {
    return std::hash<size_t>()(
        (literal.symbol << 1) ^ (literal.value ? 1ULL : 0ULL));
  }
};

struct ClauseLiteral {
  size_t symbol = 0;
  bool positive = false;

  bool operator==(const ClauseLiteral& other) const {
    return symbol == other.symbol && positive == other.positive;
  }
};

using StateCube = std::vector<CubeLiteral>;
using StateClause = std::vector<ClauseLiteral>;

void mixHashValue(size_t& seed, size_t value) {
  seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

struct StateCubeHash {
  size_t operator()(const StateCube& cube) const {
    size_t seed = 0x9e3779b97f4a7c15ULL;
    for (const auto& literal : cube) {
      mixHashValue(seed, CubeLiteralHash{}(literal));
    }
    return seed;
  }
};

size_t cubeFingerprint(const StateCube& cube) {
  return StateCubeHash{}(cube);
}

struct StateClauseHash {
  size_t operator()(const StateClause& clause) const {
    size_t seed = 0x517cc1b727220a95ULL;
    for (const auto& literal : clause) {
      mixHashValue(seed, std::hash<size_t>()(literal.symbol));
      mixHashValue(seed, std::hash<bool>()(literal.positive));
    }
    return seed;
  }
// LCOV_EXCL_START
};

bool cubeLiteralLess(const CubeLiteral& lhs, const CubeLiteral& rhs) {
  return detail::pdrCubeLiteralOrderLess(
      lhs.symbol, lhs.value, rhs.symbol, rhs.value);
}

bool clauseLiteralLess(const ClauseLiteral& lhs, const ClauseLiteral& rhs) {
  if (lhs.symbol != rhs.symbol) {
    return lhs.symbol < rhs.symbol;
  }
  return lhs.positive < rhs.positive;
}

bool stateCubeLess(const StateCube& lhs, const StateCube& rhs) {
  if (lhs.size() != rhs.size()) {
    return lhs.size() < rhs.size();
  }
  return std::lexicographical_compare(
      lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), cubeLiteralLess);
}

bool stateClauseLess(const StateClause& lhs, const StateClause& rhs) {
  if (lhs.size() != rhs.size()) {
    return lhs.size() < rhs.size();
  }
  return std::lexicographical_compare(
      lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), clauseLiteralLess);
}

void sortStateCubesDeterministically(std::vector<StateCube>& cubes) {
  std::sort(cubes.begin(), cubes.end(), stateCubeLess);
}


struct FrameClauses {
  // F[i] stores clauses known to hold for all states reachable within i steps.
  std::vector<StateClause> clauses;
  // Cached SAT queries can keep old, subsumed clauses soundly; they only need
  // to see every newly learned clause.  Keep an append-only log so they can
  // synchronize by delta instead of rescanning ASIC-size frames after each
  // local refinement.
  std::vector<StateClause> addedClauseLog;
  // Fingerprints identify an exact frame context for SAT-result caches. Frame
  // clauses change only through addClauseToFrame(), which invalidates this
  // host-side memo without changing any PDR clause or query.
  mutable std::optional<std::pair<size_t, size_t>> clauseFingerprint;
};

size_t frameClausesFingerprint(const std::vector<FrameClauses>& frames,
                               size_t level) {
  if (level >= frames.size()) {
    return 0; // LCOV_EXCL_LINE
  }
  const auto& frame = frames[level];
  if (!frame.clauseFingerprint.has_value() ||
      frame.clauseFingerprint->first != level) {
    // Preserve the original hash operation order exactly. Only the completed
    // value is memoized, so cache keys remain byte-for-byte unchanged.
    size_t seed = std::hash<size_t>()(level);
    mixHashValue(seed, std::hash<size_t>()(frame.clauses.size()));
    for (const auto& clause : frame.clauses) {
      mixHashValue(seed, StateClauseHash{}(clause));
    }
    frame.clauseFingerprint = std::pair{level, seed};
  }
  return frame.clauseFingerprint->second;
}

enum class IndexedStateRelationKind {
  Equality,
  Complement,
  DualRailValidity,
};

class OrderedStateRelationIndex {
 public:
  OrderedStateRelationIndex() = default;

  explicit OrderedStateRelationIndex(
      const std::vector<std::pair<size_t, size_t>>& pairs)
      : orderedPairs_(pairs) {
    buildIndex();
  }

  explicit OrderedStateRelationIndex(
      const std::vector<DualRailSymbolPair>& railPairs) {
    orderedPairs_.reserve(railPairs.size());
    for (const auto& rails : railPairs) {
      orderedPairs_.emplace_back(rails.mayBeOne, rails.mayBeZero);
    }
    buildIndex();
  }

  void addPartnerClosure(std::unordered_set<size_t>& symbols) const {
    std::vector<size_t> componentIndices;
    componentIndices.reserve(symbols.size());
    for (const size_t symbol : symbols) {
      const auto indexIt = pairIndicesBySymbol_.find(symbol);
      if (indexIt == pairIndicesBySymbol_.end()) {
        continue;
      }
      componentIndices.push_back(
          componentIndexByPair_[indexIt->second.front()]);
    }
    std::sort(componentIndices.begin(), componentIndices.end());
    componentIndices.erase(
        std::unique(componentIndices.begin(), componentIndices.end()),
        componentIndices.end());

    // Relation graphs are immutable for a PDR run. Reuse their exact connected
    // components instead of rediscovering the same transitive partner closure
    // for every predecessor surface. The caller still sorts the final symbols,
    // so SAT variable and clause order is unchanged.
    for (const size_t componentIndex : componentIndices) {
      symbols.insert(componentSymbols_[componentIndex].begin(),
                     componentSymbols_[componentIndex].end());
    }
  }

  void addClauses(SATSolverWrapper& solver,
                  const FrameVariableStore& variables,
                  const std::vector<size_t>& solverSymbols,
                  size_t numFrames,
                  IndexedStateRelationKind kind) const {
    // Dense F[0] surfaces already contain most state symbols. Preserve the old
    // linear pair scan there; sparse output cones use the reverse index below.
    if (solverSymbols.size() >= orderedPairs_.size()) {
      for (size_t frame = 0; frame < numFrames; ++frame) {
        for (size_t pairIndex = 0; pairIndex < orderedPairs_.size();
             ++pairIndex) {
          addPairClause(solver, variables, pairIndex, frame, kind);
        }
      }
      return;
    }

    const std::vector<size_t> pairIndices = relevantPairIndices(solverSymbols);
    for (size_t frame = 0; frame < numFrames; ++frame) {
      for (const size_t pairIndex : pairIndices) {
        addPairClause(solver, variables, pairIndex, frame, kind);
      }
    }
  }

  void addClausesForAddedSymbols(
      SATSolverWrapper& solver,
      const FrameVariableStore& variables,
      const std::vector<size_t>& addedSymbols,
      size_t frame,
      IndexedStateRelationKind kind) const {
    // Every relation that becomes encodable after a monotonic symbol-surface
    // extension touches at least one newly added symbol. Existing pairs are
    // already permanent clauses in the incremental solver.
    for (const size_t pairIndex : relevantPairIndices(addedSymbols)) {
      addPairClause(solver, variables, pairIndex, frame, kind);
    }
  }

 private:
  void addPairClause(SATSolverWrapper& solver,
                     const FrameVariableStore& variables,
                     size_t pairIndex,
                     size_t frame,
                     IndexedStateRelationKind kind) const {
    const auto& [lhs, rhs] = orderedPairs_[pairIndex];
    if (!variables.hasSymbol(lhs) || !variables.hasSymbol(rhs)) {
      return;
    }
    const int lhsLiteral = variables.getLiteral(lhs, frame);
    const int rhsLiteral = variables.getLiteral(rhs, frame);
    switch (kind) {
      case IndexedStateRelationKind::Equality:
        addLiteralEquivalence(solver, lhsLiteral, rhsLiteral);
        break;
      case IndexedStateRelationKind::Complement:
        addLiteralEquivalence(solver, rhsLiteral, -lhsLiteral);
        break;
      case IndexedStateRelationKind::DualRailValidity:
        solver.addClause({lhsLiteral, rhsLiteral});
        break;
    }
  }

  void buildIndex() {
    pairIndicesBySymbol_.reserve(orderedPairs_.size() * 2);
    for (size_t index = 0; index < orderedPairs_.size(); ++index) {
      const auto& [lhs, rhs] = orderedPairs_[index];
      pairIndicesBySymbol_[lhs].push_back(index);
      pairIndicesBySymbol_[rhs].push_back(index);
    }
    buildComponents();
  }

  void buildComponents() {
    constexpr size_t kUnassigned = std::numeric_limits<size_t>::max();
    componentIndexByPair_.assign(orderedPairs_.size(), kUnassigned);
    for (size_t start = 0; start < orderedPairs_.size(); ++start) {
      if (componentIndexByPair_[start] != kUnassigned) {
        continue;
      }

      const size_t componentIndex = componentSymbols_.size();
      componentSymbols_.emplace_back();
      std::vector<size_t> pairWorklist{start};
      componentIndexByPair_[start] = componentIndex;
      for (size_t cursor = 0; cursor < pairWorklist.size(); ++cursor) {
        const auto& [lhs, rhs] = orderedPairs_[pairWorklist[cursor]];
        componentSymbols_.back().push_back(lhs);
        componentSymbols_.back().push_back(rhs);
        for (const size_t symbol : {lhs, rhs}) {
          for (const size_t pairIndex : pairIndicesBySymbol_.at(symbol)) {
            if (componentIndexByPair_[pairIndex] == kUnassigned) {
              componentIndexByPair_[pairIndex] = componentIndex;
              pairWorklist.push_back(pairIndex);
            }
          }
        }
      }
      auto& component = componentSymbols_.back();
      std::sort(component.begin(), component.end());
      component.erase(std::unique(component.begin(), component.end()),
                      component.end());
    }
  }

  std::vector<size_t>
  relevantPairIndices(const std::vector<size_t>& symbols) const {
    std::vector<size_t> indices;
    for (const size_t symbol : symbols) {
      const auto indexIt = pairIndicesBySymbol_.find(symbol);
      if (indexIt == pairIndicesBySymbol_.end()) {
        continue;
      }
      indices.insert(indices.end(), indexIt->second.begin(),
                     indexIt->second.end());
    }
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
  }

  std::vector<std::pair<size_t, size_t>> orderedPairs_;
  std::unordered_map<size_t, std::vector<size_t>> pairIndicesBySymbol_;
  std::vector<size_t> componentIndexByPair_;
  std::vector<std::vector<size_t>> componentSymbols_;
};

// All entries originate from explicit same-design model relations.  The index
// uses combined symbol IDs only; it never relates internal elements by name.
struct ComplementPartnerIndex {
  explicit ComplementPartnerIndex(const KInductionProblem& problem)
      : complemented0(problem.complementedStatePairs0),
        complemented1(problem.complementedStatePairs1),
        sameFrameEqualities0(problem.sameFrameStateEqualityPairs0),
        sameFrameEqualities1(problem.sameFrameStateEqualityPairs1),
        dualRailPairs(problem.dualRailStatePairs) {}

  void addComplementedPartnerClosure(
      std::unordered_set<size_t>& symbols) const {
    complemented0.addPartnerClosure(symbols);
    complemented1.addPartnerClosure(symbols);
  }

  void addSameFrameEqualityPartnerClosure(
      std::unordered_set<size_t>& symbols) const {
    sameFrameEqualities0.addPartnerClosure(symbols);
    sameFrameEqualities1.addPartnerClosure(symbols);
  }

  void addDualRailPartnerClosure(std::unordered_set<size_t>& symbols) const {
    dualRailPairs.addPartnerClosure(symbols);
  }

  void addClauses(SATSolverWrapper& solver,
                  const FrameVariableStore& variables,
                  const std::vector<size_t>& solverSymbols,
                  size_t numFrames) const {
    complemented0.addClauses(solver, variables, solverSymbols, numFrames,
                             IndexedStateRelationKind::Complement);
    complemented1.addClauses(solver, variables, solverSymbols, numFrames,
                             IndexedStateRelationKind::Complement);
    sameFrameEqualities0.addClauses(solver, variables, solverSymbols, numFrames,
                                    IndexedStateRelationKind::Equality);
    sameFrameEqualities1.addClauses(solver, variables, solverSymbols, numFrames,
                                    IndexedStateRelationKind::Equality);
    dualRailPairs.addClauses(solver, variables, solverSymbols, numFrames,
                             IndexedStateRelationKind::DualRailValidity);
  }

  void addClausesForAddedSymbols(
      SATSolverWrapper& solver,
      const FrameVariableStore& variables,
      const std::vector<size_t>& addedSymbols,
      size_t frame) const {
    complemented0.addClausesForAddedSymbols(
        solver, variables, addedSymbols, frame,
        IndexedStateRelationKind::Complement);
    complemented1.addClausesForAddedSymbols(
        solver, variables, addedSymbols, frame,
        IndexedStateRelationKind::Complement);
    sameFrameEqualities0.addClausesForAddedSymbols(
        solver, variables, addedSymbols, frame,
        IndexedStateRelationKind::Equality);
    sameFrameEqualities1.addClausesForAddedSymbols(
        solver, variables, addedSymbols, frame,
        IndexedStateRelationKind::Equality);
    dualRailPairs.addClausesForAddedSymbols(
        solver, variables, addedSymbols, frame,
        IndexedStateRelationKind::DualRailValidity);
  }

  OrderedStateRelationIndex complemented0;
  OrderedStateRelationIndex complemented1;
  OrderedStateRelationIndex sameFrameEqualities0;
  OrderedStateRelationIndex sameFrameEqualities1;
  OrderedStateRelationIndex dualRailPairs;
};

struct ProofObligation {
  // "cube is bad at level" requests either a predecessor or a blocking clause.
  StateCube cube;
  size_t level = 0;
  size_t badFrame = 0;
  size_t sequence = 0;
};

struct ProofObligationKey {
  size_t level = 0;
  size_t badFrame = 0;
  StateCube cube;

  bool operator==(const ProofObligationKey& other) const {
    return level == other.level &&
           badFrame == other.badFrame &&
           cube == other.cube;
  }
};

struct ProofObligationKeyHash {
  size_t operator()(const ProofObligationKey& key) const {
    size_t seed = std::hash<size_t>()(key.level);
    mixHashValue(seed, std::hash<size_t>()(key.badFrame));
    mixHashValue(seed, StateCubeHash{}(key.cube));
    return seed;
  }
};

struct SymbolPair {
  size_t first = 0;
  size_t second = 0;

  bool operator==(const SymbolPair& other) const {
    return first == other.first && second == other.second;
  }
};

struct SymbolPairHash {
  size_t operator()(const SymbolPair& pair) const {
    // Splitmix-style mixing keeps pair lookup cheap and avoids repeatedly
    // scanning thousands of extracted startup equalities during PDR seeding.
    size_t seed = pair.first + 0x9e3779b97f4a7c15ULL;
    seed ^= pair.second + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

class InitParityRelations {
 public:
  void ensureSymbol(size_t symbol) {
    if (parent_.find(symbol) == parent_.end()) {
      parent_.emplace(symbol, symbol);
      parityToParent_.emplace(symbol, false);
    }
  }

  void addEquality(size_t lhs, size_t rhs) { unite(lhs, rhs, false); }

  void addComplement(size_t lhs, size_t rhs) { unite(lhs, rhs, true); }

  std::optional<std::pair<size_t, bool>> findWithParity(size_t symbol) const {
    const auto parentIt = parent_.find(symbol);
    if (parentIt == parent_.end()) {
      return std::nullopt;
    }
    const size_t parent = parentIt->second;
    // LCOV_EXCL_START
    const bool parity = parityToParent_.at(symbol);
    // LCOV_EXCL_STOP
    if (parent == symbol) {
      return std::pair{symbol, false};
    }
    const auto parentRoot = findWithParity(parent);
    if (!parentRoot.has_value()) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    return std::pair{parentRoot->first, parity ^ parentRoot->second};
  }

 private:
  std::pair<size_t, bool> mutableFind(size_t symbol) {
    // LCOV_EXCL_START
    ensureSymbol(symbol);
    const size_t parent = parent_[symbol];
    const bool parity = parityToParent_[symbol];
    if (parent == symbol) {
    // LCOV_EXCL_STOP
      return {symbol, false};
    }
    const auto root = mutableFind(parent);  // LCOV_EXCL_LINE
    parent_[symbol] = root.first;  // LCOV_EXCL_LINE
    parityToParent_[symbol] = parity ^ root.second;  // LCOV_EXCL_LINE
    return {parent_[symbol], parityToParent_[symbol]};  // LCOV_EXCL_LINE
  }

  void unite(size_t lhs, size_t rhs, bool inverted) {
    const auto lhsRoot = mutableFind(lhs);
    const auto rhsRoot = mutableFind(rhs);
    if (lhsRoot.first == rhsRoot.first) {
      return; // LCOV_EXCL_LINE
    }
    parent_[lhsRoot.first] = rhsRoot.first;
    // value(lhs) xor value(rhs) must equal `inverted`.
    parityToParent_[lhsRoot.first] =
        lhsRoot.second ^ rhsRoot.second ^ inverted;
  }

  std::unordered_map<size_t, size_t> parent_;
  std::unordered_map<size_t, bool> parityToParent_;
};

struct ExprPair {
  BoolExpr* first = nullptr;
  BoolExpr* second = nullptr;

  bool operator==(const ExprPair& other) const {
    return first == other.first && second == other.second;
  }
};

struct ExprPairHash {
  size_t operator()(const ExprPair& pair) const {
    size_t seed =
        reinterpret_cast<size_t>(pair.first) + 0x9e3779b97f4a7c15ULL;
    seed ^= reinterpret_cast<size_t>(pair.second) +
            0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct InitFactIndex {
  std::unordered_map<size_t, bool> assignments;
  // Preserve every literal from the source assignment vector, including a
  // malformed duplicate with the opposite value. Bit 0 records false and bit
  // 1 records true, matching the old full-vector contradiction scan exactly.
  std::unordered_map<size_t, unsigned char> assignmentValueMasks;
  std::unordered_map<size_t, bool> rootAssignments;
  std::unordered_set<SymbolPair, SymbolPairHash> equalities;
  std::unordered_set<SymbolPair, SymbolPairHash> complements;
  InitParityRelations relations;
};

struct TransitionAssumptionKey {
  size_t transitionSymbol = 0;
  bool desiredValue = false;

  bool operator==(const TransitionAssumptionKey& other) const {
    return transitionSymbol == other.transitionSymbol &&
           desiredValue == other.desiredValue;
  }
};

struct TransitionAssumptionKeyHash {
  size_t operator()(const TransitionAssumptionKey& key) const {
    size_t seed = std::hash<size_t>()(key.transitionSymbol);
    mixHashValue(seed, std::hash<bool>()(key.desiredValue));
    return seed;
  }
};

struct PredecessorQueryResultKey { // LCOV_EXCL_LINE
  const KInductionProblem* problem = nullptr;
  const TransitionExprResolver* transitionByState = nullptr;
  const BoolExpr* initFormula = nullptr;
  const BoolExpr* frameInvariant = nullptr;
  size_t level = 0;
  size_t frameFingerprint = 0;
  bool excludeTargetOnCurrentFrame = false;
  StateCube targetCube;

  bool operator==(const PredecessorQueryResultKey& other) const {
    return problem == other.problem &&
           transitionByState == other.transitionByState &&
           initFormula == other.initFormula &&
           frameInvariant == other.frameInvariant &&
           level == other.level &&
           frameFingerprint == other.frameFingerprint &&
           excludeTargetOnCurrentFrame == other.excludeTargetOnCurrentFrame &&
           targetCube == other.targetCube;
  }
};

struct PredecessorQueryResultKeyHash {
  size_t operator()(const PredecessorQueryResultKey& key) const {
    size_t seed = std::hash<const void*>()(key.problem);
    mixHashValue(seed, std::hash<const void*>()(key.transitionByState));
    mixHashValue(seed, std::hash<const void*>()(key.initFormula));
    mixHashValue(seed, std::hash<const void*>()(key.frameInvariant));
    mixHashValue(seed, std::hash<size_t>()(key.level));
    mixHashValue(seed, std::hash<size_t>()(key.frameFingerprint));
    mixHashValue(seed, std::hash<bool>()(key.excludeTargetOnCurrentFrame));
    mixHashValue(seed, StateCubeHash{}(key.targetCube));
    return seed;
  }
};

struct PredecessorQueryResultEntry {
  bool hasPredecessor = false;
  bool hasPredecessorModel = false;
  StateCube predecessor;
  bool hasUnsatCore = false;
  StateCube unsatCore;
};

struct PredecessorUnsatCoreCacheKey {
  const KInductionProblem* problem = nullptr;
  const TransitionExprResolver* transitionByState = nullptr;
  const BoolExpr* initFormula = nullptr;
  const BoolExpr* frameInvariant = nullptr;
  size_t level = 0;
  bool excludeTargetOnCurrentFrame = false;

  bool operator==(const PredecessorUnsatCoreCacheKey& other) const {
    return problem == other.problem &&
           transitionByState == other.transitionByState &&
           initFormula == other.initFormula &&
           frameInvariant == other.frameInvariant &&
           level == other.level &&
           excludeTargetOnCurrentFrame == other.excludeTargetOnCurrentFrame;
  }
};

struct PredecessorUnsatCoreCacheKeyHash {
  size_t operator()(const PredecessorUnsatCoreCacheKey& key) const {
    size_t seed = std::hash<const void*>()(key.problem);
    mixHashValue(seed, std::hash<const void*>()(key.transitionByState));
    mixHashValue(seed, std::hash<const void*>()(key.initFormula));
    mixHashValue(seed, std::hash<const void*>()(key.frameInvariant));
    mixHashValue(seed, std::hash<size_t>()(key.level));
    mixHashValue(seed, std::hash<bool>()(key.excludeTargetOnCurrentFrame));
    return seed;
  }
};

struct PredecessorQueryResultStore {
  // Exact frame fingerprints age quickly as IC3 learns clauses. Keep the most
  // recently reused exact answers instead of clearing the complete cache when
  // its entry bound is reached.
  detail::PdrWeightedLruCache<PredecessorQueryResultKey,
                              PredecessorQueryResultEntry,
                              PredecessorQueryResultKeyHash>
      queryResults{kMaxPredecessorQueryResultCacheEntries};
  std::unordered_set<PredecessorQueryResultKey,
                     PredecessorQueryResultKeyHash>
      unsatQueries;
  std::unordered_map<PredecessorUnsatCoreCacheKey,
                     std::vector<StateCube>,
                     PredecessorUnsatCoreCacheKeyHash>
      unsatCoresByContext;
};

class PdrFormulaSupportCache;

struct PredecessorFrameSymbolSurfaceKey {
  const KInductionProblem* problem = nullptr;
  BoolExpr* initFormula = nullptr;
  BoolExpr* frameInvariant = nullptr;
  const ComplementPartnerIndex* complementPartners = nullptr;
  const PdrFormulaSupportCache* supportCache = nullptr;
  size_t level = 0;
  size_t frameFingerprint = 0;

  bool operator==(const PredecessorFrameSymbolSurfaceKey& other) const { // LCOV_EXCL_LINE
    return problem == other.problem && // LCOV_EXCL_LINE
           initFormula == other.initFormula && // LCOV_EXCL_LINE
           frameInvariant == other.frameInvariant && // LCOV_EXCL_LINE
           complementPartners == other.complementPartners && // LCOV_EXCL_LINE
           supportCache == other.supportCache && // LCOV_EXCL_LINE
           level == other.level && // LCOV_EXCL_LINE
           frameFingerprint == other.frameFingerprint; // LCOV_EXCL_LINE
  }
};

struct PredecessorFrameSymbolSurface {
  bool valid = false;
  PredecessorFrameSymbolSurfaceKey key;
  std::vector<size_t> symbols;
};

struct TransitionEncodingLiteral {
  size_t transitionSymbol = 0;
  bool desiredValue = false;
  CubeLiteral originalLiteral;
};

struct TransitionEncodingLiteralGroup {
  const std::unordered_map<size_t, size_t>* symbolMap = nullptr;
  std::vector<TransitionEncodingLiteral> literals;
  std::vector<size_t> stateSymbols;
};

struct PdrTernarySimulationRoot {
  BoolExpr* formula = nullptr;
  const std::unordered_map<size_t, size_t>* symbolMap = nullptr;
  bool expectedValue = false;
};

std::vector<TransitionEncodingLiteralGroup>
groupTransitionCubeLiteralsBySymbolMap(
    const TransitionExprResolver& transitionByState,
    const StateCube& targetCube);

struct PredecessorTargetSurface { // LCOV_EXCL_LINE
  std::vector<size_t> targetSymbols;
  std::vector<size_t> encodedTargets;
  std::vector<size_t> transitionSupportSymbols;
  std::vector<size_t> predecessorSymbols;
  std::vector<size_t> closedPredecessorSymbols;
  std::vector<size_t> closedTransitionSupportSymbols;
  StateClause exclusionClause;
  // Target cubes are immutable. Keep their symbol-map grouping beside the
  // support vectors so every exact SAT query reuses the same transition roots.
  std::vector<TransitionEncodingLiteralGroup> transitionGroups;
  std::vector<PdrTernarySimulationRoot> ternaryRoots;
  size_t transitionEncodingNodes = 0;
};

using PredecessorTargetSurfaceCache =
    detail::PdrWeightedLruCache<StateCube,
                                PredecessorTargetSurface,
                                StateCubeHash>;

struct PredecessorAssumptionCacheKey {
  const KInductionProblem* problem = nullptr;
  // This is an identity token only; transition expressions are always read
  // from the resolver passed to the current exact query.
  const void* transitionModel = nullptr;
  const BoolExpr* initFormula = nullptr;
  const BoolExpr* frameInvariant = nullptr;
  size_t level = 0;
  size_t frameFingerprint = 0;
  std::vector<size_t> solverSymbols;

  bool operator==(const PredecessorAssumptionCacheKey& other) const {
    return problem == other.problem &&
           transitionModel == other.transitionModel &&
           initFormula == other.initFormula &&
           frameInvariant == other.frameInvariant &&
           level == other.level &&
           frameFingerprint == other.frameFingerprint &&
           solverSymbols == other.solverSymbols;
  }

  bool hasSameReusableContext(
      const PredecessorAssumptionCacheKey& other) const {
    return problem == other.problem &&
           transitionModel == other.transitionModel &&
           initFormula == other.initFormula &&
           frameInvariant == other.frameInvariant &&
           level == other.level;
  }
};

struct PreparedPredecessorTargetAssumptions {
  std::vector<int> assumptions;
  std::unordered_map<int, CubeLiteral> literalByAssumption;
};

struct GuardedPredecessorFrameContext {
  size_t runId = 0;
  int activationLiteral = 0;
  const BoolExpr* property = nullptr;
  const BoolExpr* frameInvariant = nullptr;
  std::unordered_set<StateClause, StateClauseHash> emittedFrameClauses;
  size_t emittedFrameFingerprint = 0;
  size_t emittedFrameLogOffset = 0;
};

struct PredecessorAssumptionSolver {
  PredecessorAssumptionCacheKey key;
  std::unique_ptr<SATSolverWrapper> solver;
  std::unique_ptr<FrameVariableStore> variables;
  // The cached SAT model is useful only if predecessor extraction can read the
  // transition-expression leaves that were encoded under assumptions.
  std::unordered_map<size_t, int> transitionLeafLits;
  std::unordered_map<TransitionAssumptionKey, int, TransitionAssumptionKeyHash>
      assumptionByTransitionLiteral;
  // The paper changes only assumptions between predecessor solves. Cache the
  // exact ordered vector and failed-core lookup for each target encoded in this
  // solver; every hit still executes the same SAT query.
  std::unordered_map<StateClause,
                     PreparedPredecessorTargetAssumptions,
                     StateClauseHash>
      preparedTargetAssumptions;
  size_t preparedTargetAssumptionHits = 0;
  // Exclusion selectors are target-specific, so retain one scratch vector for
  // appending that selector without modifying the cached base assumptions.
  std::vector<int> targetAssumptions;
  // Reuse the transition-DAG encoder together with the cached predecessor
  // solver. Neighboring dual-rail PDR targets often share most of the same
  // transition cone; keeping the encoder node cache avoids re-emitting that
  // Tseitin structure for every target literal.
  std::unordered_map<const std::unordered_map<size_t, size_t>*,
                     std::unique_ptr<FrameFormulaEncoder>>
      transitionEncoderBySymbolMap;
  // Exact F[0] bad-state queries differ from predecessor queries only by
  // assumptions. Encode their roots lazily in the same incremental solver so
  // a whole-chip Init CNF is not retained a second time.
  std::unique_ptr<FrameFormulaEncoder> badRootEncoder;
  std::unordered_map<BoolExpr*, int> encodedBadRoots;
  // Figure 6 repeatedly asks whether candidate cubes intersect exact Init.
  // Keep those answers with the same F[0] solver that answers the query.
  std::unordered_map<StateCube, bool, StateCubeHash>
      initIntersectionResultByCube;
  std::unordered_set<size_t> querySymbolSet;
  std::unordered_set<StateClause, StateClauseHash> emittedFrameClauses;
  size_t emittedFrameFingerprint = 0;
  size_t emittedFrameLogOffset = 0;
  struct Q2SelectorCacheEntry {
    int selector = 0;
    bool blockingQuery = false;
    std::list<const StateClause*>::iterator recency;
  };
  // Q2 target clauses remain exact assumptions. Cache recurring selectors,
  // but keep their count linear in the SAT state surface so one-off
  // generalization cubes cannot grow the incremental solver without bound.
  std::unordered_map<StateClause, Q2SelectorCacheEntry, StateClauseHash>
      q2SelectorByExclusionClause;
  std::list<const StateClause*> q2BlockingSelectorRecency;
  std::list<const StateClause*> q2StatusSelectorRecency;
  size_t q2SelectorReuseCount = 0;
  size_t q2SelectorEvictionCount = 0;
  // Identity of the shared exact F[0] surface used to build this solver. The
  // vector may only widen; its size therefore detects when extension is needed.
  const std::vector<size_t>* sharedFrameZeroSolverSymbols = nullptr;
  // Higher-frame output batches share this SAT instance. Property and frame
  // clauses are guarded by the active context literal, so only CaDiCaL's
  // context-independent learned clauses can affect a later batch.
  GuardedPredecessorFrameContext guardedContext;
  size_t guardedContextCount = 0;

  int q2SelectorFor(const StateClause& exclusionClause,
                    size_t frame,
                    bool blockingQuery);
  size_t retireStatusQ2Selectors();

  bool canExtendTo(const PredecessorAssumptionCacheKey& candidate) const {
    return key.hasSameReusableContext(candidate) &&
           std::includes(
               candidate.solverSymbols.begin(),
               candidate.solverSymbols.end(),
               key.solverSymbols.begin(),
               key.solverSymbols.end());
  }

  bool hasSameSharedFrameZeroContext(
      const KInductionProblem* problem,
      const void* transitionModel,
      BoolExpr* initFormula,
      const std::vector<size_t>& solverSymbols) const {
    return sharedFrameZeroSolverSymbols == &solverSymbols &&
           key.problem == problem &&
           key.transitionModel == transitionModel &&
           key.initFormula == initFormula &&
           key.frameInvariant == nullptr &&
           key.level == 0 &&
           key.solverSymbols.size() == solverSymbols.size();
  }

  void extendSymbolSurface(const ComplementPartnerIndex& stateRelations,
                           const std::vector<size_t>& solverSymbols);

  const PreparedPredecessorTargetAssumptions& prepareTargetAssumptions(
      const TransitionExprResolver& transitionByState,
      size_t frame,
      const StateClause& targetIdentity,
      const std::vector<TransitionEncodingLiteralGroup>& groups);

};

struct InitIntersectionAssumptionSolver {
  const KInductionProblem* problem = nullptr;
  const BoolExpr* initFormula = nullptr;
  KEPLER_FORMAL::Config::SolverType requestedSolverType =
      KEPLER_FORMAL::Config::SolverType::KISSAT;
  std::unique_ptr<SATSolverWrapper> solver;
  std::unique_ptr<FrameVariableStore> variables;
  // F[0] is immutable. Cache exact cube-intersection answers alongside its
  // incremental SAT solver so repeated obligation/generalization checks do not
  // solve the same assumptions again.
  std::unordered_map<StateCube, bool, StateCubeHash> resultByCube;
};

struct SharedPredecessorAssumptionSolverPool {
  struct Selection {
    std::unique_ptr<PredecessorAssumptionSolver>* solver = nullptr;
    size_t entryIndex = 0;
    bool selectedNewRun = false;
    bool cacheHit = false;
    bool evicted = false;
    size_t closestEntry = 0;
    size_t closestOverlap = 0;
    bool restarted = false;
    size_t retiredContexts = 0;
  };

  Selection select(size_t runId,
                   const std::vector<size_t>& familySolverSymbols,
                   bool usePathLocalReuse) {
    if (selectedRunId == runId && selectedEntry.has_value()) {
      return {&entries[*selectedEntry].solver,
              *selectedEntry,
              false,
              true,
              false,
              *selectedEntry,
              familySolverSymbols.size(),
              false,
              0};
    }

    size_t entryIndex = entries.size();
    bool cacheHit = false;
    size_t closestEntry = 0;
    size_t closestOverlap = 0;
    for (size_t index = 0; index < entries.size(); ++index) {
      const std::vector<size_t>& reusableFamily =
          usePathLocalReuse ? entries[index].lastPathFamilySolverSymbols
                            : entries[index].rootFamilySolverSymbols;
      const size_t overlap = sortedIntersectionSize(
          reusableFamily, familySolverSymbols);
      if (overlap > closestOverlap) {
        closestEntry = index;
        closestOverlap = overlap;
      }
      // Strict equality follows one nested depth-first path. Guarded equality
      // reuses only an exact family: distinct siblings otherwise accumulate
      // inactive property contexts that make later SAT calls much slower.
      const bool canReuseFamily =
          usePathLocalReuse
              ? overlap == familySolverSymbols.size()
              : reusableFamily == familySolverSymbols;
      if (canReuseFamily &&
          (entryIndex == entries.size() ||
           reusableFamily.size() <
               (usePathLocalReuse
                    ? entries[entryIndex].lastPathFamilySolverSymbols.size()
                    : entries[entryIndex].rootFamilySolverSymbols.size()) ||
           (reusableFamily.size() ==
                (usePathLocalReuse
                     ? entries[entryIndex].lastPathFamilySolverSymbols.size()
                     : entries[entryIndex].rootFamilySolverSymbols.size()) &&
            entries[index].lastUse > entries[entryIndex].lastUse))) {
        entryIndex = index;
        cacheHit = true;
      }
    }

    bool evicted = false;
    if (entryIndex == entries.size()) {
      if (entries.size() < kMaxSharedPredecessorSolversPerLevel) {
        entries.push_back(
            {familySolverSymbols, familySolverSymbols, nullptr, 0});
      } else {
        entryIndex = 0;
        for (size_t index = 1; index < entries.size(); ++index) {
          if (entries[index].lastUse < entries[entryIndex].lastUse) {
            entryIndex = index;
          }
        }
        evicted = entries[entryIndex].solver != nullptr;
        entries[entryIndex].rootFamilySolverSymbols = familySolverSymbols;
        entries[entryIndex].lastPathFamilySolverSymbols = familySolverSymbols;
        entries[entryIndex].solver.reset();
      }
    }

    if (usePathLocalReuse) {
      entries[entryIndex].lastPathFamilySolverSymbols = familySolverSymbols;
    }
    bool restarted = false;
    size_t retiredContexts = 0;
    if (!usePathLocalReuse && entries[entryIndex].solver != nullptr &&
        entries[entryIndex].solver->guardedContextCount >=
            kMaxRetiredGuardedPredecessorContexts) {
      // The model and property-family key remain reusable. Only retire the SAT
      // cache whose disabled guarded contexts have reached the bounded window.
      retiredContexts = entries[entryIndex].solver->guardedContextCount;
      entries[entryIndex].solver.reset();
      restarted = true;
    }
    entries[entryIndex].lastUse = ++useSequence;
    selectedRunId = runId;
    selectedEntry = entryIndex;
    return {&entries[entryIndex].solver,
            entryIndex,
            true,
            cacheHit,
            evicted,
            closestEntry,
            closestOverlap,
            restarted,
            retiredContexts};
  }

  size_t retireStatusQ2Selectors() {
    size_t retiredCount = 0;
    for (auto& entry : entries) {
      if (entry.solver != nullptr) {
        retiredCount += entry.solver->retireStatusQ2Selectors();
      }
    }
    return retiredCount;
  }

 private:
  struct Entry {
    // Guarded equality uses the exact root family. Strict equality tracks the
    // most recent nested child so sibling cones do not accumulate in one
    // solver while parent-to-child learned clauses remain available.
    std::vector<size_t> rootFamilySolverSymbols;
    std::vector<size_t> lastPathFamilySolverSymbols;
    std::unique_ptr<PredecessorAssumptionSolver> solver;
    size_t lastUse = 0;
  };

  static size_t sortedIntersectionSize(const std::vector<size_t>& lhs,
                                       const std::vector<size_t>& rhs) {
    size_t lhsIndex = 0;
    size_t rhsIndex = 0;
    size_t intersection = 0;
    while (lhsIndex < lhs.size() && rhsIndex < rhs.size()) {
      if (lhs[lhsIndex] < rhs[rhsIndex]) {
        ++lhsIndex;
      } else if (rhs[rhsIndex] < lhs[lhsIndex]) {
        ++rhsIndex;
      } else {
        ++intersection;
        ++lhsIndex;
        ++rhsIndex;
      }
    }
    return intersection;
  }

  std::vector<Entry> entries;
  size_t selectedRunId = 0;
  std::optional<size_t> selectedEntry;
  size_t useSequence = 0;
};

struct PredecessorAssumptionCache {
  // PDR level-local predecessor queries share the same frame context and
  // differ mostly by target cube.
  std::unordered_map<size_t, std::unique_ptr<PredecessorAssumptionSolver>>
      solversByLevel;
  // Figure 6 asks whether many candidate cubes intersect the same exact F[0].
  // Keep that immutable formula encoded and vary only cube assumptions.
  std::unique_ptr<InitIntersectionAssumptionSolver> initIntersectionSolver;
  // F[0] and its transition relation are identical for every output batch.
  // Clauses streamed into it are exact-Init consequences; higher-frame
  // solvers, frame vectors, and proof obligations remain local to one run.
  std::unique_ptr<PredecessorAssumptionSolver>*
      sharedFrameZeroPredecessorSolver = nullptr;
  std::vector<size_t>* sharedFrameZeroPredecessorSymbols = nullptr;
  const KInductionProblem* sharedFrameZeroPredecessorProblem = nullptr;
  const void* sharedFrameZeroTransitionModel = nullptr;
  // Higher PDR frames are property-specific, but their immutable transition
  // solver can cross serial output batches when every local fact is guarded.
  std::unordered_map<size_t, SharedPredecessorAssumptionSolverPool>*
      sharedHigherFrameSolverPools = nullptr;
  const KInductionProblem* sharedHigherFrameProblem = nullptr;
  const void* sharedHigherFrameTransitionModel = nullptr;
  size_t sharedHigherFrameRunId = 0;
  // Select a persistent solver by the complete property family, not by the
  // first predecessor cube. Unrelated outputs often share reset/control bits;
  // using that small first cube would merge their otherwise distinct cones.
  const std::vector<size_t>* sharedHigherFrameFamilySymbols = nullptr;
  bool usePathLocalHigherFrameSolverReuse = false;
  // Full predecessor-query result cache. SAT entries are keyed by the exact
  // frame fingerprint; UNSAT entries also get a fingerprint-free key because
  // PDR frames only strengthen over time, so a proven-empty predecessor set
  // remains empty after more clauses are learned.
  PredecessorQueryResultStore queryResultStore;
  // Level-zero predecessor queries depend only on exact F[0], T, and the
  // target cube. Keep their exact answers with the validated shared model;
  // property-specific higher-frame answers remain local to this PDR run.
  PredecessorQueryResultStore* sharedFrameZeroQueryResultStore = nullptr;
  const KInductionProblem* sharedFrameZeroQueryProblem = nullptr;
  const TransitionExprResolver* sharedFrameZeroQueryTransition = nullptr;
  // A predecessor UNSAT core for cube U also proves UNSAT for every later
  // target cube that contains U under the same PDR context. Keep those cores
  // separately from exact target results so neighboring dual-rail cubes can
  // reuse the proof without re-solving a wider assumption set.
  // Local dual-rail leaves repeatedly ask nearly identical predecessor
  // questions while obligations move between frames. Keep each frame's solver
  // surface monotonic without carrying symbols into a distinct frame solver.
  detail::PdrFrameSymbolSurfaceCache predecessorSolverSymbolSurfaces;
  // IC3 obligations move between levels. Retain each exact frame fingerprint
  // independently so returning to F[i] does not rebuild its unchanged symbol
  // surface after a query at F[j].
  std::unordered_map<size_t, PredecessorFrameSymbolSurface>
      currentFrameSymbolsByLevel;
  PredecessorTargetSurfaceCache targetSurfaces{
      kMaxPredecessorTargetSurfaceCacheBytes};
  PredecessorTargetSurfaceCache* sharedTargetSurfaces = nullptr;
  // Relation clauses are selected through an immutable model index. The index
  // changes query preparation cost only; it emits the same clauses in the same
  // order as the original full-vector scans.
  const ComplementPartnerIndex* stateRelations = nullptr;
  // Exact Init facts are immutable and already built once per PDR run. Reuse
  // their assignment index for the cheap pre-SAT contradiction check.
  const InitFactIndex* initFacts = nullptr;
};

void retireGeneralizationStatusQ2Selectors(
    PredecessorAssumptionCache* cache) {
  if (cache == nullptr) {
    return;
  }
  for (auto& [level, solver] : cache->solversByLevel) {
    (void)level;
    if (solver != nullptr) {
      solver->retireStatusQ2Selectors();
    }
  }
  if (cache->sharedFrameZeroPredecessorSolver != nullptr &&
      *cache->sharedFrameZeroPredecessorSolver != nullptr) {
    (*cache->sharedFrameZeroPredecessorSolver)->retireStatusQ2Selectors();
  }
  if (cache->sharedHigherFrameSolverPools != nullptr) {
    for (auto& [level, pool] : *cache->sharedHigherFrameSolverPools) {
      (void)level;
      pool.retireStatusQ2Selectors();
    }
  }
}

struct BadCubeAssumptionCacheKey {
  const KInductionProblem* problem = nullptr;
  const BoolExpr* initFormula = nullptr;
  const BoolExpr* frameInvariant = nullptr;
  size_t level = 0;
  std::vector<size_t> solverSymbols;

  bool operator==(const BadCubeAssumptionCacheKey& other) const {
    return problem == other.problem &&
           initFormula == other.initFormula &&
           frameInvariant == other.frameInvariant &&
           level == other.level &&
           solverSymbols == other.solverSymbols;
  }
};

struct BadCubeAssumptionSolver {
  BadCubeAssumptionCacheKey key;
  std::unique_ptr<SATSolverWrapper> solver;
  std::unique_ptr<FrameVariableStore> variables;
  std::unique_ptr<FrameFormulaEncoder> encoder;
  std::unordered_map<BoolExpr*, int> encodedBadRoots;
  std::unordered_set<size_t> querySymbolSet;
  std::unordered_set<StateClause, StateClauseHash> emittedFrameClauses;
  size_t emittedFrameFingerprint = 0;
  size_t emittedFrameLogOffset = 0;
};

struct BadCubeAssumptionCache {
  // Bad-cube searches repeatedly ask the same frame context with different
  // output-bad roots. Keep frame facts permanent and vary only the root
  // literal as a solver assumption.
  std::unique_ptr<BadCubeAssumptionSolver> solver;
};

}  // namespace

struct PDRExactInitCache::Impl {
  Impl(const KInductionProblem& source,
       KEPLER_FORMAL::Config::SolverType requestedSolverType)
      : sourceProblem(&source), solverType(requestedSolverType),
        transitionByState(source), stateRelations(source) {
    validatedProblems.insert(&source);
  }

  bool hasSameDualRailPairs(const KInductionProblem& candidate) const {
    if (sourceProblem->dualRailStatePairs.size() !=
        candidate.dualRailStatePairs.size()) {
      return false;
    }
    for (size_t index = 0; index < sourceProblem->dualRailStatePairs.size();
         ++index) {
      const auto& sourcePair = sourceProblem->dualRailStatePairs[index];
      const auto& candidatePair = candidate.dualRailStatePairs[index];
      if (sourcePair.mayBeOne != candidatePair.mayBeOne ||
          sourcePair.mayBeZero != candidatePair.mayBeZero) {
        return false;
      }
    }
    return true;
  }

  bool matches(const KInductionProblem& candidate,
               KEPLER_FORMAL::Config::SolverType candidateSolverType) const {
    // The SEC strategy mutates only output/property fields on two reusable
    // batch objects.  Once one of those objects has passed the complete model
    // check, its immutable transition identity does not need another ASIC-size
    // compare.
    if (solverType == candidateSolverType &&
        validatedProblems.contains(&candidate)) {
      return true;
    }
    if (sourceProblem == nullptr || solverType != candidateSolverType ||
        sourceProblem->resetBootstrapCycles != candidate.resetBootstrapCycles) {
      return false;
    }

    // Only output/property fields may differ between users of this cache.
    // Comparing the complete reset/transition model prevents stale F[0] reuse
    // if a caller accidentally passes a problem from another SEC model.
    const bool sameModel =
        sourceProblem->inputSymbols == candidate.inputSymbols &&
        sourceProblem->resetBootstrapInputs == candidate.resetBootstrapInputs &&
        sourceProblem->initialStateAssignments ==
            candidate.initialStateAssignments &&
        sourceProblem->bootstrapStateAssignments ==
            candidate.bootstrapStateAssignments &&
        sourceProblem->state0Symbols == candidate.state0Symbols &&
        sourceProblem->state1Symbols == candidate.state1Symbols &&
        sourceProblem->auxiliaryStateSymbols ==
            candidate.auxiliaryStateSymbols &&
        sourceProblem->allSymbols == candidate.allSymbols &&
        sourceProblem->complementedStatePairs0 ==
            candidate.complementedStatePairs0 &&
        sourceProblem->complementedStatePairs1 ==
            candidate.complementedStatePairs1 &&
        sourceProblem->sameFrameStateEqualityPairs0 ==
            candidate.sameFrameStateEqualityPairs0 &&
        sourceProblem->sameFrameStateEqualityPairs1 ==
            candidate.sameFrameStateEqualityPairs1 &&
        hasSameDualRailPairs(candidate) &&
        sourceProblem->transitions0 == candidate.transitions0 &&
        sourceProblem->transitions1 == candidate.transitions1 &&
        sourceProblem->auxiliaryTransitions == candidate.auxiliaryTransitions &&
        sourceProblem->lazyTransitions == candidate.lazyTransitions &&
        sourceProblem->initialCondition == candidate.initialCondition;
    if (sameModel) {
      validatedProblems.insert(&candidate);
    }
    return sameModel;
  }

  const KInductionProblem* sourceProblem = nullptr;
  KEPLER_FORMAL::Config::SolverType solverType;
  BoolExpr* initFormula = nullptr;
  std::unique_ptr<PredecessorAssumptionSolver> frameZeroPredecessorSolver;
  std::vector<size_t> frameZeroPredecessorSymbols;
  PredecessorQueryResultStore frameZeroPredecessorResults;
  std::unordered_map<size_t, SharedPredecessorAssumptionSolverPool>
      higherFramePredecessorSolverPools;
  size_t nextHigherFrameRunId = 1;
  // A frame clause from an unfinished property run is level-specific, not an
  // absolute invariant. Keep it only as an invariant-finder candidate and
  // expose solely the subset certified against this exact I and T.
  std::vector<StateClause> reusableInvariantCandidates;
  FrameClauses reusableInvariant;
  size_t reusableInvariantCandidateLiteralCount = 0;
  size_t reusableInvariantCandidateRevision = 0;
  size_t reusableInvariantCertifiedRevision = 0;
  std::optional<SATSolverWrapper::CadicalWorkBudget>
      reusableInvariantCertificationBudget;
  bool reusableInvariantCertificationDisabled = false;
  // These structures depend only on the validated transition model. SEC runs
  // output batches serially, so every batch can reuse their exact contents
  // without sharing property-specific proof state or changing query order.
  TransitionExprResolver transitionByState;
  ComplementPartnerIndex stateRelations;
  std::shared_ptr<PdrFormulaSupportCache> formulaSupportCache;
  std::optional<InitFactIndex> initFacts;
  // Target-surface entries contain only exact transition-support preparation.
  // They may cross output batches, unlike SAT answers and learned proof state.
  PredecessorTargetSurfaceCache targetSurfaces{
      kMaxPredecessorTargetSurfaceCacheBytes};
  size_t immutableMetadataUses = 0;
  mutable std::unordered_set<const KInductionProblem*> validatedProblems;
};

PDRExactInitCache::PDRExactInitCache(
    const KInductionProblem& sourceProblem,
    KEPLER_FORMAL::Config::SolverType solverType)
    : impl_(std::make_unique<Impl>(sourceProblem, solverType)) {}

PDRExactInitCache::~PDRExactInitCache() = default;

namespace {

enum class PdrBudgetExhaustion {
  None,
  LocalQuery,
};

thread_local PdrBudgetExhaustion pdrBudgetExhaustion =
    PdrBudgetExhaustion::None;
thread_local size_t pdrPredecessorQueryLimit = 0;
thread_local const PDRQueryLimits* activePdrQueryLimits = nullptr;

class ScopedPdrQueryLimits {
 public:
  explicit ScopedPdrQueryLimits(const PDRQueryLimits* limits)
      : previous_(activePdrQueryLimits) {
    activePdrQueryLimits = limits;
  }

  ~ScopedPdrQueryLimits() { activePdrQueryLimits = previous_; }

 private:
  const PDRQueryLimits* previous_ = nullptr;
};

bool pdrStatsEnabled();

void resetPdrBudgetExhaustion() {
  pdrBudgetExhaustion = PdrBudgetExhaustion::None;
}

void setPdrPredecessorQueryLimit(size_t limit) {
  pdrPredecessorQueryLimit = limit;
}

void markPdrBudgetExhausted(PdrBudgetExhaustion reason) {
  if (pdrBudgetExhaustion == PdrBudgetExhaustion::None) {
    pdrBudgetExhaustion = reason;
  }
}

bool hasPdrBudgetExhaustion() {
  return pdrBudgetExhaustion != PdrBudgetExhaustion::None;
}

bool consumePdrPredecessorQueryBudget(size_t* remainingQueries) {
  if (remainingQueries == nullptr) {
    return true;
  }
  if (*remainingQueries == 0) {
    if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC PDR stats: predecessor query-count budget exhausted limit=",
          pdrPredecessorQueryLimit);
    }  // LCOV_EXCL_LINE
    markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);  // LCOV_EXCL_LINE
    return false;  // LCOV_EXCL_LINE
  }
  --(*remainingQueries);
  return true;
}

bool pdrStatsEnabled() {
  return std::getenv("KEPLER_SEC_PDR_STATS") != nullptr;
}

size_t pdrTransitionSourceCount(const KInductionProblem& problem) {
  size_t count = problem.transitions0.size() + problem.transitions1.size() +
                 problem.auxiliaryTransitions.size();
  if (problem.lazyTransitions != nullptr) {
    count += problem.lazyTransitions->sourceByStateSymbol.size();
  }
  return count;
}

size_t envSizeLimitOrDefault(const char* name, size_t defaultValue);

size_t pdrStatsInterval() {
  const char* intervalText = std::getenv("KEPLER_SEC_PDR_STATS_INTERVAL");
  if (intervalText == nullptr || *intervalText == '\0') {
    return kDefaultPdrStatsInterval;  // LCOV_EXCL_LINE
  }

  const auto interval = std::strtoull(intervalText, nullptr, 10);
  return interval == 0 ? 1 : static_cast<size_t>(interval);
}

size_t envSizeLimitOrDefault(const char* name, size_t defaultValue) {
  const char* valueText = std::getenv(name);
  if (valueText == nullptr || *valueText == '\0') {
    return defaultValue;
  }
  const auto value = std::strtoull(valueText, nullptr, 10);
  return value == 0 ? defaultValue : static_cast<size_t>(value);
}

unsigned envUnsignedLimitOrDefaultAllowZero(const char* name,
                                            // LCOV_EXCL_START
                                            unsigned defaultValue) {
                                            // LCOV_EXCL_STOP
  const char* valueText = std::getenv(name);
  if (valueText == nullptr || *valueText == '\0') {
    return defaultValue;
  }
  const auto value = std::strtoull(valueText, nullptr, 10);
  return value > std::numeric_limits<unsigned>::max()
             ? std::numeric_limits<unsigned>::max()  // LCOV_EXCL_LINE
             : static_cast<unsigned>(value);
}

unsigned dualRailBadCubeConflictLimit() {
  return envUnsignedLimitOrDefaultAllowZero(
      "KEPLER_SEC_PDR_DUAL_RAIL_BAD_CUBE_CONFLICT_LIMIT",
      kDefaultDualRailBadCubeConflictLimit);
}

unsigned dualRailPredecessorConflictLimit(PredecessorQueryPurpose purpose) {
  if (activePdrQueryLimits != nullptr) {
    return purpose == PredecessorQueryPurpose::BlockObligation
               ? activePdrQueryLimits->blockingConflictLimit
               : activePdrQueryLimits->predecessorConflictLimit;
  }
  const unsigned configuredLimit = envUnsignedLimitOrDefaultAllowZero(
      kDualRailPredecessorConflictLimitEnv,
      kDefaultDualRailPredecessorConflictLimit);
  if (std::getenv(kDualRailPredecessorConflictLimitEnv) != nullptr ||
      purpose != PredecessorQueryPurpose::BlockObligation) {
    return configuredLimit;
  }
  return std::max(configuredLimit, kDefaultDualRailBlockingConflictLimit);
}

unsigned dualRailPredecessorDecisionLimit(PredecessorQueryPurpose purpose) {
  if (activePdrQueryLimits != nullptr) {
    return purpose == PredecessorQueryPurpose::BlockObligation
               ? activePdrQueryLimits->blockingDecisionLimit
               : activePdrQueryLimits->predecessorDecisionLimit;
  }
  const unsigned configuredLimit = envUnsignedLimitOrDefaultAllowZero(
      kDualRailPredecessorDecisionLimitEnv,
      kDefaultDualRailPredecessorDecisionLimit);
  if (std::getenv(kDualRailPredecessorDecisionLimitEnv) != nullptr ||
      purpose != PredecessorQueryPurpose::BlockObligation) {
    return configuredLimit;
  }
  return std::max(configuredLimit, kDefaultDualRailBlockingDecisionLimit);
}

unsigned boundedNarrowGeneralizationProbeLimit(unsigned fullLimit,
                                               unsigned probeLimit) {
  return fullLimit == 0 ? probeLimit : std::min(fullLimit, probeLimit);
}

size_t dualRailPredecessorEncodingNodeLimit() {
  if (activePdrQueryLimits != nullptr &&
      activePdrQueryLimits->predecessorEncodingNodeLimit != 0) {
    return activePdrQueryLimits->predecessorEncodingNodeLimit;
  }
  return envSizeLimitOrDefault(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_ENCODING_NODE_LIMIT",
      kDefaultDualRailPredecessorEncodingNodeLimit);
}

size_t dualRailPredecessorNodeHintTargetLimit() {
  if (activePdrQueryLimits != nullptr &&
      activePdrQueryLimits->predecessorNodeHintTargetLimit != 0) {
    return activePdrQueryLimits->predecessorNodeHintTargetLimit;
  }
  return std::numeric_limits<size_t>::max();
}

size_t dualRailPredecessorEncodingSupportLimit() {
  return envSizeLimitOrDefault(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_ENCODING_SUPPORT_LIMIT",
      kDefaultDualRailPredecessorEncodingSupportLimit);
}

size_t nextPdrPredecessorQueryNumber() {
  // The stats path is intentionally process-local and diagnostic-only. PDR is
  // currently run serially per SEC output slice, so a simple counter gives a
  // stable view of where a long proof is spending time without touching the
  // proof algorithm or adding synchronization overhead to normal runs.
  static size_t queryNumber = 0;
  return ++queryNumber;
}


size_t nextPdrBadCubeQueryNumber() {
  static size_t queryNumber = 0;
  return ++queryNumber;
}

size_t nextPdrDualRailPredecessorCoreSkipNumber() {
  static size_t skipNumber = 0;
  return ++skipNumber;
}

bool shouldEmitPdrStats(size_t queryNumber) {
  if (!pdrStatsEnabled()) {
    return false;
  }
  return queryNumber <= kInitialPdrStatsQueries ||
         queryNumber % pdrStatsInterval() == 0;
}

bool shouldEmitFrequentPdrStats() {
  if (!pdrStatsEnabled()) {
    return false;
  }
  // Full diagnostics retain representative cache/allocation events without a
  // synchronous write for every SAT query. Interval 1 remains fully verbose
  // for focused unit tests.
  static size_t eventNumber = 0;
  return shouldEmitPdrStats(++eventNumber);
}

class PdrFormulaSupportCache {
  // LCOV_EXCL_START
 public:
  // LCOV_EXCL_STOP
  struct TernaryNode {
    Op op = Op::NONE;
    size_t symbol = 0;
    size_t left = static_cast<size_t>(-1);
    size_t right = static_cast<size_t>(-1);
  };

  struct TernaryEvaluationMemoEntry {
    size_t generation = 0;
    size_t queuedGeneration = 0;
    std::optional<bool> value;
  };

  using TernaryEvaluationMemo = std::vector<TernaryEvaluationMemoEntry>;

  static constexpr size_t kInvalidTernaryNode = static_cast<size_t>(-1);

  PdrFormulaSupportCache() = default;

  const std::set<size_t>& support(BoolExpr* formula) {
    static const std::set<size_t> emptySupport;
    if (formula == nullptr) {
      return emptySupport;  // LCOV_EXCL_LINE
    }
    if (const auto it = supportByExpr_.find(formula);
        it != supportByExpr_.end()) {
      return it->second;
    }
    const auto [it, _] =
        supportByExpr_.emplace(formula, formula->getSupportVars());
    return it->second;
  }

  const std::vector<size_t>& mappedTernarySupport(
      BoolExpr* formula,
      const std::unordered_map<size_t, size_t>* symbolMap) {
    static const std::vector<size_t> emptySupport;
    if (formula == nullptr) {
      return emptySupport;
    }

    const TernaryMappedSupportKey key{formula, symbolMap};
    if (const auto existing = mappedTernarySupportByRoot_.find(key);
        existing != mappedTernarySupportByRoot_.end()) {
      ++mappedTernarySupportReuseCount_;
      return existing->second;
    }

    std::vector<size_t> mappedSupport;
    mappedSupport.reserve(support(formula).size());
    for (const size_t localSymbol : support(formula)) {
      size_t mappedSymbol = localSymbol;
      if (localSymbol >= 2 && symbolMap != nullptr) {
        const auto mapped = symbolMap->find(localSymbol);
        if (mapped == symbolMap->end()) {
          continue;
        }
        mappedSymbol = mapped->second;
      }
      if (mappedSymbol >= 2) {
        mappedSupport.push_back(mappedSymbol);
      }
    }
    std::sort(mappedSupport.begin(), mappedSupport.end());
    mappedSupport.erase(
        std::unique(mappedSupport.begin(), mappedSupport.end()),
        mappedSupport.end());
    const auto [inserted, insertedNew] = mappedTernarySupportByRoot_.emplace(
        key, std::move(mappedSupport));
    (void)insertedNew;
    return inserted->second;
  }

  const std::vector<size_t>& relationClosedSupport(
      BoolExpr* formula,
      const ComplementPartnerIndex& complementPartners);

  size_t ternaryNodeIndex(BoolExpr* formula) {
    if (formula == nullptr) {
      return kInvalidTernaryNode;
    }
    if (const auto existing = ternaryNodeIndexByExpr_.find(formula);
        existing != ternaryNodeIndexByExpr_.end()) {
      return existing->second;
    }

    // BoolExpr DAGs are immutable. Compile each node once so the paper's
    // repeated ternary trials can use dense vector slots instead of hashing
    // the same expression pointer at every recursive visit.
    const size_t index = ternaryNodes_.size();
    ternaryNodeIndexByExpr_.emplace(formula, index);
    ternaryNodes_.push_back(
        {formula->getOp(), formula->getId(), kInvalidTernaryNode,
         kInvalidTernaryNode});
    ternaryParentNodes_.emplace_back();
    if (formula->getOp() == Op::VAR) {
      ternaryVariableNodesBySymbol_[formula->getId()].push_back(index);
    }
    const size_t left = ternaryNodeIndex(formula->getLeft());
    const size_t right = ternaryNodeIndex(formula->getRight());
    ternaryNodes_[index].left = left;
    ternaryNodes_[index].right = right;
    if (left != kInvalidTernaryNode) {
      ternaryParentNodes_[left].push_back(index);
    }
    if (right != kInvalidTernaryNode && right != left) {
      ternaryParentNodes_[right].push_back(index);
    }
    return index;
  }

  const TernaryNode& ternaryNode(size_t index) const {
    return ternaryNodes_[index];
  }

  const std::vector<size_t>& ternaryParentNodes(size_t index) const {
    return ternaryParentNodes_[index];
  }

  const std::vector<size_t>& ternaryVariableNodes(size_t symbol) const {
    static const std::vector<size_t> noNodes;
    const auto nodes = ternaryVariableNodesBySymbol_.find(symbol);
    return nodes == ternaryVariableNodesBySymbol_.end()
               ? noNodes
               : nodes->second;
  }

  size_t ternaryNodeCount() const { return ternaryNodes_.size(); }

  TernaryEvaluationMemo& ternaryEvaluationMemo(
      const std::unordered_map<size_t, size_t>* symbolMap) {
    auto& memo = ternaryEvaluationMemosBySymbolMap_[symbolMap];
    if (memo.size() < ternaryNodes_.size()) {
      memo.resize(ternaryNodes_.size());
    }
    return memo;
  }

  size_t nextTernaryEvaluationGeneration() {
    return ++ternaryEvaluationGeneration_;
  }

  void emitMappedTernarySupportStats() const {
    emitSecDiag(
        "SEC PDR stats: ternary mapped support cache reused=",
        mappedTernarySupportReuseCount_,
        " entries=",
        mappedTernarySupportByRoot_.size());
  }

 private:
  struct TernaryMappedSupportKey {
    BoolExpr* formula = nullptr;
    const std::unordered_map<size_t, size_t>* symbolMap = nullptr;

    bool operator==(const TernaryMappedSupportKey& other) const {
      return formula == other.formula && symbolMap == other.symbolMap;
    }
  };

  struct TernaryMappedSupportKeyHash {
    size_t operator()(const TernaryMappedSupportKey& key) const {
      size_t seed = std::hash<const void*>()(key.formula);
      mixHashValue(seed, std::hash<const void*>()(key.symbolMap));
      return seed;
    }
  };

  // PDR rebuilds many local SAT queries over the same frame/property formulas.
  // Memoizing formula support avoids repeatedly walking large BoolExpr DAGs
  // while keeping each query's selected symbol set unchanged.
  std::unordered_map<BoolExpr*, std::set<size_t>> supportByExpr_;
  // Bad-state queries combine immutable formula support with changing frame
  // clauses. Cache the formula side after applying the same state-relation
  // closure so ASIC-size invariants are not copied into a hash set per output.
  const ComplementPartnerIndex* closedSupportRelations_ = nullptr;
  std::unordered_map<BoolExpr*, std::vector<size_t>> closedSupportByExpr_;
  std::unordered_map<BoolExpr*, size_t> ternaryNodeIndexByExpr_;
  std::vector<TernaryNode> ternaryNodes_;
  // Parent links let a changed model literal propagate its exact ternary value
  // only through computed ancestors. This is a cache dependency graph over
  // immutable BoolExpr nodes; it does not alter the paper's simulation or
  // literal-trial order.
  std::vector<std::vector<size_t>> ternaryParentNodes_;
  std::unordered_map<size_t, std::vector<size_t>>
      ternaryVariableNodesBySymbol_;
  // Transition roots and same-design symbol maps are immutable for the life
  // of a PDR run. Retain their exact mapped support instead of rebuilding an
  // unordered set for every SAT predecessor model.
  std::unordered_map<TernaryMappedSupportKey,
                     std::vector<size_t>,
                     TernaryMappedSupportKeyHash>
      mappedTernarySupportByRoot_;
  // PDR output batches run serially. Retain the dense memo allocations between
  // exact predecessor queries and distinguish every trial by generation.
  std::unordered_map<const std::unordered_map<size_t, size_t>*,
                     TernaryEvaluationMemo>
      ternaryEvaluationMemosBySymbolMap_;
  size_t ternaryEvaluationGeneration_ = 0;
  size_t mappedTernarySupportReuseCount_ = 0;
};

void addCubeAssumptions(SATSolverWrapper& solver,
                        const FrameVariableStore& variables,
                        const StateCube& cube,
                        size_t frame);

void addNegatedCubeClause(SATSolverWrapper& solver,
                          const FrameVariableStore& variables,
                          const StateCube& cube,
                          size_t frame);

StateClause clauseFromCube(const StateCube& cube);

std::vector<size_t> cubeStateSymbols(const StateCube& cube);

void addPostBootstrapResetInputConstraints(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t frame);

void addFormulaSymbols(BoolExpr* formula,
                       std::unordered_set<size_t>& symbols,
                       PdrFormulaSupportCache* supportCache = nullptr);

bool predecessorSourceFrameIsKnownSafe(size_t level);

void normalizeCube(StateCube& cube);

void addRelevantComplementedStatePartners(
    const ComplementPartnerIndex& complementPartners,
    std::unordered_set<size_t>& symbols);

// LCOV_EXCL_START
std::vector<size_t> sortUniqueSymbols(std::unordered_set<size_t> symbols) {
// LCOV_EXCL_STOP
  std::vector<size_t> ordered(symbols.begin(), symbols.end());
  std::sort(ordered.begin(), ordered.end());
  ordered.erase(std::unique(ordered.begin(), ordered.end()), ordered.end());
  return ordered;
}

std::optional<std::vector<size_t>> collectBoundedStateSupportSymbols(
    BoolExpr* formula,
    size_t maxVisitedNodes,
    size_t maxStateSymbols,
    const std::unordered_set<size_t>& stateSymbolSet) {
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return {};  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  std::unordered_set<size_t> stateSupport;
  std::unordered_set<const BoolExpr*> visited;
  std::vector<const BoolExpr*> stack{formula};
  while (!stack.empty()) {
    const BoolExpr* node = stack.back();
    stack.pop_back();
    if (node == nullptr || !visited.insert(node).second) {
      continue;  // LCOV_EXCL_LINE
    }
    if (visited.size() > maxVisitedNodes) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    if (node->getOp() == Op::VAR) {
      if (stateSymbolSet.find(node->getId()) != stateSymbolSet.end()) {
        stateSupport.insert(node->getId());
        if (maxStateSymbols != 0 && stateSupport.size() > maxStateSymbols) {
          return std::nullopt;
        }
      }
      continue;
    }
    if (node->getRight() != nullptr) {
      stack.push_back(node->getRight());
    }
    if (node->getLeft() != nullptr) {
      stack.push_back(node->getLeft());
    }
  }
  return sortUniqueSymbols(std::move(stateSupport));
}

std::vector<size_t> retainPdrStateSymbols(
    const std::vector<size_t>& symbols,
    const std::unordered_set<size_t>& stateSymbols) {
  std::vector<size_t> retained;
  retained.reserve(symbols.size());
  for (const size_t symbol : symbols) {
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      retained.push_back(symbol);
    }
  }
  return retained;
}

// LCOV_EXCL_START
std::vector<size_t> expandTransitionTargets(
    const KInductionProblem& problem,
    const std::vector<size_t>& requestedTargets,
    const TransitionExprResolver& transitionByState) {
  const auto& primaryByComplement = transitionByState.primaryByComplement();
  // LCOV_EXCL_STOP
  std::unordered_set<size_t> targets;
  targets.reserve(requestedTargets.size());

  for (const auto symbol : requestedTargets) {
    if (transitionByState.contains(symbol)) {
      targets.insert(symbol);
      continue;
    }
    if (const auto primaryIt = primaryByComplement.find(symbol);  // LCOV_EXCL_LINE
        primaryIt != primaryByComplement.end() &&  // LCOV_EXCL_LINE
        transitionByState.contains(primaryIt->second)) {  // LCOV_EXCL_LINE
      targets.insert(primaryIt->second);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }

  return sortUniqueSymbols(std::move(targets));
}

std::vector<size_t> collectTransitionSupportSymbols(
    const TransitionExprResolver& transitionByState,
    // LCOV_EXCL_START
    const std::vector<size_t>& encodedTargets) {
    // LCOV_EXCL_STOP
  std::unordered_set<size_t> supportSymbols;
  for (const auto stateSymbol : encodedTargets) {
    const auto& support = transitionByState.support(stateSymbol);
    supportSymbols.insert(support.begin(), support.end());
  }
  return sortUniqueSymbols(std::move(supportSymbols));
}

size_t estimateTransitionEncodingNodes(
    const TransitionExprResolver& transitionByState,
    const std::vector<size_t>& encodedTargets,
    size_t targetCountLimit = kMaxExactTransitionNodeCountHintTargets) {
  if (encodedTargets.size() > targetCountLimit) {
    return 0;  // LCOV_EXCL_LINE
  }
  size_t estimate = 0;
  for (const auto stateSymbol : encodedTargets) {
    estimate += transitionByState.nodeCount(stateSymbol);
  }
  return estimate;
}

std::vector<size_t> sortClosedCurrentFrameSymbols(
    const ComplementPartnerIndex& complementPartners,
    std::unordered_set<size_t> symbols);

PredecessorTargetSurface buildPredecessorTargetSurface(
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    const ComplementPartnerIndex& complementPartners,
    const StateCube& targetCube) {
  PredecessorTargetSurface surface;
  surface.targetSymbols = cubeStateSymbols(targetCube);
  surface.encodedTargets =
      expandTransitionTargets(problem, surface.targetSymbols, transitionByState);
  surface.transitionSupportSymbols =
      collectTransitionSupportSymbols(transitionByState, surface.encodedTargets);
  surface.predecessorSymbols = retainPdrStateSymbols(
      surface.transitionSupportSymbols, transitionByState.stateSymbols());
  surface.closedPredecessorSymbols = sortClosedCurrentFrameSymbols(
      complementPartners,
      std::unordered_set<size_t>(surface.predecessorSymbols.begin(),
                                 surface.predecessorSymbols.end()));
  std::unordered_set<size_t> transitionSymbols;
  transitionSymbols.reserve(surface.transitionSupportSymbols.size());
  for (const size_t symbol : surface.transitionSupportSymbols) {
    if (symbol >= 2) {
      transitionSymbols.insert(symbol);
    }
  }
  surface.closedTransitionSupportSymbols = sortClosedCurrentFrameSymbols(
      complementPartners, std::move(transitionSymbols));
  surface.exclusionClause = clauseFromCube(targetCube);
  surface.transitionGroups =
      groupTransitionCubeLiteralsBySymbolMap(transitionByState, targetCube);
  surface.ternaryRoots.reserve(targetCube.size());
  for (const auto& group : surface.transitionGroups) {
    for (const auto& literal : group.literals) {
      const TransitionExprView view =
          transitionByState.expressionView(literal.transitionSymbol);
      surface.ternaryRoots.push_back(
          {view.expr, view.symbolMap, literal.desiredValue});
    }
  }
  surface.transitionEncodingNodes =
      estimateTransitionEncodingNodes(
          transitionByState,
          surface.encodedTargets,
          dualRailPredecessorNodeHintTargetLimit());
  return surface;
}

bool shouldUsePredecessorSolverCache(const KInductionProblem& problem,
                                     size_t level,
                                     size_t querySymbolCount) {
  return level == 0 || !problem.usesDualRailStateEncoding ||
         querySymbolCount <= kMaxDualRailPredecessorSolverCacheStateSymbols;
}

bool shouldUseBadCubeSolverCache(const KInductionProblem& problem,
                                 size_t level,
                                 size_t querySymbolCount) {
  return level == 0 || !problem.usesDualRailStateEncoding ||
         querySymbolCount <= kMaxDualRailBadCubeSolverCacheStateSymbols;
}

size_t predecessorTargetSurfaceBytes(const StateCube& targetCube,
                                     const PredecessorTargetSurface& surface) {
  size_t bytes = targetCube.size() * sizeof(CubeLiteral) +
                 (surface.targetSymbols.size() +
                  surface.encodedTargets.size() +
                  surface.transitionSupportSymbols.size() +
                  surface.predecessorSymbols.size() +
                  surface.closedPredecessorSymbols.size() +
                  surface.closedTransitionSupportSymbols.size()) *
                     sizeof(size_t) +
                 surface.exclusionClause.size() * sizeof(ClauseLiteral);
  for (const auto& group : surface.transitionGroups) {
    bytes += group.literals.size() * sizeof(TransitionEncodingLiteral) +
             group.stateSymbols.size() * sizeof(size_t);
  }
  bytes += surface.ternaryRoots.size() * sizeof(PdrTernarySimulationRoot);
  return bytes;
}

const PredecessorTargetSurface& predecessorTargetSurfaceFor(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    const ComplementPartnerIndex& complementPartners,
    const StateCube& targetCube,
    PredecessorTargetSurface& uncachedSurface) {
  auto& targetSurfaces = cache.sharedTargetSurfaces != nullptr
                             ? *cache.sharedTargetSurfaces
                             : cache.targetSurfaces;
  if (PredecessorTargetSurface* existing = targetSurfaces.find(targetCube);
      existing != nullptr) {
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag("SEC PDR stats: predecessor target surface reused target=",
                  targetCube.size(), " entries=", targetSurfaces.size());
    }
    return *existing;
  }
  uncachedSurface =
      buildPredecessorTargetSurface(
          problem, transitionByState, complementPartners, targetCube);
  const size_t entryBytes =
      predecessorTargetSurfaceBytes(targetCube, uncachedSurface);
  if (entryBytes > kMaxPredecessorTargetSurfaceCacheBytes) {
    return uncachedSurface; // LCOV_EXCL_LINE
  }
  auto inserted = targetSurfaces.insert(
      targetCube, std::move(uncachedSurface), entryBytes);
  if (inserted.evictedEntries != 0 && shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: predecessor target surface cache evicted entries=",
        inserted.evictedEntries,
        " retained_entries=",
        targetSurfaces.size(),
        " bytes=",
        targetSurfaces.retainedWeight());
  }
  if (shouldEmitFrequentPdrStats()) {
    emitSecDiag("SEC PDR stats: predecessor target surface cached target=",
                targetCube.size(),
                " encoded_targets=", inserted.value->encodedTargets.size(),
                " transition_support=",
                inserted.value->transitionSupportSymbols.size(),
                " entries=", targetSurfaces.size(),
                " bytes=", targetSurfaces.retainedWeight());
  }
  return *inserted.value;
}

void appendTransitionEncodingLiteralGroup(
    std::vector<TransitionEncodingLiteralGroup>& groups,
    const std::unordered_map<size_t, size_t>* symbolMap,
    TransitionEncodingLiteral literal) {
  for (auto& group : groups) {
    if (group.symbolMap == symbolMap) {
      group.stateSymbols.push_back(literal.transitionSymbol);
      group.literals.push_back(std::move(literal));
      return;
    }
  }
  TransitionEncodingLiteralGroup group;
  group.symbolMap = symbolMap;
  group.stateSymbols.push_back(literal.transitionSymbol);
  group.literals.push_back(std::move(literal));
  // LCOV_EXCL_START
  groups.push_back(std::move(group));
}

std::vector<TransitionEncodingLiteralGroup> groupTransitionCubeLiteralsBySymbolMap(
// LCOV_EXCL_STOP
    const TransitionExprResolver& transitionByState,
    // LCOV_EXCL_START
    const StateCube& targetCube) {
  const auto& primaryByComplement = transitionByState.primaryByComplement();
  std::vector<TransitionEncodingLiteralGroup> groups;
  // LCOV_EXCL_STOP
  groups.reserve(3);
  for (const auto& literal : targetCube) {
    size_t transitionSymbol = literal.symbol;
    bool desiredValue = literal.value;
    if (!transitionByState.contains(transitionSymbol)) {
      const auto primaryIt = primaryByComplement.find(transitionSymbol);  // LCOV_EXCL_LINE
      if (primaryIt == primaryByComplement.end() ||  // LCOV_EXCL_LINE
          !transitionByState.contains(primaryIt->second)) {  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      transitionSymbol = primaryIt->second;  // LCOV_EXCL_LINE
      desiredValue = !desiredValue;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE

    const TransitionExprView view =
        transitionByState.expressionView(transitionSymbol);
    appendTransitionEncodingLiteralGroup(
        groups,
        view.symbolMap,
        TransitionEncodingLiteral{transitionSymbol, desiredValue, literal});
  }
  for (auto& group : groups) {
    std::sort(group.stateSymbols.begin(), group.stateSymbols.end());
    group.stateSymbols.erase(
        std::unique(group.stateSymbols.begin(), group.stateSymbols.end()),
        group.stateSymbols.end());
  }
  return groups;
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
std::vector<size_t> cubeStateSymbols(const StateCube& cube) {
  std::unordered_set<size_t> symbols;
  symbols.reserve(cube.size());
  for (const auto& literal : cube) {
    symbols.insert(literal.symbol);
  }
  return sortUniqueSymbols(std::move(symbols));
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
bool cubeContainsCube(const StateCube& cube, const StateCube& core) {
  return std::includes(
      cube.begin(),
      cube.end(),
      core.begin(),
      core.end(),
      [](const CubeLiteral& lhs, const CubeLiteral& rhs) {
        if (lhs.symbol != rhs.symbol) {
          return lhs.symbol < rhs.symbol;
        }
        // LCOV_EXCL_START
        return lhs.value < rhs.value;
        // LCOV_EXCL_STOP
      });
}

void addSupportSymbols(const std::set<size_t>& support,
                       std::unordered_set<size_t>& symbols) {
  for (const auto symbol : support) {
    if (symbol >= 2) {
      symbols.insert(symbol);
    }
  }
}


void addFormulaSymbols(BoolExpr* formula,
                       std::unordered_set<size_t>& symbols,
                       PdrFormulaSupportCache* supportCache) {
  if (formula == nullptr) {
    return;  // LCOV_EXCL_LINE
  }
  if (supportCache != nullptr) {
    addSupportSymbols(supportCache->support(formula), symbols);
    return;
  }
  addSupportSymbols(formula->getSupportVars(), symbols);
}


void addRelevantComplementedStatePartners(
    const ComplementPartnerIndex& complementPartners,
    std::unordered_set<size_t>& symbols) {
  complementPartners.addComplementedPartnerClosure(symbols);
}

void addRelevantSameFrameStateEqualityPartners(
    const ComplementPartnerIndex& complementPartners,
    std::unordered_set<size_t>& symbols) {
  complementPartners.addSameFrameEqualityPartnerClosure(symbols);
}

void addRelevantDualRailPartners(
    const ComplementPartnerIndex& complementPartners,
    std::unordered_set<size_t>& symbols) {
  complementPartners.addDualRailPartnerClosure(symbols);
}

const std::vector<size_t>& PdrFormulaSupportCache::relationClosedSupport(
    BoolExpr* formula,
    const ComplementPartnerIndex& complementPartners) {
  static const std::vector<size_t> emptySupport;
  if (formula == nullptr) {
    return emptySupport;  // LCOV_EXCL_LINE
  }
  if (closedSupportRelations_ != &complementPartners) {
    closedSupportByExpr_.clear();  // LCOV_EXCL_LINE
    closedSupportRelations_ = &complementPartners;
  }
  if (const auto existing = closedSupportByExpr_.find(formula);
      existing != closedSupportByExpr_.end()) {
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: formula closed support reused symbols=",
          existing->second.size());
    }
    return existing->second;
  }

  std::unordered_set<size_t> symbols;
  addSupportSymbols(support(formula), symbols);
  addRelevantComplementedStatePartners(complementPartners, symbols);
  addRelevantSameFrameStateEqualityPartners(complementPartners, symbols);
  addRelevantDualRailPartners(complementPartners, symbols);
  auto [inserted, insertedNew] = closedSupportByExpr_.emplace(
      formula, sortUniqueSymbols(std::move(symbols)));
  (void)insertedNew;
  if (shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: formula closed support cached symbols=",
        inserted->second.size());
  }
  return inserted->second;
}

void addCubeSymbols(const StateCube& cube, std::unordered_set<size_t>& symbols) {
  for (const auto& literal : cube) {
    symbols.insert(literal.symbol);
  }
}

void addClauseSymbols(const StateClause& clause, std::unordered_set<size_t>& symbols) {
  for (const auto& literal : clause) {
    symbols.insert(literal.symbol);
  }
}


void addAllFrameClauseSymbols(const FrameClauses& frame,
                              std::unordered_set<size_t>& symbols) {
  for (const auto& clause : frame.clauses) {
    addClauseSymbols(clause, symbols);
  }
}


void addFrameConstraintSymbols(BoolExpr* initFormula,
                               BoolExpr* frameInvariant,
                               const std::vector<FrameClauses>& frames,
                               size_t level,
                               const ComplementPartnerIndex& complementPartners,
                               std::unordered_set<size_t>& symbols,
  PdrFormulaSupportCache* supportCache) {
  if (level == 0) {
    addFormulaSymbols(initFormula, symbols, supportCache);
    addAllFrameClauseSymbols(frames[0], symbols);
  } else {
    addFormulaSymbols(frameInvariant, symbols, supportCache);
    addAllFrameClauseSymbols(frames[level], symbols);
  }
  addRelevantComplementedStatePartners(complementPartners, symbols);
  addRelevantSameFrameStateEqualityPartners(complementPartners, symbols);
  addRelevantDualRailPartners(complementPartners, symbols);
}

std::vector<size_t> findBadQuerySymbols(BoolExpr* initFormula,
                                        BoolExpr* frameInvariant,
                                        const std::vector<FrameClauses>& frames,
                                        BoolExpr* badFormula,
                                        size_t level,
                                        const ComplementPartnerIndex& complementPartners,
                                        PdrFormulaSupportCache* supportCache) {
  if (supportCache != nullptr) {
    // Relation closure distributes over set union. Close the immutable formula
    // once, close only the current bad root and learned frame clauses here, and
    // merge the sorted sets to produce the exact original query surface.
    const std::vector<size_t>& stableFormulaSymbols =
        supportCache->relationClosedSupport(
            level == 0 ? initFormula : frameInvariant,
            complementPartners);
    std::unordered_set<size_t> dynamicSymbols;
    addFormulaSymbols(badFormula, dynamicSymbols, supportCache);
    addAllFrameClauseSymbols(frames[level], dynamicSymbols);
    addRelevantComplementedStatePartners(
        complementPartners, dynamicSymbols);
    addRelevantSameFrameStateEqualityPartners(
        complementPartners, dynamicSymbols);
    addRelevantDualRailPartners(complementPartners, dynamicSymbols);
    return detail::mergeSortedPdrSymbolVectors(
        stableFormulaSymbols,
        sortUniqueSymbols(std::move(dynamicSymbols)));
  }

  std::unordered_set<size_t> symbols;
  addFormulaSymbols(badFormula, symbols, supportCache);
  addFrameConstraintSymbols(
      initFormula,
      frameInvariant,
      frames,
      level,
      complementPartners,
      symbols,
      supportCache);
  return sortUniqueSymbols(std::move(symbols));
}

void prepareSharedExactInitQueries(
    PDRExactInitCache::Impl& cache,
    BoolExpr* initFormula,
    const ComplementPartnerIndex& complementPartners,
    PdrFormulaSupportCache* supportCache) {
  if (!cache.frameZeroPredecessorSymbols.empty()) {
    return;
  }

  // The full source output surface covers every guarded, strict, and split bad
  // root that this top-level SEC call may ask about. F[0] itself supplies the
  // reset-history symbols. Build this stable union once so the incremental SAT
  // solver never has to be rebuilt merely because the next batch is wider.
  std::unordered_set<size_t> symbols;
  addFormulaSymbols(initFormula, symbols, supportCache);
  addFormulaSymbols(cache.sourceProblem->bad, symbols, supportCache);
  for (BoolExpr* output : cache.sourceProblem->observedOutputExprs0) {
    addFormulaSymbols(output, symbols, supportCache);
  }
  for (BoolExpr* output : cache.sourceProblem->observedOutputExprs1) {
    addFormulaSymbols(output, symbols, supportCache);
  }
  for (BoolExpr* strictEquality :
       cache.sourceProblem->dualRailOutputStrictEqualityExprs) {
    addFormulaSymbols(strictEquality, symbols, supportCache);
  }
  addRelevantComplementedStatePartners(complementPartners, symbols);
  addRelevantSameFrameStateEqualityPartners(complementPartners, symbols);
  addRelevantDualRailPartners(complementPartners, symbols);
  symbols.insert(cache.sourceProblem->allSymbols.begin(),
                 cache.sourceProblem->allSymbols.end());

  cache.frameZeroPredecessorSymbols = sortUniqueSymbols(std::move(symbols));
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: shared exact F[0] query surface symbols=",
        cache.frameZeroPredecessorSymbols.size());
  }
}

std::vector<size_t> sortClosedCurrentFrameSymbols(
    const ComplementPartnerIndex& complementPartners,
    std::unordered_set<size_t> symbols) {
  addRelevantComplementedStatePartners(complementPartners, symbols);
  addRelevantSameFrameStateEqualityPartners(complementPartners, symbols);
  addRelevantDualRailPartners(complementPartners, symbols);
  return sortUniqueSymbols(std::move(symbols));
} // LCOV_EXCL_LINE

PredecessorFrameSymbolSurfaceKey makePredecessorFrameSymbolSurfaceKey(
    const KInductionProblem& problem,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const ComplementPartnerIndex& complementPartners,
    PdrFormulaSupportCache* supportCache) {
  PredecessorFrameSymbolSurfaceKey key;
  key.problem = &problem;
  key.initFormula = initFormula;
  key.frameInvariant = frameInvariant;
  key.complementPartners = &complementPartners;
  key.supportCache = supportCache;
  key.level = level;
  key.frameFingerprint = frameClausesFingerprint(frames, level);
  return key;
}

std::vector<size_t> buildStablePredecessorCurrentFrameSymbols(
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const ComplementPartnerIndex& complementPartners,
    PdrFormulaSupportCache* supportCache) {
  std::unordered_set<size_t> symbols;
  if (level == 0) {
    addFormulaSymbols(initFormula, symbols, supportCache);
    addAllFrameClauseSymbols(frames[0], symbols);
  } else {
    addFormulaSymbols(frameInvariant, symbols, supportCache); // LCOV_EXCL_LINE
    addAllFrameClauseSymbols(frames[level], symbols); // LCOV_EXCL_LINE
  }

  // The relation closures below are independent of the target cube. Closing
  // this stable frame side once is equivalent to closing it together with each
  // query's dynamic symbols, because the closures only add partner/equality
  // symbols and do not inspect SAT polarity or clause state.
  return sortClosedCurrentFrameSymbols(
      complementPartners, std::move(symbols));
}

const std::vector<size_t>& cachedStablePredecessorCurrentFrameSymbols(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const ComplementPartnerIndex& complementPartners,
    PdrFormulaSupportCache* supportCache) {
  const PredecessorFrameSymbolSurfaceKey key =
      makePredecessorFrameSymbolSurfaceKey(
          problem,
          initFormula,
          frameInvariant,
          frames,
          level,
          complementPartners,
          supportCache);
  auto& currentFrameSymbols = cache.currentFrameSymbolsByLevel[level];
  if (!currentFrameSymbols.valid ||
      !(currentFrameSymbols.key == key)) { // LCOV_EXCL_LINE
    currentFrameSymbols.symbols =
        buildStablePredecessorCurrentFrameSymbols(
            initFormula,
            frameInvariant,
            frames,
            level,
            complementPartners,
            supportCache);
    currentFrameSymbols.key = key;
    currentFrameSymbols.valid = true;
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: predecessor frame symbol cache built level=",
          level,
          " symbols=",
          currentFrameSymbols.symbols.size(),
          " frame_fingerprint=",
          key.frameFingerprint);
    }
  }
  return currentFrameSymbols.symbols;
}

std::vector<size_t> mergePredecessorSymbolAddition(
    std::vector<size_t> base,
    const std::vector<size_t>& addition) {
  if (addition.empty()) {
    return base;
  }
  return detail::mergeSortedPdrSymbolVectors(base, addition);
}

std::vector<size_t> predecessorCurrentFrameQuerySymbolsFromCachedSurface(
    const KInductionProblem& problem,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    bool excludeTargetOnCurrentFrame,
    const PredecessorTargetSurface& targetSurface,
    const ComplementPartnerIndex& complementPartners,
    PredecessorAssumptionCache& predecessorAssumptionCache,
    PdrFormulaSupportCache* supportCache) {
  const std::vector<size_t>& stableSymbols =
      cachedStablePredecessorCurrentFrameSymbols(
          predecessorAssumptionCache,
          problem,
          initFormula,
          frameInvariant,
          frames,
          level,
          complementPartners,
          supportCache);
  std::vector<size_t> merged = stableSymbols;
  merged = mergePredecessorSymbolAddition(
      std::move(merged), targetSurface.closedPredecessorSymbols);

  // Partner/equality closure distributes over set union. Closing the immutable
  // target and property supports once, then merging their sorted vectors,
  // produces the same exact symbol set as closing their union on every query.
  if (predecessorSourceFrameIsKnownSafe(level)) {
    if (supportCache != nullptr) {
      merged = mergePredecessorSymbolAddition(
          std::move(merged),
          supportCache->relationClosedSupport(
              problem.property, complementPartners));
    } else {
      std::unordered_set<size_t> propertySymbols;
      addFormulaSymbols(problem.property, propertySymbols, nullptr);
      merged = mergePredecessorSymbolAddition(
          std::move(merged),
          sortClosedCurrentFrameSymbols(
              complementPartners, std::move(propertySymbols)));
    }
  } // LCOV_EXCL_LINE
  merged = mergePredecessorSymbolAddition(
      std::move(merged), targetSurface.closedTransitionSupportSymbols);
  return excludeTargetOnCurrentFrame
             ? mergePredecessorSymbolAddition(
                   std::move(merged), targetSurface.targetSymbols)
             : merged;
}

std::vector<size_t> predecessorCurrentFrameQuerySymbols(
    const KInductionProblem& problem,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const StateCube& targetCube,
    bool excludeTargetOnCurrentFrame,
    const PredecessorTargetSurface& targetSurface,
    const ComplementPartnerIndex& complementPartners,
    PredecessorAssumptionCache* predecessorAssumptionCache,
    PdrFormulaSupportCache* supportCache) {
  if (predecessorAssumptionCache != nullptr) {
    // The cache key includes the complete frame identity and fingerprint, so
    // the exact stable symbol closure is reusable for every PDR problem size.
    return predecessorCurrentFrameQuerySymbolsFromCachedSurface(
        problem,
        initFormula,
        frameInvariant,
        frames,
        level,
        excludeTargetOnCurrentFrame,
        targetSurface,
        complementPartners,
        *predecessorAssumptionCache,
        supportCache);
  }

  std::unordered_set<size_t> symbols;
  symbols.reserve(
      targetSurface.predecessorSymbols.size() +
      targetSurface.transitionSupportSymbols.size() +
      (excludeTargetOnCurrentFrame ? targetCube.size() : 0));
  symbols.insert(targetSurface.predecessorSymbols.begin(),
                 targetSurface.predecessorSymbols.end());
  addFrameConstraintSymbols(
      initFormula,
      frameInvariant,
      frames,
      level,
      complementPartners,
      symbols,
      supportCache);
  if (predecessorSourceFrameIsKnownSafe(level)) {
    // The safe-frame property is encoded below, so include its support in the
    // exact query surface.
    addFormulaSymbols(problem.property, symbols, supportCache);
  }
  for (const auto symbol : targetSurface.transitionSupportSymbols) {
    if (symbol >= 2) {
      symbols.insert(symbol);
    }
  }
  addRelevantComplementedStatePartners(complementPartners, symbols);
  addRelevantSameFrameStateEqualityPartners(complementPartners, symbols);
  addRelevantDualRailPartners(complementPartners, symbols);
  if (excludeTargetOnCurrentFrame) {
    addCubeSymbols(targetCube, symbols);
  }
  return sortUniqueSymbols(std::move(symbols));
}

const std::vector<size_t>& predecessorAssumptionCacheSymbols(
    const TransitionExprResolver& transitionByState,
    size_t level,
    const std::vector<size_t>& solverSymbols,
    PredecessorAssumptionCache* cache) {
  if (cache == nullptr) {
    return solverSymbols;
  }
  if (level == 0 &&
      cache->sharedFrameZeroPredecessorSymbols != nullptr) {
    if (cache->sharedFrameZeroPredecessorSymbols != &solverSymbols) {
      detail::widenSortedPdrSymbolSurface(
          *cache->sharedFrameZeroPredecessorSymbols, solverSymbols);
    }
    return *cache->sharedFrameZeroPredecessorSymbols;
  }

  // Section V uses one incremental SAT instance per frame. Preserve each
  // frame's monotonic surface when obligations move to another level and back.
  bool surfaceWidened = false;
  const auto& stableSurface = cache->predecessorSolverSymbolSurfaces.widen(
      &transitionByState, level, solverSymbols, &surfaceWidened);
  if (surfaceWidened) {
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: predecessor cached solver surface widened symbols=",
          stableSurface.size(),
          " requested=",
          solverSymbols.size());
    }
  }
  return stableSurface;
}

std::vector<size_t> initIntersectionSymbols(const KInductionProblem& problem,
                                            BoolExpr* initFormula) {
  // One incremental solver serves every candidate cube in this PDR run, so its
  // symbol surface includes every state bit that a later cube may assume.
  std::unordered_set<size_t> symbols;
  addFormulaSymbols(initFormula, symbols);
  const auto stateSymbols = problem.combinedStateSymbols();
  symbols.insert(stateSymbols.begin(), stateSymbols.end());
  // All relation endpoints are state symbols and are already present above.
  return sortUniqueSymbols(std::move(symbols));
}

std::optional<bool> findCubeLiteralValue(const StateCube& cube, size_t symbol) {
  const auto it = std::lower_bound(
      cube.begin(),
      cube.end(),
      symbol,
      [](const CubeLiteral& literal, size_t requestedSymbol) {
        return literal.symbol < requestedSymbol;
      // LCOV_EXCL_START
      });
      // LCOV_EXCL_STOP
  if (it == cube.end() || it->symbol != symbol) {
    return std::nullopt;
  }
  return it->value;
}

bool contradictsAssignments(
    const StateCube& cube,
    const std::vector<std::pair<size_t, bool>>& initAssignments) {
  for (const auto& [symbol, value] : initAssignments) {
    if (const auto cubeValue = findCubeLiteralValue(cube, symbol);
        cubeValue.has_value() && *cubeValue != value) {
      // LCOV_EXCL_START
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
  }
  return false;
}

bool contradictsIndexedAssignments(const StateCube& cube,
                                    const InitFactIndex& initFacts) {
  for (const auto& literal : cube) {
    const auto assignment =
        initFacts.assignmentValueMasks.find(literal.symbol);
    if (assignment == initFacts.assignmentValueMasks.end()) {
      continue;
    }
    const unsigned char oppositeValueMask =
        static_cast<unsigned char>(literal.value ? 1U : 2U);
    if ((assignment->second & oppositeValueMask) != 0) {
      return true;
    }
  }
  return false;
}

bool contradictsComplements(
    const StateCube& cube,
    const std::vector<std::pair<size_t, size_t>>& complements) {
  for (const auto& [primarySymbol, complementedSymbol] : complements) {
    const auto primaryValue = findCubeLiteralValue(cube, primarySymbol);  // LCOV_EXCL_LINE
    const auto complementedValue = findCubeLiteralValue(cube, complementedSymbol);  // LCOV_EXCL_LINE
    if (primaryValue.has_value() && complementedValue.has_value() &&  // LCOV_EXCL_LINE
        *primaryValue == *complementedValue) {  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
  }
  return false;
}

void reservePdrTransitionEncodingVars(SATSolverWrapper& solver,
                                      size_t estimatedNodes) {
  if (estimatedNodes < kMinPdrTransitionSolverReserveNodes) {
    return;
  }
  solver.reserveAdditionalVars( // LCOV_EXCL_LINE
      std::min(estimatedNodes, kMaxPdrTransitionSolverReserveHint)); // LCOV_EXCL_LINE
}

bool cubeContradictsKnownInitFacts(
    const KInductionProblem& problem,
    const StateCube& cube,
    const InitFactIndex* initFacts) {
  const bool usesBootstrapFrontier = problem.resetBootstrapCycles != 0;
  if (!usesBootstrapFrontier) {
    const bool contradictsInit =
        initFacts != nullptr
            ? contradictsIndexedAssignments(cube, *initFacts)
            : contradictsAssignments(cube, problem.initialStateAssignments);
    if (contradictsInit) {
      return true;
    }
  }
  if (problem.complementedStatePairs0.size() <=
      kMaxComplementPairsForCheapInitCheck &&
      contradictsComplements(cube, problem.complementedStatePairs0)) {
    return true;
  }
  if (problem.complementedStatePairs1.size() <=
      kMaxComplementPairsForCheapInitCheck &&
      contradictsComplements(cube, problem.complementedStatePairs1)) {
    return true;
  }
  return false;
}


void addTransitionConstraintsForTargetGroups(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    // LCOV_EXCL_START
    const TransitionExprResolver& transitionByState,
    size_t frame,
    const std::vector<TransitionEncodingLiteralGroup>& groups,
    const std::vector<size_t>& supportSymbols,
    std::unordered_map<size_t, int>* encodedLeafLits = nullptr) {
  // LCOV_EXCL_STOP
  for (const auto& group : groups) {
    std::unordered_map<size_t, int> leafLits =
        variables.makeLeafLits(frame, supportSymbols);
    const size_t estimatedNodes =
        estimateTransitionEncodingNodes(transitionByState, group.stateSymbols);
    reservePdrTransitionEncodingVars(solver, estimatedNodes);
    FrameFormulaEncoder encoder(
        solver,
        std::move(leafLits),
        // LCOV_EXCL_START
        group.symbolMap,
        false,
        estimatedNodes);
    for (const auto& literal : group.literals) {
      const TransitionExprView view =
      // LCOV_EXCL_STOP
          transitionByState.expressionView(literal.transitionSymbol);
      if (view.symbolMap != group.symbolMap) {
        throw std::runtime_error("Inconsistent transition symbol map");  // LCOV_EXCL_LINE
      }
      const int transitionLit = encoder.encode(
          view.expr,
          transitionByState.encodingPostorder(literal.transitionSymbol));
      solver.addClause({literal.desiredValue ? transitionLit : -transitionLit});
    }
    if (encodedLeafLits != nullptr) {
      const auto& groupLeafLits = encoder.leafLits();
      encodedLeafLits->insert(groupLeafLits.begin(), groupLeafLits.end());
    }
  }
}

FrameFormulaEncoder& cachedPredecessorTransitionEncoder(
    PredecessorAssumptionSolver& cachedSolver,
    const std::unordered_map<size_t, size_t>* symbolMap,
    size_t frame,
    size_t estimatedNodes) {
  const auto existing =
      cachedSolver.transitionEncoderBySymbolMap.find(symbolMap);
  if (existing != cachedSolver.transitionEncoderBySymbolMap.end()) {
    return *existing->second;
  }

  // Use the cached solver's complete symbol surface for this encoder. It is
  // built once per reusable predecessor solver, and it prevents a later target
  // in the same surface from missing a leaf that was outside the first target's
  // transition support slice.
  auto encoder = std::make_unique<FrameFormulaEncoder>(
      *cachedSolver.solver,
      cachedSolver.variables->makeLeafLits(frame),
      symbolMap,
      false,
      estimatedNodes);
  cachedSolver.transitionLeafLits.insert(
      encoder->leafLits().begin(), encoder->leafLits().end());
  auto [inserted, insertedNew] =
      cachedSolver.transitionEncoderBySymbolMap.emplace(
          symbolMap, std::move(encoder));
  (void)insertedNew;
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: predecessor transition encoder cached symbols=",
        inserted->second->leafLits().size(),
        " estimated_nodes=",
        estimatedNodes);
  }
  return *inserted->second;
}

const PreparedPredecessorTargetAssumptions&
PredecessorAssumptionSolver::prepareTargetAssumptions(
    const TransitionExprResolver& transitionByState,
    size_t frame,
    const StateClause& targetIdentity,
    const std::vector<TransitionEncodingLiteralGroup>& groups) {
  const auto cachedTarget = preparedTargetAssumptions.find(targetIdentity);
  if (cachedTarget != preparedTargetAssumptions.end()) {
    ++preparedTargetAssumptionHits;
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: predecessor target assumptions reused hits=",
          preparedTargetAssumptionHits,
          " entries=",
          preparedTargetAssumptions.size(),
          " assumptions=",
          cachedTarget->second.assumptions.size());
    }
    return cachedTarget->second;
  }

  size_t assumptionCount = 0;
  for (const auto& group : groups) {
    assumptionCount += group.literals.size();
  }
  PreparedPredecessorTargetAssumptions prepared;
  prepared.assumptions.reserve(assumptionCount);
  prepared.literalByAssumption.reserve(assumptionCount);
  for (const auto& group : groups) {
    FrameFormulaEncoder* encoder = nullptr;
    for (const auto& literal : group.literals) {
      const TransitionAssumptionKey key{
          literal.transitionSymbol,
          literal.desiredValue};
      const auto cachedIt =
          assumptionByTransitionLiteral.find(key);
      if (cachedIt != assumptionByTransitionLiteral.end()) {
        prepared.literalByAssumption.emplace(
            cachedIt->second, literal.originalLiteral);
        prepared.assumptions.push_back(cachedIt->second);
        continue;
      }

      if (encoder == nullptr) {
        const size_t estimatedNodes =
            estimateTransitionEncodingNodes(
                transitionByState, group.stateSymbols);
        reservePdrTransitionEncodingVars(*solver, estimatedNodes);
        encoder = &cachedPredecessorTransitionEncoder(
            *this,
            group.symbolMap,
            frame,
            estimatedNodes);
      }
      const TransitionExprView view =
          transitionByState.expressionView(literal.transitionSymbol);
      if (view.symbolMap != group.symbolMap) {
        throw std::runtime_error("Inconsistent transition symbol map");  // LCOV_EXCL_LINE
      }
      const int transitionLit = encoder->encode(
          view.expr,
          transitionByState.encodingPostorder(literal.transitionSymbol));
      // Store both polarities once the transition root is encoded. Neighboring
      // PDR cubes often ask for the opposite value of the same next-state bit;
      // reusing the root literal avoids rebuilding the same transition cone.
      assumptionByTransitionLiteral.emplace(
          TransitionAssumptionKey{literal.transitionSymbol, true},
          transitionLit);
      assumptionByTransitionLiteral.emplace(
          TransitionAssumptionKey{literal.transitionSymbol, false},
          -transitionLit);
      const int assumptionLit =
          literal.desiredValue ? transitionLit : -transitionLit;
      prepared.literalByAssumption.emplace(
          assumptionLit, literal.originalLiteral);
      prepared.assumptions.push_back(assumptionLit);
    }
  }
  auto [inserted, insertedNew] = preparedTargetAssumptions.emplace(
      targetIdentity, std::move(prepared));
  (void)insertedNew;
  return inserted->second;
}

StateCube failedAssumptionCubeFromTargetContext(
    const SATSolverWrapper& solver,
    const PreparedPredecessorTargetAssumptions& targetContext) {
  StateCube core;
  for (const int failedLit : solver.failedAssumptions()) {
    const auto literalIt =
        targetContext.literalByAssumption.find(failedLit);
    if (literalIt == targetContext.literalByAssumption.end()) {
      continue;
    }
    core.push_back(literalIt->second);
  }
  normalizeCube(core);
  return core;
}

StateCube cachedPredecessorUnsatCoreFromTargetContext(
    SATSolverWrapper& solver,
    const PreparedPredecessorTargetAssumptions& targetContext) {
  // Section V takes the failed target assumptions directly from the one
  // incremental solver. Figure 7 performs any further reduction explicitly.
  return failedAssumptionCubeFromTargetContext(solver, targetContext);
}

int PredecessorAssumptionSolver::q2SelectorFor(
    const StateClause& exclusionClause,
    size_t frame,
    bool blockingQuery) {
  const auto cachedIt =
      q2SelectorByExclusionClause.find(exclusionClause);
  if (cachedIt != q2SelectorByExclusionClause.end()) {
    auto& entry = cachedIt->second;
    if (blockingQuery && !entry.blockingQuery) {
      q2StatusSelectorRecency.erase(entry.recency);
      q2BlockingSelectorRecency.push_front(&cachedIt->first);
      entry.recency = q2BlockingSelectorRecency.begin();
      entry.blockingQuery = true;
    } else {
      auto& recency = entry.blockingQuery
                          ? q2BlockingSelectorRecency
                          : q2StatusSelectorRecency;
      recency.splice(recency.begin(), recency, entry.recency);
    }
    ++q2SelectorReuseCount;
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: Q2 selector cache reused count=",
          q2SelectorReuseCount,
          " cube_literals=",
          exclusionClause.size());
    }
    return entry.selector;
  }

  // SAT literals reserve 0/1 for constants; raw solver variable indices do
  // not. This selector guards only the exact temporary Q2 clause.
  const int selector = solver->newVar() + 2;
  std::vector<int> satClause;
  satClause.reserve(exclusionClause.size() + 1);
  satClause.push_back(-selector);
  for (const auto& literal : exclusionClause) {
    if (!variables->hasSymbol(literal.symbol)) {
      throw std::runtime_error( // LCOV_EXCL_LINE
          "PDR cached negated-cube encoding missing symbol " + // LCOV_EXCL_LINE
          std::to_string(literal.symbol) + " at frame " + // LCOV_EXCL_LINE
          std::to_string(frame) + " in cube of size " + // LCOV_EXCL_LINE
          std::to_string(exclusionClause.size())); // LCOV_EXCL_LINE
    }
    const int satLiteral =
        variables->getLiteral(literal.symbol, frame);
    satClause.push_back(literal.positive ? satLiteral : -satLiteral);
  }
  solver->addClause(satClause);

  auto [inserted, insertedNew] = q2SelectorByExclusionClause.emplace(
      exclusionClause, Q2SelectorCacheEntry{selector, blockingQuery, {}});
  (void)insertedNew;
  auto& insertedRecency = blockingQuery
                              ? q2BlockingSelectorRecency
                              : q2StatusSelectorRecency;
  insertedRecency.push_front(&inserted->first);
  inserted->second.recency = insertedRecency.begin();

  // Each state symbol has two literal polarities. Retaining at most one exact
  // target context per possible literal keeps this accelerator linear while
  // preserving the small recurring cubes that dominate blocking queries.
  const size_t cacheLimit =
      2 * std::max<size_t>(key.solverSymbols.size(), 1);
  while (q2SelectorByExclusionClause.size() > cacheLimit) {
    // Figure 6 blocking targets recur while obligations move through frames.
    // Prefer retiring least-recently-used status-only targets from Figures 7
    // and 9; this changes cache retention only, never the exact SAT query.
    const bool retireStatusOnly = !q2StatusSelectorRecency.empty();
    auto& retiredRecency = retireStatusOnly
                               ? q2StatusSelectorRecency
                               : q2BlockingSelectorRecency;
    const StateClause* retiredClause = retiredRecency.back();
    const auto retired = q2SelectorByExclusionClause.find(*retiredClause);
    solver->addClause({-retired->second.selector});
    retiredRecency.pop_back();
    q2SelectorByExclusionClause.erase(retired);
    ++q2SelectorEvictionCount;
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: Q2 selector cache evicted count=",
          q2SelectorEvictionCount,
          " retained=",
          q2SelectorByExclusionClause.size(),
          " limit=",
          cacheLimit,
          " class=",
          retireStatusOnly ? "status" : "blocking");
    }
  }
  return selector;
}

size_t PredecessorAssumptionSolver::retireStatusQ2Selectors() {
  size_t retiredCount = 0;
  while (!q2StatusSelectorRecency.empty()) {
    const StateClause* retiredClause = q2StatusSelectorRecency.back();
    const auto retired = q2SelectorByExclusionClause.find(*retiredClause);
    solver->addClause({-retired->second.selector});
    q2StatusSelectorRecency.pop_back();
    q2SelectorByExclusionClause.erase(retired);
    ++retiredCount;
  }
  q2SelectorEvictionCount += retiredCount;
  if (retiredCount != 0 && shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: Q2 status selectors retired after generalization count=",
        retiredCount,
        " retained_blocking=",
        q2BlockingSelectorRecency.size());
  }
  return retiredCount;
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP



// LCOV_DISABLED_STOP

void normalizeCube(StateCube& cube) {
  // Canonical ordering lets us compare cubes structurally and avoid learning
  // the same obligation more than once with a different literal order.
  if (!std::is_sorted(cube.begin(), cube.end(), cubeLiteralLess)) {
    std::sort(cube.begin(), cube.end(), cubeLiteralLess);
  }
  cube.erase(std::unique(cube.begin(), cube.end()), cube.end());
}

void normalizeClause(StateClause& clause) {
  // Clauses are canonicalized for the same reason: later subsumption and
  // LCOV_DISABLED_START
  // convergence checks depend on stable ordering and deduplication.
  if (!std::is_sorted(clause.begin(), clause.end(), clauseLiteralLess)) {
    std::sort(clause.begin(), clause.end(), clauseLiteralLess);
  }
  // LCOV_DISABLED_STOP
  clause.erase(std::unique(clause.begin(), clause.end()), clause.end());
}

SymbolPair canonicalPair(size_t lhs, size_t rhs) {
  if (rhs < lhs) {
    std::swap(lhs, rhs);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return SymbolPair{lhs, rhs};
}

InitFactIndex buildInitFactIndex(const KInductionProblem& problem) {
  const bool usesBootstrapFrontier = problem.resetBootstrapCycles != 0;
  InitFactIndex index;
  if (!usesBootstrapFrontier) {
    index.assignments.reserve(problem.initialStateAssignments.size());
    index.assignmentValueMasks.reserve(
        problem.initialStateAssignments.size());
    for (const auto& [symbol, value] : problem.initialStateAssignments) {
      index.assignments.emplace(symbol, value);
      index.assignmentValueMasks[symbol] |=
          static_cast<unsigned char>(value ? 2U : 1U);
      index.relations.ensureSymbol(symbol);
    }
  }
  for (const auto& [lhsSymbol, rhsSymbol] :
       problem.sameFrameStateEqualityPairs0) {
    index.equalities.insert(canonicalPair(lhsSymbol, rhsSymbol));
    index.relations.addEquality(lhsSymbol, rhsSymbol);
  }
  for (const auto& [lhsSymbol, rhsSymbol] :
       problem.sameFrameStateEqualityPairs1) {
    index.equalities.insert(canonicalPair(lhsSymbol, rhsSymbol));
    index.relations.addEquality(lhsSymbol, rhsSymbol);
  }
  index.complements.reserve(
      problem.complementedStatePairs0.size() +
      problem.complementedStatePairs1.size());
  for (const auto& [primarySymbol, complementedSymbol] :
       // LCOV_DISABLED_START
       problem.complementedStatePairs0) {
       // LCOV_DISABLED_STOP
    index.complements.insert(canonicalPair(primarySymbol, complementedSymbol));
    index.relations.addComplement(primarySymbol, complementedSymbol);
  }
  for (const auto& [primarySymbol, complementedSymbol] :
       problem.complementedStatePairs1) {
    index.complements.insert(canonicalPair(primarySymbol, complementedSymbol));
    index.relations.addComplement(primarySymbol, complementedSymbol);
  }
  index.rootAssignments.reserve(index.assignments.size());
  std::vector<std::pair<size_t, bool>> orderedAssignments(
      index.assignments.begin(), index.assignments.end());
  std::sort(orderedAssignments.begin(), orderedAssignments.end());
  for (const auto& [symbol, value] : orderedAssignments) {
    const auto root = index.relations.findWithParity(symbol);
    if (!root.has_value()) {
      continue;  // LCOV_EXCL_LINE
    }
    const bool rootValue = value ^ root->second;
    if (const auto it = index.rootAssignments.find(root->first);
        it == index.rootAssignments.end()) {
      index.rootAssignments.emplace(root->first, rootValue);
    }
  }
  return index;
}

std::optional<StateCube> knownInitConflictCube(const InitFactIndex& facts,
                                               const StateCube& cube) {
  // PDR frequently reaches a level-0 cube that is impossible only because it
  // violates a startup equality such as "state0 == state1".  Learning the full
  // LCOV_DISABLED_START
  // 100+ literal cube makes the engine enumerate many adjacent impossible
  // LCOV_DISABLED_STOP
  // startup states.  This extractor turns the visible conflict into the
  // smallest safe cube:
  // LCOV_DISABLED_START
  //   - one literal for an init assignment conflict;
  //   - two literals for equality/complement conflicts.
  // The learned clause is still exactly an Init consequence, but much stronger.
  std::unordered_map<size_t, std::pair<bool, CubeLiteral>> cubeValueByRoot;
  // LCOV_DISABLED_STOP
  cubeValueByRoot.reserve(cube.size());
  for (const auto& literal : cube) {
    const auto root = facts.relations.findWithParity(literal.symbol);
    if (!root.has_value()) {
      const auto assignment = facts.assignments.find(literal.symbol);
      // LCOV_DISABLED_START
      if (assignment == facts.assignments.end() ||
          assignment->second == literal.value) {  // LCOV_EXCL_LINE
        continue;
      }
      // LCOV_DISABLED_STOP
      StateCube conflict{literal};  // LCOV_EXCL_LINE
      normalizeCube(conflict);  // LCOV_EXCL_LINE
      return conflict;  // LCOV_EXCL_LINE
    // LCOV_DISABLED_START
    }  // LCOV_EXCL_LINE

    const bool rootValue = literal.value ^ root->second;
    const auto assignment = facts.rootAssignments.find(root->first);
    if (assignment != facts.rootAssignments.end() &&
        assignment->second != rootValue) {
        // LCOV_DISABLED_STOP
      StateCube conflict{literal};  // LCOV_EXCL_LINE
      normalizeCube(conflict);  // LCOV_EXCL_LINE
      return conflict;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE

    if (const auto it = cubeValueByRoot.find(root->first);
        it != cubeValueByRoot.end()) {
      if (it->second.first != rootValue) {  // LCOV_EXCL_LINE
        StateCube conflict{it->second.second, literal};  // LCOV_EXCL_LINE
        normalizeCube(conflict);  // LCOV_EXCL_LINE
        return conflict;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_DISABLED_START
    cubeValueByRoot.emplace(root->first, std::pair{rootValue, literal});
  }
  // LCOV_DISABLED_STOP

  return std::nullopt;
}

// LCOV_DISABLED_START


StateClause clauseFromCube(const StateCube& cube) {
  StateClause clause;
  clause.reserve(cube.size());
  for (const auto& literal : cube) {
    clause.push_back({literal.symbol, !literal.value});
  }
  normalizeClause(clause);
  return clause;
}

StateCube cubeFromClauseNegation(const StateClause& clause) {
  StateCube cube;
  cube.reserve(clause.size());
  for (const auto& literal : clause) {
    cube.push_back({literal.symbol, !literal.positive});
  }
  normalizeCube(cube);
  return cube;
}

bool clauseSubsumes(const StateClause& lhs, const StateClause& rhs) {
  return std::includes(rhs.begin(), rhs.end(), lhs.begin(), lhs.end(),
                       [](const ClauseLiteral& a, const ClauseLiteral& b) {
                         if (a.symbol != b.symbol) {
                           return a.symbol < b.symbol;
                         }
                         return a.positive < b.positive;
                       });
}

bool frameHasSubsumingClause(const FrameClauses& frame,
                             const StateClause& clause) {
  for (const auto& existingClause : frame.clauses) {
    // Frames are sorted by clause size first. A larger clause cannot subsume
    // this candidate, so the remaining suffix cannot contain a match either.
    if (existingClause.size() > clause.size()) {
      break;
    }
    if (clauseSubsumes(existingClause, clause)) {
      return true;
    }
  }
  return false;
}


bool addClauseToFrame(FrameClauses& frame, StateClause clause) {
  normalizeClause(clause);
  if (frameHasSubsumingClause(frame, clause)) {
    return false;
  }

  // Keep each frame minimal so later SAT queries do not carry redundant facts.
  frame.clauses.erase(
      std::remove_if(
          frame.clauses.begin(),
          frame.clauses.end(),
          [&](const StateClause& existingClause) {
            return clauseSubsumes(clause, existingClause);
          }),
      frame.clauses.end());
  frame.addedClauseLog.push_back(clause);
  // The remaining clauses stay sorted after erase(), so a lower_bound insert
  // preserves the deterministic frame order without resorting the whole frame
  // for every learned clause.
  auto insertPosition =
      std::lower_bound(frame.clauses.begin(), frame.clauses.end(), clause,
                       stateClauseLess);
  frame.clauses.insert(insertPosition, std::move(clause));
  frame.clauseFingerprint.reset();
  return true;
}

bool addClauseToFrames(std::vector<FrameClauses>& frames,
                       const StateClause& clause,
                       size_t maxLevel) {
  bool addedAny = false;
  for (size_t level = 1; level <= maxLevel; ++level) {
    addedAny = addClauseToFrame(frames[level], clause) || addedAny;
  }
  return addedAny;
}  // LCOV_EXCL_LINE

bool hasExactClause(const std::vector<StateClause>& clauses,
                    const StateClause& clause) {
  const auto position =
      std::lower_bound(clauses.begin(), clauses.end(), clause, stateClauseLess);
  return position != clauses.end() && *position == clause;
}

bool addReusableInvariantCandidate(
    std::vector<StateClause>& candidates,
    StateClause clause) {
  normalizeClause(clause);
  const auto position = std::lower_bound(
      candidates.begin(), candidates.end(), clause, stateClauseLess);
  if (position != candidates.end() && *position == clause) {
    return false;
  }
  candidates.insert(position, std::move(clause));
  return true;
}

size_t injectReusableInvariantClauses(
    const PDRExactInitCache::Impl& cache,
    FrameClauses& target) {
  size_t imported = 0;
  for (const StateClause& clause : cache.reusableInvariant.clauses) {
    if (addClauseToFrame(target, clause)) {
      ++imported;
    }
  }
  if (imported != 0 && pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: reusable invariant clauses injected count=",
        imported,
        " retained=",
        cache.reusableInvariant.clauses.size());
  }
  return imported;
}

void storeReusableInvariantCandidates(
    PDRExactInitCache::Impl& cache,
    const std::vector<FrameClauses>& frames) {
  if (cache.reusableInvariantCertificationDisabled) {
    return;
  }
  std::vector<StateClause> additions;
  size_t additionalLiterals = 0;
  for (size_t level = 1; level < frames.size(); ++level) {
    for (const StateClause& clause : frames[level].clauses) {
      if (hasExactClause(cache.reusableInvariantCandidates, clause) ||
          hasExactClause(additions, clause)) {
        continue;
      }
      const auto position = std::lower_bound(
          additions.begin(), additions.end(), clause, stateClauseLess);
      additions.insert(position, clause);
      additionalLiterals += clause.size();
    }
  }
  if (additions.empty()) {
    return;
  }
  if (cache.reusableInvariantCandidateLiteralCount + additionalLiterals >
      kMaxReusableInvariantCandidateLiterals) {
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: reusable invariant candidates skipped literals=",
          additionalLiterals,
          " retained_literal_budget=",
          cache.reusableInvariantCandidateLiteralCount);
    }
    return;
  }

  size_t stored = 0;
  for (StateClause& clause : additions) {
    if (addReusableInvariantCandidate(
            cache.reusableInvariantCandidates, std::move(clause))) {
      ++stored;
    }
  }
  cache.reusableInvariantCandidateLiteralCount += additionalLiterals;
  ++cache.reusableInvariantCandidateRevision;
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: reusable invariant candidates stored count=",
        stored,
        " total=",
        cache.reusableInvariantCandidates.size(),
        " literal_budget_used=",
        cache.reusableInvariantCandidateLiteralCount);
  }
}

class ReusableInvariantCandidateRecorder {
 public:
  ReusableInvariantCandidateRecorder(
      PDRExactInitCache::Impl* cache,
      const std::vector<FrameClauses>& frames)
      : cache_(cache), frames_(frames) {}

  ~ReusableInvariantCandidateRecorder() {
    if (cache_ == nullptr) {
      return;
    }
    // Candidate reuse is optional. Allocation failure must not change PDR's
    // property verdict or discard an already certified invariant.
    try {
      storeReusableInvariantCandidates(*cache_, frames_);
    } catch (...) {  // LCOV_EXCL_LINE
    }
  }

 private:
  PDRExactInitCache::Impl* cache_ = nullptr;
  const std::vector<FrameClauses>& frames_;
};

void addStateClause(SATSolverWrapper& solver,
                    const FrameVariableStore& variables,
                    const StateClause& clause,
                    size_t frame) {
  std::vector<int> satClause;
  satClause.reserve(clause.size());
  for (const auto& literal : clause) {
    if (!variables.hasSymbol(literal.symbol)) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR frame-clause encoding missing symbol " +  // LCOV_EXCL_LINE
          std::to_string(literal.symbol) + " at frame " +  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          std::to_string(frame) + " in clause of size " +  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
          std::to_string(clause.size()));  // LCOV_EXCL_LINE
    }
    const int satLiteral = variables.getLiteral(literal.symbol, frame);
    satClause.push_back(literal.positive ? satLiteral : -satLiteral);
  }
  solver.addClause(satClause);
}

void addGuardedStateClause(SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const StateClause& clause,
                           size_t frame,
                           int activationLiteral) {
  std::vector<int> satClause;
  satClause.reserve(clause.size() + 1);
  satClause.push_back(-activationLiteral);
  for (const auto& literal : clause) {
    if (!variables.hasSymbol(literal.symbol)) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR guarded frame-clause encoding missing symbol " +  // LCOV_EXCL_LINE
          std::to_string(literal.symbol));  // LCOV_EXCL_LINE
    }
    const int satLiteral = variables.getLiteral(literal.symbol, frame);
    satClause.push_back(literal.positive ? satLiteral : -satLiteral);
  }
  solver.addClause(satClause);
}

bool clauseCoveredByVariables(const FrameVariableStore& variables,
                              const StateClause& clause) {
  for (const auto& literal : clause) {
    if (!variables.hasSymbol(literal.symbol)) {
      return false;  // LCOV_EXCL_LINE
    }
  }
  // LCOV_EXCL_START
  return true;
}



void addAllFrameClauses(SATSolverWrapper& solver,
                        const FrameVariableStore& variables,
                        const FrameClauses& frameClauses,
                        size_t frame) {
  // Every PDR query sees the complete learned frame.
  for (const auto& clause : frameClauses.clauses) {
    // LCOV_EXCL_START
    if (!clauseCoveredByVariables(variables, clause)) {
      continue;  // LCOV_EXCL_LINE
    }
    addStateClause(solver, variables, clause, frame);
  }
  // LCOV_EXCL_STOP
}

void addCubeAssumptions(SATSolverWrapper& solver,
                        const FrameVariableStore& variables,
                        const StateCube& cube,
                        size_t frame) {
  for (const auto& literal : cube) {
    if (!variables.hasSymbol(literal.symbol)) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR cube-assumption encoding missing symbol " +  // LCOV_EXCL_LINE
          std::to_string(literal.symbol) + " at frame " +  // LCOV_EXCL_LINE
          std::to_string(frame) + " in cube of size " +  // LCOV_EXCL_LINE
          std::to_string(cube.size()));  // LCOV_EXCL_LINE
    }
    solver.addClause(
        // LCOV_EXCL_START
        {literal.value ? variables.getLiteral(literal.symbol, frame)
                       : -variables.getLiteral(literal.symbol, frame)});
  }
}


// LCOV_EXCL_STOP
void addNegatedCubeClause(SATSolverWrapper& solver,
                          const FrameVariableStore& variables,
                          const StateCube& cube,
                          size_t frame) {
  std::vector<int> satClause;
  satClause.reserve(cube.size());
  for (const auto& literal : cube) {
    if (!variables.hasSymbol(literal.symbol)) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR negated-cube encoding missing symbol " +  // LCOV_EXCL_LINE
          std::to_string(literal.symbol) + " at frame " +  // LCOV_EXCL_LINE
          std::to_string(frame) + " in cube of size " +  // LCOV_EXCL_LINE
          std::to_string(cube.size()));  // LCOV_EXCL_LINE
    }
    const int satLiteral = variables.getLiteral(literal.symbol, frame);
    satClause.push_back(literal.value ? -satLiteral : satLiteral);
  }
  solver.addClause(satClause);
}

void addPostBootstrapResetInputConstraints(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t frame) {
  if (problem.resetBootstrapCycles == 0) {
    return;
  }

  // PDR frames are already positioned after the concrete reset prefix. The
  // reset controls are therefore no longer free environment inputs in one-step
  // predecessor queries; they must stay at their deasserted value on every PDR
  // transition, exactly as the concrete base solver constrains them.
  for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
    if (!variables.hasSymbol(symbol)) {
      continue;
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
    solver.addClause(
        {assertedValue ? -variables.getLiteral(symbol, frame)
                       : variables.getLiteral(symbol, frame)});
  }
}

struct ReusableInvariantFinderResult {
  FrameClauses invariant;
  size_t initialRejected = 0;
  size_t inductiveRejected = 0;
  size_t inductiveQueries = 0;
};

// Chockler et al.'s FMCAD'11 invariant finder converts arbitrary saved IC3
// clauses into the maximum subset H satisfying Init => H and H /\ T => H'.
// Only H may cross property runs; unfinished frame levels remain candidates.
class ReusableInvariantFinder {
 public:
  ReusableInvariantFinder(PDRExactInitCache::Impl& cache,
                          BoolExpr* initFormula)
      : cache_(cache),
        problem_(*cache.sourceProblem),
        initFormula_(initFormula),
        baseInvariant_(cache.reusableInvariant) {
    for (const StateClause& candidate :
         cache.reusableInvariantCandidates) {
      if (!hasExactClause(baseInvariant_.clauses, candidate)) {
        candidates_.push_back(candidate);
      }
    }
    collectCandidateSymbols();
  }

  std::optional<ReusableInvariantFinderResult> run() {
    ReusableInvariantFinderResult result;
    result.invariant = baseInvariant_;
    if (candidates_.empty()) {
      return result;
    }

    const auto initiallyValid = findInitiallyValidCandidates();
    if (!initiallyValid.has_value()) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    active_ = *initiallyValid;
    result.initialRejected =
        static_cast<size_t>(std::count(active_.begin(), active_.end(), false));
    if (result.initialRejected == candidates_.size()) {
      return result;
    }

    if (!findMaximumInductiveSubset(result)) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    return result;
  }

 private:
  void collectCandidateSymbols() {
    for (const StateClause& clause : candidates_) {
      addClauseSymbols(clause, candidateSymbols_);
      addClauseSymbols(clause, currentConstraintSymbols_);
    }
    for (const StateClause& clause : baseInvariant_.clauses) {
      addClauseSymbols(clause, currentConstraintSymbols_);
    }
  }

  std::vector<size_t> initialQuerySymbols() {
    std::unordered_set<size_t> symbols = candidateSymbols_;
    addFormulaSymbols(
        initFormula_, symbols, cache_.formulaSupportCache.get());
    addRelevantComplementedStatePartners(cache_.stateRelations, symbols);
    addRelevantSameFrameStateEqualityPartners(cache_.stateRelations, symbols);
    addRelevantDualRailPartners(cache_.stateRelations, symbols);
    return sortUniqueSymbols(std::move(symbols));
  }

  std::optional<std::vector<bool>> findInitiallyValidCandidatesUsing(
      SATSolverWrapper& solver,
      const FrameVariableStore& variables,
      std::vector<bool> initiallyValid,
      const std::vector<size_t>& unresolvedCandidates) const {
    std::vector<int> assumptions;
    for (const size_t index : unresolvedCandidates) {
      assumptions.clear();
      assumptions.reserve(candidates_[index].size());
      for (const ClauseLiteral& literal : candidates_[index]) {
        const int stateLiteral =
            variables.getLiteral(literal.symbol, 0);
        // Init violates a clause only when every one of its literals is false.
        assumptions.push_back(
            literal.positive ? -stateLiteral : stateLiteral);
      }
      const SATSolverWrapper::SolveStatus status =
          solveCertificationQuery(solver, assumptions);
      if (status == SATSolverWrapper::SolveStatus::Unknown) {
        return std::nullopt;  // LCOV_EXCL_LINE
      }
      initiallyValid[index] =
          status == SATSolverWrapper::SolveStatus::Unsat;
    }
    return initiallyValid;
  }

  std::optional<std::vector<bool>> findInitiallyValidCandidates() {
    std::vector<bool> initiallyValid(candidates_.size(), false);
    std::vector<size_t> unresolvedCandidates;
    unresolvedCandidates.reserve(candidates_.size());
    for (size_t index = 0; index < candidates_.size(); ++index) {
      const StateCube violatingCube =
          cubeFromClauseNegation(candidates_[index]);
      const bool hasKnownInitConflict =
          cache_.initFacts.has_value()
              ? knownInitConflictCube(
                    *cache_.initFacts, violatingCube).has_value()
              : cubeContradictsKnownInitFacts(
                    problem_, violatingCube, nullptr);
      if (hasKnownInitConflict) {
        // This is the same exact Init-conflict shortcut used by Figure 6 cube
        // blocking. It proves Init => clause without entering SAT.
        initiallyValid[index] = true;
      } else {
        unresolvedCandidates.push_back(index);
      }
    }
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: reusable invariant initial facts resolved=",
          candidates_.size() - unresolvedCandidates.size(),
          " unresolved=",
          unresolvedCandidates.size());
    }
    if (unresolvedCandidates.empty()) {
      return initiallyValid;
    }

    if (cache_.frameZeroPredecessorSolver != nullptr) {
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: reusable invariant initial checks reused "
            "shared exact F[0] solver unresolved=",
            unresolvedCandidates.size());
      }
      return findInitiallyValidCandidatesUsing(
          *cache_.frameZeroPredecessorSolver->solver,
          *cache_.frameZeroPredecessorSolver->variables,
          std::move(initiallyValid),
          unresolvedCandidates);
    }

    // This fallback preserves standalone PDR use before a shared F[0] owner
    // exists. Normal SEC output batches always take the reuse path above.
    SATSolverWrapper solver(
        SATSolverWrapper::assumptionSolverTypeFor(cache_.solverType));
    const std::vector<size_t> querySymbols = initialQuerySymbols();
    solver.configureForSecPdrPersistentQuery(querySymbols.size());
    FrameVariableStore variables(solver, querySymbols, 1);
    cache_.stateRelations.addClauses(
        solver, variables, querySymbols, 1);
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(0));
    solver.addClause({encoder.encode(initFormula_)});
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: reusable invariant local init solver created "
          "symbols=",
          querySymbols.size());
    }
    return findInitiallyValidCandidatesUsing(
        solver,
        variables,
        std::move(initiallyValid),
        unresolvedCandidates);
  }

  std::vector<size_t> inductiveQuerySymbols(
      std::vector<size_t>& transitionTargets) {
    std::unordered_set<size_t> symbols = currentConstraintSymbols_;
    transitionTargets = expandTransitionTargets(
        problem_,
        sortUniqueSymbols(candidateSymbols_),
        cache_.transitionByState);
    symbols.insert(transitionTargets.begin(), transitionTargets.end());
    for (const size_t target : transitionTargets) {
      const std::set<size_t>& support =
          cache_.transitionByState.support(target);
      symbols.insert(support.begin(), support.end());
    }
    for (const auto& [symbol, /*assertedValue*/ _] :
         problem_.resetBootstrapInputs) {
      symbols.insert(symbol);
    }
    addRelevantComplementedStatePartners(cache_.stateRelations, symbols);
    addRelevantSameFrameStateEqualityPartners(cache_.stateRelations, symbols);
    addRelevantDualRailPartners(cache_.stateRelations, symbols);
    return sortUniqueSymbols(std::move(symbols));
  }

  struct TransitionTargetGroup {
    const std::unordered_map<size_t, size_t>* symbolMap = nullptr;
    std::vector<size_t> targets;
  };

  void appendTransitionTarget(
      std::vector<TransitionTargetGroup>& groups,
      const std::unordered_map<size_t, size_t>* symbolMap,
      size_t target) {
    for (TransitionTargetGroup& group : groups) {
      if (group.symbolMap == symbolMap) {
        group.targets.push_back(target);
        return;
      }
    }
    TransitionTargetGroup group;
    group.symbolMap = symbolMap;
    group.targets.push_back(target);
    groups.push_back(std::move(group));
  }

  void addTransitionEquations(
      SATSolverWrapper& solver,
      const FrameVariableStore& variables,
      const std::vector<size_t>& transitionTargets) {
    std::vector<TransitionTargetGroup> groups;
    groups.reserve(3);
    for (const size_t target : transitionTargets) {
      const TransitionExprView view =
          cache_.transitionByState.expressionView(target);
      appendTransitionTarget(groups, view.symbolMap, target);
    }

    for (const TransitionTargetGroup& group : groups) {
      size_t estimatedNodes = 0;
      for (const size_t target : group.targets) {
        estimatedNodes += cache_.transitionByState.nodeCount(target);
      }
      reservePdrTransitionEncodingVars(solver, estimatedNodes);
      FrameFormulaEncoder encoder(
          solver,
          variables.makeLeafLits(0),
          group.symbolMap,
          false,
          estimatedNodes);
      for (const size_t target : group.targets) {
        const TransitionExprView view =
            cache_.transitionByState.expressionView(target);
        if (view.symbolMap != group.symbolMap) {
          throw std::runtime_error(  // LCOV_EXCL_LINE
              "Inconsistent invariant-finder transition symbol map");
        }
        const int transitionLiteral = encoder.encode(
            view.expr,
            cache_.transitionByState.encodingPostorder(target));
        addLiteralEquivalence(
            solver,
            variables.getLiteral(target, 1),
            transitionLiteral);
      }
    }
  }

  void addCandidateSelectors(
      SATSolverWrapper& solver,
      const FrameVariableStore& variables,
      std::vector<int>& currentSelectors,
      std::vector<int>& nextHolds) {
    currentSelectors.reserve(candidates_.size());
    nextHolds.reserve(candidates_.size());
    std::vector<int> someNextClauseFails;
    someNextClauseFails.reserve(candidates_.size());
    for (const StateClause& clause : candidates_) {
      const int currentSelector = solver.newVar() + 2;
      const int nextClauseHolds = solver.newVar() + 2;
      currentSelectors.push_back(currentSelector);
      nextHolds.push_back(nextClauseHolds);
      addGuardedStateClause(
          solver, variables, clause, 0, currentSelector);
      // For every literal a' in c', add (y OR !a'). Thus c' => y,
      // and a model with y=false identifies a clause violated after one step.
      for (const ClauseLiteral& literal : clause) {
        const int nextLiteral =
            variables.getLiteral(literal.symbol, 1);
        solver.addClause(
            {nextClauseHolds,
             literal.positive ? -nextLiteral : nextLiteral});
      }
      someNextClauseFails.push_back(-nextClauseHolds);
    }
    solver.addClause(someNextClauseFails);
  }

  std::vector<int> activeAssumptions(
      const std::vector<int>& currentSelectors,
      const std::vector<int>& nextHolds) const {
    std::vector<int> assumptions;
    assumptions.reserve(candidates_.size() * 2);
    for (size_t index = 0; index < candidates_.size(); ++index) {
      if (active_[index]) {
        assumptions.push_back(currentSelectors[index]);
        continue;
      }
      // Disabled candidates must constrain neither the current conjunction nor
      // the disjunction that asks for an active next-state clause violation.
      assumptions.push_back(-currentSelectors[index]);
      assumptions.push_back(nextHolds[index]);
    }
    return assumptions;
  }

  size_t removeViolatedCandidates(
      const SATSolverWrapper& solver,
      const std::vector<int>& nextHolds) {
    size_t removed = 0;
    for (size_t index = 0; index < candidates_.size(); ++index) {
      if (active_[index] && !solver.getLiteralValue(nextHolds[index])) {
        active_[index] = false;
        ++removed;
      }
    }
    return removed;
  }

  int64_t certificationConflictLimit() const {
    return problem_.usesDualRailStateEncoding
               ? dualRailPredecessorConflictLimit(
                     PredecessorQueryPurpose::GeneralizeBlocker)
               : -1;
  }

  int64_t certificationDecisionLimit() const {
    return problem_.usesDualRailStateEncoding
               ? dualRailPredecessorDecisionLimit(
                     PredecessorQueryPurpose::GeneralizeBlocker)
               : -1;
  }

  int64_t certificationTickLimit() const {
    return problem_.usesDualRailStateEncoding
               ? kDualRailInvariantCertificationPerQueryTickLimit
               : -1;
  }

  SATSolverWrapper::SolveStatus solveCertificationQuery(
      SATSolverWrapper& solver,
      const std::vector<int>& assumptions) const {
    if (!cache_.reusableInvariantCertificationBudget.has_value()) {
      return solver.solveWithAssumptionsStatus(
          assumptions,
          certificationConflictLimit(),
          certificationDecisionLimit(),
          certificationTickLimit());
    }
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(
        *cache_.reusableInvariantCertificationBudget);
    return solver.solveWithAssumptionsStatus(
        assumptions,
        certificationConflictLimit(),
        certificationDecisionLimit(),
        certificationTickLimit());
  }

  void buildInvariant(ReusableInvariantFinderResult& result) const {
    for (size_t index = 0; index < candidates_.size(); ++index) {
      if (active_[index]) {
        addClauseToFrame(result.invariant, candidates_[index]);
      }
    }
    result.invariant.addedClauseLog.clear();
  }

  bool findMaximumInductiveSubset(
      ReusableInvariantFinderResult& result) {
    SATSolverWrapper solver(
        SATSolverWrapper::assumptionSolverTypeFor(cache_.solverType));
    std::vector<size_t> transitionTargets;
    const std::vector<size_t> querySymbols =
        inductiveQuerySymbols(transitionTargets);
    solver.configureForSecPdrPersistentQuery(querySymbols.size());
    FrameVariableStore variables(solver, querySymbols, 2);
    cache_.stateRelations.addClauses(
        solver, variables, querySymbols, 2);
    addPostBootstrapResetInputConstraints(
        solver, variables, problem_, 0);
    addTransitionEquations(solver, variables, transitionTargets);
    // H is already inductive, so it only constrains the current state. The
    // exact extension query needs next-state selectors for new candidates.
    for (const StateClause& clause : baseInvariant_.clauses) {
      addStateClause(solver, variables, clause, 0);
    }
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: reusable invariant inductive local solver created "
          "symbols=",
          querySymbols.size(),
          " transition_targets=",
          transitionTargets.size(),
          " conflict_limit=",
          certificationConflictLimit(),
          " decision_limit=",
          certificationDecisionLimit(),
          " tick_limit=",
          certificationTickLimit());
    }

    std::vector<int> currentSelectors;
    std::vector<int> nextHolds;
    addCandidateSelectors(
        solver, variables, currentSelectors, nextHolds);
    while (std::find(active_.begin(), active_.end(), true) != active_.end()) {
      const std::vector<int> assumptions =
          activeAssumptions(currentSelectors, nextHolds);
      ++result.inductiveQueries;
      const SATSolverWrapper::SolveStatus status =
          solveCertificationQuery(solver, assumptions);
      if (status == SATSolverWrapper::SolveStatus::Unknown) {
        const auto& budget = cache_.reusableInvariantCertificationBudget;
        if (pdrStatsEnabled()) {
          emitSecDiag(
              "SEC PDR stats: reusable invariant certification budget "
              "exhausted query=",
              result.inductiveQueries,
              " active_candidates=",
              std::count(active_.begin(), active_.end(), true),
              " conflict_limit=",
              certificationConflictLimit(),
              " decision_limit=",
              certificationDecisionLimit(),
              " tick_limit=",
              certificationTickLimit(),
              " total_ticks=",
              budget.has_value() ? budget->ticksUsed() : 0,
              "/",
              budget.has_value() ? budget->tickLimit() : 0);
        }
        // Candidate reuse is optional. UNKNOWN leaves the previously certified
        // invariant unchanged and must not make the exact PDR property query
        // inconclusive.
        return false;  // LCOV_EXCL_LINE
      }
      if (status == SATSolverWrapper::SolveStatus::Unsat) {
        buildInvariant(result);
        return true;
      }
      const size_t removed =
          removeViolatedCandidates(solver, nextHolds);
      if (removed == 0) {
        return false;  // LCOV_EXCL_LINE
      }
      result.inductiveRejected += removed;
    }
    return true;
  }

  PDRExactInitCache::Impl& cache_;
  const KInductionProblem& problem_;
  BoolExpr* initFormula_ = nullptr;
  const FrameClauses& baseInvariant_;
  std::vector<StateClause> candidates_;
  std::unordered_set<size_t> candidateSymbols_;
  std::unordered_set<size_t> currentConstraintSymbols_;
  std::vector<bool> active_;
};

uint64_t reusableInvariantCertificationTotalTickLimit() {
  if (activePdrQueryLimits != nullptr &&
      activePdrQueryLimits->invariantCertificationTotalTickLimit >= 0) {
    return static_cast<uint64_t>(
        activePdrQueryLimits->invariantCertificationTotalTickLimit);
  }
  return kDualRailInvariantCertificationTotalTickLimit;
}

void retainOnlyCertifiedReusableInvariant(
    PDRExactInitCache::Impl& cache) {
  cache.reusableInvariantCandidates = cache.reusableInvariant.clauses;
  cache.reusableInvariantCandidateLiteralCount = 0;
  for (const StateClause& clause : cache.reusableInvariantCandidates) {
    cache.reusableInvariantCandidateLiteralCount += clause.size();
  }
  cache.reusableInvariantCertifiedRevision =
      cache.reusableInvariantCandidateRevision;
}

void initializeReusableInvariantCertificationBudget(
    PDRExactInitCache::Impl& cache) {
  if (cache.reusableInvariantCertificationBudget.has_value() ||
      !cache.sourceProblem->usesDualRailStateEncoding ||
      SATSolverWrapper::assumptionSolverTypeFor(cache.solverType) !=
          KEPLER_FORMAL::Config::SolverType::CADICAL) {
    return;
  }
  const uint64_t tickLimit =
      reusableInvariantCertificationTotalTickLimit();
  cache.reusableInvariantCertificationBudget.emplace(
      std::numeric_limits<uint64_t>::max(),
      std::numeric_limits<uint64_t>::max(),
      tickLimit);
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: reusable invariant cumulative certification budget "
        "created ticks=",
        tickLimit);
  }
}

void disableReusableInvariantCertification(
    PDRExactInitCache::Impl& cache) {
  size_t discarded = 0;
  for (const StateClause& candidate :
       cache.reusableInvariantCandidates) {
    if (!hasExactClause(cache.reusableInvariant.clauses, candidate)) {
      ++discarded;
    }
  }
  retainOnlyCertifiedReusableInvariant(cache);
  cache.reusableInvariantCertificationDisabled = true;
  if (pdrStatsEnabled()) {
    const auto& budget = cache.reusableInvariantCertificationBudget;
    emitSecDiag(
        "SEC PDR stats: reusable invariant cumulative certification budget "
        "exhausted; disabling optional certification ticks=",
        budget.has_value() ? budget->ticksUsed() : 0,
        "/",
        budget.has_value() ? budget->tickLimit() : 0,
        " discarded_candidates=",
        discarded,
        " retained=",
        cache.reusableInvariant.clauses.size());
  }
}

void refreshReusableInvariant(
    PDRExactInitCache::Impl& cache,
    BoolExpr* initFormula) {
  if (cache.reusableInvariantCertificationDisabled) {
    return;
  }
  if (cache.reusableInvariantCandidateRevision ==
      cache.reusableInvariantCertifiedRevision) {
    return;
  }
  initializeReusableInvariantCertificationBudget(cache);

  // Invariant reuse is optional. Any inability to certify the new candidates
  // leaves the previous proven H intact and cannot alter this PDR run.
  try {
    ReusableInvariantFinder finder(cache, initFormula);
    const auto result = finder.run();
    if (!result.has_value()) {
      if (cache.reusableInvariantCertificationBudget.has_value() &&
          cache.reusableInvariantCertificationBudget->exhausted()) {
        disableReusableInvariantCertification(cache);
      }
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: reusable invariant certification unavailable");
      }
      return;
    }
    const size_t candidateCount =
        cache.reusableInvariantCandidates.size();
    cache.reusableInvariant = result->invariant;
    // The FMCAD'11 removal loop permanently drops rejected candidates. Keep
    // only certified H so later batches extend that invariant instead of
    // repeatedly solving every clause already rejected by an exact query.
    retainOnlyCertifiedReusableInvariant(cache);
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: reusable invariant clauses certified candidates=",
          candidateCount,
          " retained=",
          cache.reusableInvariant.clauses.size(),
          " initial_rejected=",
          result->initialRejected,
          " inductive_rejected=",
          result->inductiveRejected,
          " queries=",
          result->inductiveQueries);
    }
    if (cache.reusableInvariantCertificationBudget.has_value() &&
        cache.reusableInvariantCertificationBudget->exhausted()) {
      disableReusableInvariantCertification(cache);
    }
  } catch (...) {  // LCOV_EXCL_LINE
    if (cache.reusableInvariantCertificationBudget.has_value() &&  // LCOV_EXCL_LINE
        cache.reusableInvariantCertificationBudget->exhausted()) {  // LCOV_EXCL_LINE
      disableReusableInvariantCertification(cache);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC PDR stats: reusable invariant certification failed");  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
}

void addFrameConstraints(SATSolverWrapper& solver,
                         const FrameVariableStore& variables,
                         BoolExpr* initFormula,
                         BoolExpr* frameInvariant,
                         const std::vector<FrameClauses>& frames,
                         size_t level,
                         size_t frame) {
  if (level == 0) {
    // Figure 6 uses F[0] = I. Every level-zero query therefore receives the
    // complete exact startup formula; no cone-local assignment projection is
    // substituted for I.
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    solver.addClause({encoder.encode(initFormula)});
    addAllFrameClauses(solver, variables, frames[0], frame);
    return;
  }

  // For higher frames, materialize the currently learned blocking clauses and
  // LCOV_EXCL_START
  // any validated strengthening invariant the strategy handed to PDR.
  addAllFrameClauses(solver, variables, frames[level], frame);
  if (frameInvariant != nullptr) {
    // The optional strengthening is treated exactly like a frame fact, but it
    // is validated before we allow the engine to rely on it.
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));  // LCOV_EXCL_LINE
    solver.addClause({encoder.encode(frameInvariant)});  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
}

bool predecessorSourceFrameIsKnownSafe(size_t level) {
  // Predecessor queries are only issued from frames that were already checked
  // safe in an earlier PDR phase: blocking a bad cube at F[i+1] asks from F[i],
  // and propagation runs after the current frontier has been exhausted. F[0]
  // is the startup frontier and is handled separately by Init/reset facts.
  return level > 0;
}

void addSafeFramePropertyConstraint(SATSolverWrapper& solver,
                                    const FrameVariableStore& variables,
                                    const KInductionProblem& problem,
                                    size_t level,
                                    PdrFormulaSupportCache* supportCache,
                                    // LCOV_EXCL_START
                                    size_t frame) {
  if (!predecessorSourceFrameIsKnownSafe(level) || problem.property == nullptr) {
    return;
  }
  // LCOV_EXCL_STOP
  // The property is logically redundant for an exact safe PDR frame, but
  // keeping it explicit strengthens the predecessor SAT query. Encode only
  // its exact support so a local property does not reserve Tseitin storage from
  // an unrelated ASIC-sized predecessor surface.
  const std::set<size_t> uncachedSupport =
      supportCache == nullptr ? problem.property->getSupportVars()
                              : std::set<size_t>{};
  const std::set<size_t>& propertySupport =
      supportCache != nullptr ? supportCache->support(problem.property)
                              : uncachedSupport;
  FrameFormulaEncoder encoder(
      solver, variables.makeLeafLits(frame, propertySupport));  // LCOV_EXCL_LINE
  solver.addClause({encoder.encode(problem.property)});  // LCOV_EXCL_LINE
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: safe frame property support symbols=",
        propertySupport.size());
  }
}

bool predecessorFrameClauseApplies(
    const PredecessorAssumptionSolver& cachedSolver,
    const StateClause& clause);

void addGuardedFormulaConstraint(
    PredecessorAssumptionSolver& cachedSolver,
    BoolExpr* formula,
    PdrFormulaSupportCache* supportCache,
    int activationLiteral) {
  if (formula == nullptr) {
    return;
  }
  const std::set<size_t> uncachedSupport =
      supportCache == nullptr ? formula->getSupportVars() : std::set<size_t>{};
  const std::set<size_t>& formulaSupport =
      supportCache != nullptr ? supportCache->support(formula)
                              : uncachedSupport;
  FrameFormulaEncoder encoder(
      *cachedSolver.solver,
      cachedSolver.variables->makeLeafLits(0, formulaSupport));
  cachedSolver.solver->addClause(
      {-activationLiteral, encoder.encode(formula)});
}

int prepareGuardedPredecessorFrameContext(
    PredecessorAssumptionSolver& cachedSolver,
    size_t runId,
    const KInductionProblem& problem,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    PdrFormulaSupportCache* supportCache,
    bool scanCompleteFrame) {
  auto& context = cachedSolver.guardedContext;
  const bool newContext = context.runId != runId;
  if (newContext) {
    if (context.activationLiteral != 0) {
      // The previous batch has completed. Permanently retire its guarded facts;
      // learned clauses independent of that activation remain available.
      cachedSolver.solver->addClause({-context.activationLiteral});
    }
    context = GuardedPredecessorFrameContext{};
    context.runId = runId;
    context.activationLiteral = cachedSolver.solver->newVar() + 2;
    context.property = problem.property;
    context.frameInvariant = frameInvariant;
    ++cachedSolver.guardedContextCount;

    addGuardedFormulaConstraint(
        cachedSolver,
        frameInvariant,
        supportCache,
        context.activationLiteral);
    if (predecessorSourceFrameIsKnownSafe(level)) {
      addGuardedFormulaConstraint(
          cachedSolver,
          problem.property,
          supportCache,
          context.activationLiteral);
    }
    scanCompleteFrame = true;
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: shared predecessor context activated run=",
          runId,
          " level=",
          level,
          " retained_contexts=",
          cachedSolver.guardedContextCount);
    }
  } else if (context.property != problem.property ||
             context.frameInvariant != frameInvariant) {
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "PDR shared predecessor context changed within one run");  // LCOV_EXCL_LINE
  }

  const size_t frameFingerprint = frameClausesFingerprint(frames, level);
  if (!scanCompleteFrame &&
      context.emittedFrameFingerprint == frameFingerprint) {
    return context.activationLiteral;
  }

  const std::vector<StateClause>* clauses = &frames[level].addedClauseLog;
  size_t beginIndex = context.emittedFrameLogOffset;
  if (scanCompleteFrame || beginIndex > clauses->size()) {
    clauses = &frames[level].clauses;
    beginIndex = 0;
  }
  for (size_t clauseIndex = beginIndex; clauseIndex < clauses->size();
       ++clauseIndex) {
    const StateClause& clause = (*clauses)[clauseIndex];
    if (!predecessorFrameClauseApplies(cachedSolver, clause) ||
        !context.emittedFrameClauses.insert(clause).second) {
      continue;
    }
    addGuardedStateClause(
        *cachedSolver.solver,
        *cachedSolver.variables,
        clause,
        0,
        context.activationLiteral);
  }
  context.emittedFrameLogOffset = frames[level].addedClauseLog.size();
  context.emittedFrameFingerprint = frameFingerprint;
  return context.activationLiteral;
}

bool predecessorFrameClauseApplies(
    const PredecessorAssumptionSolver& cachedSolver,
    const StateClause& clause) {
  return clauseCoveredByVariables(*cachedSolver.variables, clause);
}

void rememberPredecessorFrameClauses(
    PredecessorAssumptionSolver& cachedSolver,
    const FrameClauses& frameClauses) {
  for (const auto& clause : frameClauses.clauses) {
    if (predecessorFrameClauseApplies(cachedSolver, clause)) {
      cachedSolver.emittedFrameClauses.insert(clause);
    }
  }
  cachedSolver.emittedFrameLogOffset = frameClauses.addedClauseLog.size();
}

size_t addNewPredecessorFrameClauses(
    PredecessorAssumptionSolver& cachedSolver,
    const FrameClauses& frameClauses,
    size_t frame,
    size_t frameFingerprint,
    bool scanCompleteFrame) {
  if (!scanCompleteFrame &&
      cachedSolver.emittedFrameFingerprint == frameFingerprint) {
    return 0;
  }

  const std::vector<StateClause>* clauses = &frameClauses.addedClauseLog;
  size_t beginIndex = cachedSolver.emittedFrameLogOffset;
  const char* source = "frame_log";
  if (scanCompleteFrame || beginIndex > clauses->size()) {
    // A widened symbol surface can make an older clause encodable for the
    // first time. Rescan only on that rare extension; ordinary frame updates
    // consume the append-only clause log below.
    clauses = &frameClauses.clauses;
    beginIndex = 0;
    source = "full_frame";
  }

  size_t addedClauses = 0;
  for (size_t clauseIndex = beginIndex; clauseIndex < clauses->size();
       ++clauseIndex) {
    const auto& clause = (*clauses)[clauseIndex];
    if (!predecessorFrameClauseApplies(cachedSolver, clause) ||
        !cachedSolver.emittedFrameClauses.insert(clause).second) {
      continue;
    }
    addStateClause(*cachedSolver.solver, *cachedSolver.variables, clause, frame);
    ++addedClauses;
  }
  cachedSolver.emittedFrameLogOffset = frameClauses.addedClauseLog.size();
  cachedSolver.emittedFrameFingerprint = frameFingerprint;
  if (shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: predecessor cached solver frame sync source=",
        source,
        " scanned=",
        clauses->size() - beginIndex,
        " added=",
        addedClauses,
        " level=",
        cachedSolver.key.level);
  }
  return addedClauses;
}

void PredecessorAssumptionSolver::extendSymbolSurface(
    const ComplementPartnerIndex& stateRelations,
    const std::vector<size_t>& solverSymbols) {
  std::vector<size_t> addedSymbols;
  std::set_difference(
      solverSymbols.begin(),
      solverSymbols.end(),
      key.solverSymbols.begin(),
      key.solverSymbols.end(),
      std::back_inserter(addedSymbols));
  if (addedSymbols.empty()) {
    return;
  }

  variables->addSymbols(*solver, addedSymbols);
  querySymbolSet.insert(addedSymbols.begin(), addedSymbols.end());
  for (const size_t symbol : addedSymbols) {
    const int literal = variables->getLiteral(symbol, 0);
    transitionLeafLits.emplace(symbol, literal);
    for (auto& [symbolMap, encoder] : transitionEncoderBySymbolMap) {
      (void)symbolMap;
      encoder->addLeafLiteral(symbol, literal);
    }
    if (badRootEncoder != nullptr) {
      badRootEncoder->addLeafLiteral(symbol, literal);
    }
  }

  // Existing transition roots and clauses remain valid. Extend each encoder's
  // leaf table so later roots reuse the exact same Tseitin DAG cache.
  if (!transitionEncoderBySymbolMap.empty() &&
      shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: predecessor transition encoder cache extended "
        "encoders=",
        transitionEncoderBySymbolMap.size(),
        " added_symbols=",
        addedSymbols.size());
  }
  key.solverSymbols = solverSymbols;

  // The symbol-surface builder closes every relation pair. Emit only pairs
  // made newly representable by this extension; old clauses remain permanent.
  stateRelations.addClausesForAddedSymbols(
      *solver, *variables, addedSymbols, 0);
  if (shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: predecessor relation surface extended added_symbols=",
        addedSymbols.size());
  }
}

PredecessorAssumptionSolver& getOrCreatePredecessorAssumptionSolver(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    const ComplementPartnerIndex& stateRelations,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const std::vector<size_t>& solverSymbols,
    PdrFormulaSupportCache* supportCache) {
  const bool useSharedFrameZeroSolver =
      level == 0 && cache.sharedFrameZeroPredecessorSolver != nullptr;
  const bool useSharedHigherFrameSolver =
      level > 0 && cache.sharedHigherFrameSolverPools != nullptr &&
      cache.sharedHigherFrameRunId != 0;
  SharedPredecessorAssumptionSolverPool::Selection sharedSelection;
  if (useSharedHigherFrameSolver) {
    if (cache.sharedHigherFrameFamilySymbols == nullptr) {  // LCOV_EXCL_LINE
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR shared predecessor property family is unavailable");  // LCOV_EXCL_LINE
    }
    sharedSelection =
        (*cache.sharedHigherFrameSolverPools)[level].select(
            cache.sharedHigherFrameRunId,
            *cache.sharedHigherFrameFamilySymbols,
            cache.usePathLocalHigherFrameSolverReuse);
    if (sharedSelection.selectedNewRun && pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: shared predecessor solver pool selected level=",
          level,
          " run=",
          cache.sharedHigherFrameRunId,
          " entry=",
          sharedSelection.entryIndex,
          " cache_hit=",
          sharedSelection.cacheHit ? 1 : 0,
          " evicted=",
          sharedSelection.evicted ? 1 : 0,
          " family_symbols=",
          cache.sharedHigherFrameFamilySymbols->size(),
          " initial_symbols=",
          solverSymbols.size(),
          " closest_entry=",
          sharedSelection.closestEntry,
          " closest_overlap=",
          sharedSelection.closestOverlap,
          " path_local=",
          cache.usePathLocalHigherFrameSolverReuse ? 1 : 0,
          " restarted=",
          sharedSelection.restarted ? 1 : 0,
          " retired_contexts=",
          sharedSelection.retiredContexts);
    }
  }
  auto& solver =
      useSharedFrameZeroSolver
          ? *cache.sharedFrameZeroPredecessorSolver
          : useSharedHigherFrameSolver
                ? *sharedSelection.solver
                : cache.solversByLevel[level];
  const KInductionProblem* keyProblem =
      useSharedFrameZeroSolver
          ? cache.sharedFrameZeroPredecessorProblem
          : useSharedHigherFrameSolver
                ? cache.sharedHigherFrameProblem
                : &problem;
  const void* keyTransitionModel =
      useSharedFrameZeroSolver
          ? cache.sharedFrameZeroTransitionModel
          : useSharedHigherFrameSolver
                ? cache.sharedHigherFrameTransitionModel
                : static_cast<const void*>(&transitionByState);
  const size_t currentFrameFingerprint =
      frameClausesFingerprint(frames, level);
  if (useSharedFrameZeroSolver && solver != nullptr &&
      solver->hasSameSharedFrameZeroContext(
          keyProblem,
          keyTransitionModel,
          initFormula,
          solverSymbols)) {
    // The caller borrowed the same complete F[0] vector. Skip copying and a
    // multi-million-element set difference; only frame-clause synchronization
    // remains, exactly as in the regular incremental-solver path below.
    const size_t addedClauses = addNewPredecessorFrameClauses(
        *solver,
        frames[level],
        0,
        currentFrameFingerprint,
        /*scanCompleteFrame=*/false);
    solver->key.frameFingerprint = currentFrameFingerprint;
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: shared exact F[0] predecessor solver reused "
          "without symbol surface comparison");
    }
    if (addedClauses != 0 && pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: predecessor cached solver frame clauses added=",
          addedClauses,
          " level=",
          level,
          " symbols=",
          solverSymbols.size());
    }
    return *solver;
  }
  std::vector<size_t> sharedSolverSymbols;
  const std::vector<size_t>* keySolverSymbols = &solverSymbols;
  if (useSharedHigherFrameSolver && solver != nullptr) {
    sharedSolverSymbols = detail::mergeSortedPdrSymbolVectors(
        solver->key.solverSymbols, solverSymbols);
    keySolverSymbols = &sharedSolverSymbols;
  }
  PredecessorAssumptionCacheKey key{
      keyProblem,
      keyTransitionModel,
      initFormula,
      level == 0 || useSharedHigherFrameSolver ? nullptr : frameInvariant,
      level,
      currentFrameFingerprint,
      *keySolverSymbols};
  if (solver != nullptr && solver->canExtendTo(key)) {
    const size_t previousSymbolCount = solver->key.solverSymbols.size();
    solver->extendSymbolSurface(stateRelations, key.solverSymbols);
    const bool surfaceWidened =
        solver->key.solverSymbols.size() != previousSymbolCount;
    if (surfaceWidened && shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: predecessor cached solver surface extended added=",
          solver->key.solverSymbols.size() - previousSymbolCount,
          " symbols=",
          solver->key.solverSymbols.size(),
          " level=",
          level);
    }
    size_t addedClauses = 0;
    if (useSharedHigherFrameSolver) {
      prepareGuardedPredecessorFrameContext(
          *solver,
          cache.sharedHigherFrameRunId,
          problem,
          frameInvariant,
          frames,
          level,
          supportCache,
          surfaceWidened);
    } else {
      // PDR frames strengthen monotonically. Reuse the expensive transition
      // solver, then stream only newly learned frame clauses into it.
      addedClauses = addNewPredecessorFrameClauses(
          *solver,
          frames[level],
          0,
          currentFrameFingerprint,
          surfaceWidened);
      solver->key.frameFingerprint = key.frameFingerprint;
    }
    if (useSharedFrameZeroSolver && shouldEmitFrequentPdrStats()) {
      emitSecDiag("SEC PDR stats: shared exact F[0] predecessor solver reused");
    }
    if (addedClauses != 0 && shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: predecessor cached solver frame clauses added=",
          addedClauses,
          " level=",
          level,
          " symbols=",
          solverSymbols.size());
    }
    return *solver;
  }

  auto next = std::make_unique<PredecessorAssumptionSolver>();
  next->key = std::move(key);
  next->sharedFrameZeroSolverSymbols =
      useSharedFrameZeroSolver
          ? cache.sharedFrameZeroPredecessorSymbols
          : nullptr;
  next->solver = std::make_unique<SATSolverWrapper>(
      SATSolverWrapper::assumptionSolverTypeFor(solverType));
  if (useSharedHigherFrameSolver) {
    next->solver->configureForSecPdrPersistentQuery(
        next->key.solverSymbols.size());
  } else {
    next->solver->configureForSecPdrQuery(next->key.solverSymbols.size());
  }
  next->variables =
      std::make_unique<FrameVariableStore>(
          *next->solver, next->key.solverSymbols, 1);
  next->querySymbolSet.insert(
      next->key.solverSymbols.begin(), next->key.solverSymbols.end());
  stateRelations.addClauses(
      *next->solver, *next->variables, next->key.solverSymbols, 1);
  if (useSharedHigherFrameSolver) {
    prepareGuardedPredecessorFrameContext(
        *next,
        cache.sharedHigherFrameRunId,
        problem,
        frameInvariant,
        frames,
        level,
        supportCache,
        /*scanCompleteFrame=*/true);
  } else {
    addFrameConstraints(*next->solver, *next->variables, initFormula,
                        frameInvariant, frames, level, 0);
    addSafeFramePropertyConstraint(
        *next->solver, *next->variables, problem, level, supportCache, 0);
  }
  addPostBootstrapResetInputConstraints(*next->solver, *next->variables,
                                        problem, 0);
  if (!useSharedHigherFrameSolver && level < frames.size()) {
    rememberPredecessorFrameClauses(*next, frames[level]);
    next->emittedFrameFingerprint = currentFrameFingerprint;
  }
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: predecessor cached solver created level=",
        level,
        " symbols=",
        next->key.solverSymbols.size(),
        " frame_clauses=",
        level < frames.size() ? frames[level].clauses.size() : 0,
        " shared_batches=",
        useSharedHigherFrameSolver ? 1 : 0);
  }
  solver = std::move(next);
  return *solver;
}

int64_t resourceLimitOrUnbounded(unsigned limit);

int encodeCachedFrameZeroBadRoot(
    PredecessorAssumptionSolver& cachedSolver,
    BoolExpr* badFormula) {
  const auto existing = cachedSolver.encodedBadRoots.find(badFormula);
  if (existing != cachedSolver.encodedBadRoots.end()) {
    return existing->second;
  }
  if (cachedSolver.badRootEncoder == nullptr) {
    cachedSolver.badRootEncoder = std::make_unique<FrameFormulaEncoder>(
        *cachedSolver.solver, cachedSolver.variables->makeLeafLits(0));
  }
  const int root = cachedSolver.badRootEncoder->encode(badFormula);
  cachedSolver.encodedBadRoots.emplace(badFormula, root);
  return root;
}

SATSolverWrapper::SolveStatus solveFrameZeroBadCubeWithSharedSolver(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    const ComplementPartnerIndex& stateRelations,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    BoolExpr* badFormula,
    const std::vector<size_t>& solverSymbols,
    unsigned badCubeConflictLimit,
    PdrFormulaSupportCache* supportCache,
    PredecessorAssumptionSolver** solvedCache) {
  auto& cachedSolver = getOrCreatePredecessorAssumptionSolver(
      cache, problem, solverType, transitionByState, stateRelations,
      initFormula, frameInvariant, frames,
      /*level=*/0, solverSymbols, supportCache);
  const int badRoot =
      encodeCachedFrameZeroBadRoot(cachedSolver, badFormula);
  *solvedCache = &cachedSolver;
  if (shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: shared exact F[0] solver used for bad cube");
  }
  return cachedSolver.solver->solveWithAssumptionsStatus(
      {badRoot}, resourceLimitOrUnbounded(badCubeConflictLimit),
      /*propagationLimit=*/-1);
}

int64_t resourceLimitOrUnbounded(unsigned limit) {
  return limit == 0 || limit == std::numeric_limits<unsigned>::max()
             ? -1
             : static_cast<int64_t>(limit);
}

PredecessorQueryResultKey makePredecessorQueryResultKey(
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    size_t level,
    size_t frameFingerprint,
    bool excludeTargetOnCurrentFrame,
    const StateCube& targetCube) {
  PredecessorQueryResultKey key;
  key.problem = &problem;
  key.transitionByState = &transitionByState;
  key.initFormula = initFormula;
  key.frameInvariant = frameInvariant;
  key.level = level;
  key.frameFingerprint = frameFingerprint;
  key.excludeTargetOnCurrentFrame = excludeTargetOnCurrentFrame;
  key.targetCube = targetCube;
  return key;
}

PredecessorQueryResultKey makeCachedPredecessorQueryResultKey(
    const PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    size_t level,
    size_t frameFingerprint,
    bool excludeTargetOnCurrentFrame,
    const StateCube& targetCube) {
  const bool usesSharedFrameZero =
      level == 0 && cache.sharedFrameZeroQueryResultStore != nullptr;
  const KInductionProblem& keyProblem =
      usesSharedFrameZero && cache.sharedFrameZeroQueryProblem != nullptr
          ? *cache.sharedFrameZeroQueryProblem
          : problem;
  const TransitionExprResolver& keyTransition =
      usesSharedFrameZero && cache.sharedFrameZeroQueryTransition != nullptr
          ? *cache.sharedFrameZeroQueryTransition
          : transitionByState;
  return makePredecessorQueryResultKey(
      keyProblem,
      keyTransition,
      initFormula,
      level == 0 ? nullptr : frameInvariant,
      level,
      frameFingerprint,
      excludeTargetOnCurrentFrame,
      targetCube);
}

const PredecessorQueryResultStore& predecessorQueryResultStoreFor(
    const PredecessorAssumptionCache& cache,
    size_t level) {
  return level == 0 && cache.sharedFrameZeroQueryResultStore != nullptr
             ? *cache.sharedFrameZeroQueryResultStore
             : cache.queryResultStore;
}

PredecessorQueryResultStore& predecessorQueryResultStoreFor(
    PredecessorAssumptionCache& cache,
    size_t level) {
  return level == 0 && cache.sharedFrameZeroQueryResultStore != nullptr
             ? *cache.sharedFrameZeroQueryResultStore
             : cache.queryResultStore;
}

std::optional<PredecessorQueryResultEntry> cachedPredecessorQueryResult(
    PredecessorAssumptionCache& cache,
    const PredecessorQueryResultKey& exactKey,
    const PredecessorQueryResultKey& stableUnsatKey) {
  auto& store = predecessorQueryResultStoreFor(cache, exactKey.level);
  if (const auto* exact = store.queryResults.find(exactKey);
      exact != nullptr) {
    return *exact;
  }
  if (store.unsatQueries.find(stableUnsatKey) != store.unsatQueries.end()) {
    return PredecessorQueryResultEntry{}; // LCOV_EXCL_LINE
  }
  return std::nullopt;
}

PredecessorUnsatCoreCacheKey makePredecessorUnsatCoreCacheKey(
    const PredecessorQueryResultKey& key) {
  PredecessorUnsatCoreCacheKey coreKey;
  coreKey.problem = key.problem;
  coreKey.transitionByState = key.transitionByState;
  coreKey.initFormula = key.initFormula;
  coreKey.frameInvariant = key.frameInvariant;
  coreKey.level = key.level;
  coreKey.excludeTargetOnCurrentFrame = key.excludeTargetOnCurrentFrame;
  return coreKey;
}

bool predecessorUnsatCoreCacheable(
    const PredecessorQueryResultKey& stableUnsatKey) {
  // Failed target-assumption cores are globally reusable only for the base
  // predecessor context. Temporary relative-induction assumptions can be part
  // of the UNSAT proof, so keep those answers in the exact target cache only.
  return detail::shouldSharePredecessorUnsatCore(
      stableUnsatKey.frameFingerprint,
      stableUnsatKey.excludeTargetOnCurrentFrame);
}

void rememberPredecessorUnsatCore(
    PredecessorAssumptionCache& cache,
    const PredecessorQueryResultKey& stableUnsatKey,
    StateCube core) {
  if (!predecessorUnsatCoreCacheable(stableUnsatKey)) {
    return;
  }
  normalizeCube(core);
  if (core.empty()) {
    return; // LCOV_EXCL_LINE
  }

  auto& store =
      predecessorQueryResultStoreFor(cache, stableUnsatKey.level);
  auto& cores =
      store.unsatCoresByContext[
          makePredecessorUnsatCoreCacheKey(stableUnsatKey)];
  for (const auto& existing : cores) {
    if (cubeContainsCube(core, existing)) {
      return;
    }
  }

  std::vector<StateCube> retained;
  retained.reserve(cores.size() + 1);
  for (auto& existing : cores) {
    if (!cubeContainsCube(existing, core)) {
      retained.push_back(std::move(existing));
    }
  }
  retained.push_back(std::move(core));
  sortStateCubesDeterministically(retained);
  if (retained.size() > kMaxPredecessorUnsatCoresPerContext) {
    retained.pop_back(); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  cores = std::move(retained);
}

std::optional<StateCube> cachedPredecessorUnsatCoreForTarget(
    const PredecessorAssumptionCache& cache,
    const PredecessorQueryResultKey& stableUnsatKey,
    const StateCube& targetCube) {
  if (!predecessorUnsatCoreCacheable(stableUnsatKey)) {
    return std::nullopt;
  }
  const auto& store =
      predecessorQueryResultStoreFor(cache, stableUnsatKey.level);
  const auto coreIt =
      store.unsatCoresByContext.find(
          makePredecessorUnsatCoreCacheKey(stableUnsatKey));
  if (coreIt == store.unsatCoresByContext.end()) {
    return std::nullopt;
  }
  for (const auto& core : coreIt->second) {
    if (cubeContainsCube(targetCube, core)) {
      return core;
    }
  }
  return std::nullopt;
}

void trimPredecessorQueryResultCache(PredecessorAssumptionCache& cache,
                                     size_t level) {
  auto& store = predecessorQueryResultStoreFor(cache, level);
  if (!detail::shouldResetPdrStableUnsatCache(
          store.unsatQueries.size(),
          kMaxPredecessorQueryResultCacheEntries)) {
    return;
  }
  // Stable UNSAT entries have their own bound. Filling the separate exact LRU
  // must not discard these still-valid answers for strengthened PDR frames.
  if (pdrStatsEnabled()) { // LCOV_EXCL_LINE
    emitSecDiag( // LCOV_EXCL_LINE
        "SEC PDR stats: predecessor stable UNSAT cache reset entries=",
        store.unsatQueries.size()); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  store.unsatQueries.clear(); // LCOV_EXCL_LINE
  store.unsatCoresByContext.clear(); // LCOV_EXCL_LINE
}

void rememberPredecessorQueryResult(
    PredecessorAssumptionCache& cache,
    const PredecessorQueryResultKey& exactKey,
    const PredecessorQueryResultKey& stableUnsatKey,
    const std::optional<StateCube>& predecessor,
    const StateCube* unsatCore = nullptr,
    bool predecessorExistsWithoutModel = false) {
  trimPredecessorQueryResultCache(cache, exactKey.level);
  auto& store = predecessorQueryResultStoreFor(cache, exactKey.level);
  PredecessorQueryResultEntry entry;
  if (predecessor.has_value() || predecessorExistsWithoutModel) {
    entry.hasPredecessor = true;
    entry.hasPredecessorModel = predecessor.has_value();
    if (predecessor.has_value()) {
      entry.predecessor = *predecessor;
    }
  } else {
    if (unsatCore != nullptr && !unsatCore->empty()) {
      entry.hasUnsatCore = true;
      entry.unsatCore = *unsatCore;
      normalizeCube(entry.unsatCore);
    }
    store.unsatQueries.insert(stableUnsatKey);
  }
  const auto inserted =
      store.queryResults.insert(exactKey, std::move(entry), /*weight=*/1);
  if (inserted.evictedEntries != 0 && shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: predecessor exact result cache evicted entries=",
        inserted.evictedEntries,
        " retained_entries=",
        store.queryResults.size());
  }
  if (unsatCore != nullptr && !unsatCore->empty()) {
    rememberPredecessorUnsatCore(cache, stableUnsatKey, *unsatCore);
  }
}

std::optional<StateCube> cachedPredecessorUnsatCoreForCube(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t sourceLevel,
    const StateCube& targetCube,
    bool excludeTargetOnCurrentFrame) {
  if (sourceLevel >= frames.size()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const auto exactKey = makeCachedPredecessorQueryResultKey(
      cache,
      problem,
      transitionByState,
      initFormula,
      frameInvariant,
      sourceLevel,
      frameClausesFingerprint(frames, sourceLevel),
      excludeTargetOnCurrentFrame,
      targetCube);
  auto& store = predecessorQueryResultStoreFor(cache, sourceLevel);
  const auto* result = store.queryResults.find(exactKey);
  if (result != nullptr && result->hasUnsatCore &&
      !result->unsatCore.empty()) {
    return result->unsatCore;
  }
  const auto stableUnsatKey = makeCachedPredecessorQueryResultKey( // LCOV_EXCL_LINE
      cache, // LCOV_EXCL_LINE
      problem, // LCOV_EXCL_LINE
      transitionByState, // LCOV_EXCL_LINE
      initFormula, // LCOV_EXCL_LINE
      frameInvariant, // LCOV_EXCL_LINE
      sourceLevel, // LCOV_EXCL_LINE
      /*frameFingerprint=*/0,
      excludeTargetOnCurrentFrame, // LCOV_EXCL_LINE
      targetCube); // LCOV_EXCL_LINE
  return cachedPredecessorUnsatCoreForTarget( // LCOV_EXCL_LINE
      cache, stableUnsatKey, targetCube); // LCOV_EXCL_LINE
}

bool badCubeFrameClauseApplies(const BadCubeAssumptionSolver& cachedSolver,
                               const StateClause& clause) {
  return clauseCoveredByVariables(*cachedSolver.variables, clause);
}

void rememberBadCubeFrameClauses(BadCubeAssumptionSolver& cachedSolver,
                                 const FrameClauses& frameClauses) {
  for (const auto& clause : frameClauses.clauses) {
    if (badCubeFrameClauseApplies(cachedSolver, clause)) {
      cachedSolver.emittedFrameClauses.insert(clause);
    }
  }
}

void addNewBadCubeFrameClauses(BadCubeAssumptionSolver& cachedSolver,
                               const std::vector<StateClause>& frameClauses,
                               size_t beginIndex,
                               size_t frame,
                               const char* source) {
  size_t addedClauses = 0;
  for (size_t clauseIndex = beginIndex; clauseIndex < frameClauses.size();
       ++clauseIndex) {
    const auto& clause = frameClauses[clauseIndex];
    if (!badCubeFrameClauseApplies(cachedSolver, clause) ||
        !cachedSolver.emittedFrameClauses.insert(clause).second) {
      continue;
    }
    // Frame vectors are compacted by subsumption, so a stronger learned clause
    // can replace a weaker one without increasing the vector size. Track by
    // clause identity instead of append index to keep cached bad-cube solvers
    // synchronized with the logical frame.
    addStateClause(*cachedSolver.solver, *cachedSolver.variables, clause, frame);
    ++addedClauses;
  }
  if (addedClauses != 0 && pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: bad cube cached frame clauses added=",
        addedClauses,
        " frame=",
        frame,
        " source=",
        source,
        " scanned=",
        frameClauses.size() - beginIndex);
  }
}

void syncBadCubeFrameClauses(BadCubeAssumptionSolver& cachedSolver,
                             const FrameClauses& frameClauses,
                             size_t frame,
                             size_t frameFingerprint) {
  if (cachedSolver.emittedFrameFingerprint == frameFingerprint) {
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: bad cube cached frame clauses unchanged frame=",
          frame,
          " fingerprint=",
          frameFingerprint);
    }
    return;
  }
  if (cachedSolver.emittedFrameLogOffset <=
      frameClauses.addedClauseLog.size()) {
    addNewBadCubeFrameClauses(
        cachedSolver,
        frameClauses.addedClauseLog,
        cachedSolver.emittedFrameLogOffset,
        frame,
        "frame_log");
    cachedSolver.emittedFrameLogOffset = frameClauses.addedClauseLog.size();
  } else {
    addNewBadCubeFrameClauses( // LCOV_EXCL_LINE
        cachedSolver, // LCOV_EXCL_LINE
        frameClauses.clauses, // LCOV_EXCL_LINE
        0,
        frame, // LCOV_EXCL_LINE
        "full_frame");
    cachedSolver.emittedFrameLogOffset = frameClauses.addedClauseLog.size(); // LCOV_EXCL_LINE
  }
  cachedSolver.emittedFrameFingerprint = frameFingerprint;
}

std::optional<SATSolverWrapper::SolveStatus>
solvePredecessorCubeWithCachedAssumptions(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    const ComplementPartnerIndex& stateRelations,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const PredecessorTargetSurface& targetSurface,
    const std::vector<size_t>& solverSymbols,
    bool excludeTargetOnCurrentFrame,
    PredecessorQueryPurpose purpose,
    unsigned predecessorConflictLimit,
    unsigned predecessorDecisionLimit,
    PdrFormulaSupportCache* supportCache,
    PredecessorAssumptionSolver** solvedCache = nullptr,
    StateCube* solvedUnsatCore = nullptr) {
  auto& cachedSolver = getOrCreatePredecessorAssumptionSolver(
      cache, problem, solverType, transitionByState, stateRelations,
      initFormula, frameInvariant, frames, level, solverSymbols, supportCache);
  const auto& preparedTarget = cachedSolver.prepareTargetAssumptions(
      transitionByState,
      0,
      targetSurface.exclusionClause,
      targetSurface.transitionGroups);
  const std::vector<int>* queryAssumptions = &preparedTarget.assumptions;
  const bool useGuardedBatchContext =
      level > 0 && cache.sharedHigherFrameSolverPools != nullptr &&
      cachedSolver.guardedContext.runId == cache.sharedHigherFrameRunId;
  if (useGuardedBatchContext || excludeTargetOnCurrentFrame) {
    cachedSolver.targetAssumptions = preparedTarget.assumptions;
    if (useGuardedBatchContext) {
      cachedSolver.targetAssumptions.push_back(
          cachedSolver.guardedContext.activationLiteral);
    }
    if (excludeTargetOnCurrentFrame) {
      cachedSolver.targetAssumptions.push_back(cachedSolver.q2SelectorFor(
          targetSurface.exclusionClause,
          0,
          predecessorQueryNeedsModel(purpose)));
    }
    queryAssumptions = &cachedSolver.targetAssumptions;
  }
  if (queryAssumptions->empty()) {
    return std::nullopt; // LCOV_EXCL_LINE
  }
  if (solvedCache != nullptr) {
    *solvedCache = &cachedSolver;
  }
  // Section V of the paper keeps one incremental SAT instance and changes only
  // its assumptions between predecessor queries. A resource-limit hit is
  // UNKNOWN and must make this output inconclusive; it is never a proof result.
  const int64_t cachedPropagationLimit =
      resourceLimitOrUnbounded(predecessorDecisionLimit);
  const auto status = cachedSolver.solver->solveWithAssumptionsStatus(
      *queryAssumptions,
      resourceLimitOrUnbounded(predecessorConflictLimit),
      cachedPropagationLimit);
  if (status == SATSolverWrapper::SolveStatus::Unsat &&
      solvedUnsatCore != nullptr) {
    // Only target-cube assumptions are mapped back. Temporary selector
    // assumptions may participate in the SAT proof, but they are not state
    // literals that can form a learned PDR blocker.
    *solvedUnsatCore = cachedPredecessorUnsatCoreFromTargetContext(
        *cachedSolver.solver, preparedTarget);
  }
  return status;
}

BadCubeAssumptionSolver& getOrCreateBadCubeAssumptionSolver(
    BadCubeAssumptionCache& cache,
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const ComplementPartnerIndex& stateRelations,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const std::vector<size_t>& solverSymbols) {
  BadCubeAssumptionCacheKey key{
      &problem,
      initFormula,
      level == 0 ? nullptr : frameInvariant,
      level,
      solverSymbols};
  const size_t currentFrameFingerprint =
      frameClausesFingerprint(frames, level);
  if (cache.solver != nullptr && cache.solver->key == key) {
    syncBadCubeFrameClauses(
        *cache.solver,
        frames[level],
        0,
        currentFrameFingerprint);
    return *cache.solver;
  }

  auto next = std::make_unique<BadCubeAssumptionSolver>();
  next->key = std::move(key);
  next->solver = std::make_unique<SATSolverWrapper>(
      SATSolverWrapper::assumptionSolverTypeFor(solverType));
  next->solver->configureForSecPdrQuery(solverSymbols.size());
  next->variables =
      std::make_unique<FrameVariableStore>(*next->solver, solverSymbols, 1);
  next->querySymbolSet.insert(solverSymbols.begin(), solverSymbols.end());
  stateRelations.addClauses(*next->solver, *next->variables, solverSymbols, 1);
  addFrameConstraints(*next->solver, *next->variables, initFormula,
                      frameInvariant, frames, level, 0);
  addPostBootstrapResetInputConstraints(*next->solver, *next->variables,
                                        problem, 0);
  next->encoder = std::make_unique<FrameFormulaEncoder>(
      *next->solver, next->variables->makeLeafLits(0));
  if (level < frames.size()) {
    rememberBadCubeFrameClauses(*next, frames[level]);
    next->emittedFrameFingerprint = currentFrameFingerprint;
    next->emittedFrameLogOffset = frames[level].addedClauseLog.size();
  }
  cache.solver = std::move(next);
  return *cache.solver;
}

int encodeCachedBadRoot(BadCubeAssumptionSolver& cachedSolver,
                        BoolExpr* badFormula) {
  const auto existing = cachedSolver.encodedBadRoots.find(badFormula);
  if (existing != cachedSolver.encodedBadRoots.end()) {
    return existing->second;
  }
  const int root = cachedSolver.encoder->encode(badFormula);
  cachedSolver.encodedBadRoots.emplace(badFormula, root);
  return root;
}

SATSolverWrapper::SolveStatus solveBadCubeWithCachedAssumption(
    BadCubeAssumptionCache& cache,
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const ComplementPartnerIndex& stateRelations,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    BoolExpr* badFormula,
    const std::vector<size_t>& solverSymbols,
    unsigned badCubeConflictLimit,
    BadCubeAssumptionSolver** solvedCache) {
  auto& cachedSolver = getOrCreateBadCubeAssumptionSolver(
      cache, problem, solverType, stateRelations, initFormula, frameInvariant,
      frames, level, solverSymbols);
  const int badRoot = encodeCachedBadRoot(cachedSolver, badFormula);
  *solvedCache = &cachedSolver;
  // The cached solver keeps learned clauses across monotonic frame updates.
  // Keep the conflict cap for workflow safety, but do not cap decisions here:
  // on wide dual-rail datapaths CaDiCaL otherwise returns UNKNOWN before those
  // learned clauses can pay back the reused frame context.
  return cachedSolver.solver->solveWithAssumptionsStatus(
      {badRoot},
      resourceLimitOrUnbounded(badCubeConflictLimit),
      /*propagationLimit=*/-1);
} // LCOV_EXCL_LINE

StateCube extractStateCube(const SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const std::vector<size_t>& stateSymbols,
                           size_t frame) {
  StateCube cube;
  cube.reserve(stateSymbols.size());
  for (const auto symbol : stateSymbols) {
    cube.push_back({symbol, solver.getLiteralValue(variables.getLiteral(symbol, frame))});
  }
  normalizeCube(cube);
  return cube;
}

class PdrTernaryModelReducer {
 public:
  PdrTernaryModelReducer(
      const SATSolverWrapper& solver,
      const FrameVariableStore& variables,
      const std::vector<PdrTernarySimulationRoot>& roots,
      const StateCube& modelCube,
      PdrFormulaSupportCache* supportCache) {
    supportCache_ = supportCache != nullptr ? supportCache : &localSupportCache_;
    memoGeneration_ = supportCache_->nextTernaryEvaluationGeneration();
    roots_.reserve(roots.size());
    assignments_.reserve(modelCube.size());
    for (const auto& spec : roots) {
      Root root;
      root.formula = supportCache_->ternaryNodeIndex(spec.formula);
      root.symbolMap = spec.symbolMap;
      root.expectedValue = spec.expectedValue;
      root.localSupport = &supportCache_->support(spec.formula);
      root.support = &supportCache_->mappedTernarySupport(
          spec.formula, spec.symbolMap);
      roots_.push_back(std::move(root));
      const size_t rootIndex = roots_.size() - 1;
      for (const size_t symbol : *roots_.back().support) {
        rootIndicesBySymbol_[symbol].push_back(rootIndex);
      }
      Root& insertedRoot = roots_.back();
      insertedRoot.memo =
          &supportCache_->ternaryEvaluationMemo(insertedRoot.symbolMap);
      auto& dependencies =
          memoDependenciesBySymbolMap_[insertedRoot.symbolMap];
      dependencies.memo = insertedRoot.memo;
      for (const size_t localSymbol : *insertedRoot.localSupport) {
        const auto symbol = mappedSymbol(insertedRoot, localSymbol);
        if (symbol.has_value() && *symbol >= 2) {
          dependencies.localSymbolsByMappedSymbol[*symbol].push_back(
              localSymbol);
        }
      }
    }
    for (auto& [symbolMap, dependencies] :
         memoDependenciesBySymbolMap_) {
      (void)symbolMap;
      for (auto& [mappedSymbol, localSymbols] :
           dependencies.localSymbolsByMappedSymbol) {
        (void)mappedSymbol;
        std::sort(localSymbols.begin(), localSymbols.end());
        localSymbols.erase(
            std::unique(localSymbols.begin(), localSymbols.end()),
            localSymbols.end());
      }
    }

    for (const auto& literal : modelCube) {
      assignments_.emplace(
          literal.symbol, TernaryAssignment{literal.value, false});
    }
    for (const auto& root : roots_) {
      for (const size_t symbol : *root.support) {
        if (assignments_.find(symbol) == assignments_.end() &&
            variables.hasSymbol(symbol)) {
          assignments_.emplace(
              symbol,
              TernaryAssignment{
                  solver.getLiteralValue(variables.getLiteral(symbol, 0)),
                  false});
        }
      }
    }
  }

  StateCube reduce(const StateCube& modelCube) {
    if (roots_.empty() || !rootsHaveExpectedValues(std::nullopt)) {
      return modelCube;
    }

    StateCube reduced;
    reduced.reserve(modelCube.size());
    for (const auto& literal : modelCube) {
      if (!anyRootUses(literal.symbol)) {
        continue;
      }

      // FMCAD'11 Section III-B keeps successfully removed literals at X while
      // probing later literals. Store that cumulative X state beside the
      // concrete value so variable evaluation needs only one hash lookup.
      setAssignmentUnknown(literal.symbol, true);
      if (rootsHaveExpectedValues(literal.symbol)) {
        continue;
      }
      setAssignmentUnknown(literal.symbol, false);
      reduced.push_back(literal);
    }
    if (reusedEvaluationMemoEntries_ != 0 &&
        shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: ternary evaluation memo storage reused entries=",
          reusedEvaluationMemoEntries_,
          " root_index_symbols=",
          rootIndicesBySymbol_.size());
    }
    if (recomputedEvaluationParents_ != 0 &&
        shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: ternary incremental propagation changed_values=",
          propagatedEvaluationValueChanges_,
          " recomputed_parents=",
          recomputedEvaluationParents_,
          " stable_parents=",
          stableEvaluationParents_);
    }
    return reduced;
  }

 private:
  struct TernaryAssignment {
    bool value = false;
    bool unknown = false;
  };

  using EvaluationMemo = PdrFormulaSupportCache::TernaryEvaluationMemo;
  using EvaluationMemoEntry =
      PdrFormulaSupportCache::TernaryEvaluationMemoEntry;

  struct Root {
    size_t formula = PdrFormulaSupportCache::kInvalidTernaryNode;
    const std::unordered_map<size_t, size_t>* symbolMap = nullptr;
    bool expectedValue = false;
    const std::set<size_t>* localSupport = nullptr;
    const std::vector<size_t>* support = nullptr;
    EvaluationMemo* memo = nullptr;
  };

  struct MemoDependencies {
    EvaluationMemo* memo = nullptr;
    std::unordered_map<size_t, std::vector<size_t>>
        localSymbolsByMappedSymbol;
  };

  void setAssignmentUnknown(size_t symbol, bool unknown) {
    auto& assignment = assignments_.at(symbol);
    if (assignment.unknown == unknown) {
      return;
    }
    assignment.unknown = unknown;
    if (std::find(
            changedAssignmentSymbols_.begin(),
            changedAssignmentSymbols_.end(),
            symbol) == changedAssignmentSymbols_.end()) {
      changedAssignmentSymbols_.push_back(symbol);
    }
  }

  std::optional<bool> assignmentValue(size_t symbol) const {
    if (symbol < 2) {
      return symbol == 1;
    }
    const auto assignment = assignments_.find(symbol);
    if (assignment == assignments_.end() || assignment->second.unknown) {
      return std::nullopt;
    }
    return assignment->second.value;
  }

  std::optional<bool> evaluatedValue(
      size_t nodeIndex,
      const EvaluationMemo& memo) const {
    if (nodeIndex == PdrFormulaSupportCache::kInvalidTernaryNode) {
      return std::nullopt;
    }
    const EvaluationMemoEntry& entry = memo[nodeIndex];
    if (entry.generation != memoGeneration_) {
      return std::nullopt;
    }
    return entry.value;
  }

  std::optional<bool> evaluateOperator(
      Op op,
      std::optional<bool> lhs,
      std::optional<bool> rhs) const {
    switch (op) {
      case Op::NOT:
        return lhs.has_value() ? std::optional<bool>(!*lhs) : std::nullopt;
      case Op::AND:
        if ((lhs.has_value() && !*lhs) ||
            (rhs.has_value() && !*rhs)) {
          return false;
        }
        return lhs.has_value() && rhs.has_value()
                   ? std::optional<bool>(*lhs && *rhs)
                   : std::nullopt;
      case Op::OR:
        if ((lhs.has_value() && *lhs) ||
            (rhs.has_value() && *rhs)) {
          return true;
        }
        return lhs.has_value() && rhs.has_value()
                   ? std::optional<bool>(*lhs || *rhs)
                   : std::nullopt;
      case Op::XOR:
        return lhs.has_value() && rhs.has_value()
                   ? std::optional<bool>(*lhs != *rhs)
                   : std::nullopt;
      case Op::VAR:
      case Op::NONE:
      default:
        return std::nullopt;
    }
  }

  std::optional<bool> recomputeNodeValue(
      size_t nodeIndex,
      const EvaluationMemo& memo) const {
    const auto& node = supportCache_->ternaryNode(nodeIndex);
    return evaluateOperator(
        node.op,
        evaluatedValue(node.left, memo),
        evaluatedValue(node.right, memo));
  }

  void enqueueComputedParents(
      MemoDependencies& dependencies,
      size_t nodeIndex) {
    for (const size_t parent :
         supportCache_->ternaryParentNodes(nodeIndex)) {
      EvaluationMemoEntry& parentEntry = (*dependencies.memo)[parent];
      if (parentEntry.generation != memoGeneration_ ||
          parentEntry.queuedGeneration == memoGeneration_) {
        continue;
      }
      parentEntry.queuedGeneration = memoGeneration_;
      propagationWorklist_.push(parent);
    }
  }

  void propagateChangedAssignments() {
    if (changedAssignmentSymbols_.empty()) {
      return;
    }

    for (auto& [symbolMap, dependencies] :
         memoDependenciesBySymbolMap_) {
      (void)symbolMap;
      for (const size_t changedSymbol : changedAssignmentSymbols_) {
        const auto localSymbols =
            dependencies.localSymbolsByMappedSymbol.find(changedSymbol);
        if (localSymbols ==
            dependencies.localSymbolsByMappedSymbol.end()) {
          continue;
        }
        for (const size_t localSymbol : localSymbols->second) {
          const auto& variableNodes =
              supportCache_->ternaryVariableNodes(localSymbol);
          for (const size_t variableNode : variableNodes) {
            EvaluationMemoEntry& entry =
                (*dependencies.memo)[variableNode];
            if (entry.generation != memoGeneration_) {
              continue;
            }
            const auto value = assignmentValue(changedSymbol);
            if (entry.value == value) {
              continue;
            }
            entry.value = value;
            ++propagatedEvaluationValueChanges_;
            enqueueComputedParents(dependencies, variableNode);
          }
        }
      }

      // Initial root evaluation computes both children of every operator. The
      // worklist can therefore update parents directly from their current
      // child values. A parent is requeued if another child changes later, and
      // propagation stops as soon as the exact ternary value remains stable.
      while (!propagationWorklist_.empty()) {
        const size_t nodeIndex = propagationWorklist_.front();
        propagationWorklist_.pop();
        EvaluationMemoEntry& entry = (*dependencies.memo)[nodeIndex];
        entry.queuedGeneration = 0;
        const auto value = recomputeNodeValue(nodeIndex, *dependencies.memo);
        ++recomputedEvaluationParents_;
        if (entry.value == value) {
          ++stableEvaluationParents_;
          continue;
        }
        entry.value = value;
        ++propagatedEvaluationValueChanges_;
        enqueueComputedParents(dependencies, nodeIndex);
      }
    }
    changedAssignmentSymbols_.clear();
  }

  std::optional<size_t> mappedSymbol(const Root& root,
                                     size_t symbol) const {
    if (symbol < 2 || root.symbolMap == nullptr) {
      return symbol;
    }
    const auto mapped = root.symbolMap->find(symbol);
    if (mapped == root.symbolMap->end()) {
      return std::nullopt;
    }
    return mapped->second;
  }

  std::optional<bool> evaluate(
      size_t nodeIndex,
      const Root& root,
      EvaluationMemo& memo) {
    if (nodeIndex == PdrFormulaSupportCache::kInvalidTernaryNode) {
      return std::nullopt;
    }
    EvaluationMemoEntry& entry = memo[nodeIndex];
    if (entry.generation == memoGeneration_) {
      ++reusedEvaluationMemoEntries_;
      return entry.value;
    }

    const auto& node = supportCache_->ternaryNode(nodeIndex);
    std::optional<bool> result;
    switch (node.op) {
      case Op::VAR: {
        const auto symbol = mappedSymbol(root, node.symbol);
        if (!symbol.has_value()) {
          break;
        }
        result = assignmentValue(*symbol);
        break;
      }
      case Op::NOT:
      case Op::AND:
      case Op::OR:
      case Op::XOR: {
        // Evaluate the complete immutable DAG once. Later literal trials use
        // parent links to update only values changed by the new X assignment.
        const auto lhs = evaluate(node.left, root, memo);
        const auto rhs = evaluate(node.right, root, memo);
        result = evaluateOperator(node.op, lhs, rhs);
        break;
      }
      case Op::NONE:
      default:
        break;
    }
    entry.generation = memoGeneration_;
    entry.value = result;
    return result;
  }

  bool rootHasExpectedValue(const Root& root) {
    const auto value = evaluate(root.formula, root, *root.memo);
    return value.has_value() && *value == root.expectedValue;
  }

  bool rootsHaveExpectedValues(std::optional<size_t> changedSymbol) {
    propagateChangedAssignments();
    if (!changedSymbol.has_value()) {
      for (const auto& root : roots_) {
        if (!rootHasExpectedValue(root)) {
          return false;
        }
      }
      return true;
    }

    // Indexing roots by support changes only host-side lookup work. It keeps
    // the original root order and evaluates exactly the roots whose value can
    // change when this state literal is tentatively replaced by X.
    const auto rootsIt = rootIndicesBySymbol_.find(*changedSymbol);
    if (rootsIt == rootIndicesBySymbol_.end()) {
      return true;
    }
    for (const size_t rootIndex : rootsIt->second) {
      const auto& root = roots_[rootIndex];
      if (!rootHasExpectedValue(root)) {
        return false;
      }
    }
    return true;
  }

  bool anyRootUses(size_t symbol) const {
    // The constructor already records every root using each symbol. Reuse that
    // exact index instead of rescanning every root for each model literal.
    return rootIndicesBySymbol_.contains(symbol);
  }

  std::vector<Root> roots_;
  PdrFormulaSupportCache localSupportCache_;
  PdrFormulaSupportCache* supportCache_ = nullptr;
  std::unordered_map<size_t, TernaryAssignment> assignments_;
  std::unordered_map<size_t, std::vector<size_t>> rootIndicesBySymbol_;
  std::unordered_map<const std::unordered_map<size_t, size_t>*,
                     MemoDependencies>
      memoDependenciesBySymbolMap_;
  std::vector<size_t> changedAssignmentSymbols_;
  std::queue<size_t> propagationWorklist_;
  size_t memoGeneration_ = 0;
  size_t reusedEvaluationMemoEntries_ = 0;
  size_t propagatedEvaluationValueChanges_ = 0;
  size_t recomputedEvaluationParents_ = 0;
  size_t stableEvaluationParents_ = 0;
};

StateCube reduceSolvedCubeByTernarySimulation(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<PdrTernarySimulationRoot>& roots,
    const StateCube& modelCube,
    PdrFormulaSupportCache* supportCache) {
  // Section III-B of the FMCAD'11 PDR paper removes a state literal only when
  // replacing it by X leaves every target value concrete and unchanged.
  const auto setupStart = std::chrono::steady_clock::now();
  PdrTernaryModelReducer reducer(
      solver, variables, roots, modelCube, supportCache);
  const auto reductionStart = std::chrono::steady_clock::now();
  StateCube reduced = reducer.reduce(modelCube);
  if (shouldEmitFrequentPdrStats()) {
    const auto reductionEnd = std::chrono::steady_clock::now();
    emitSecDiag(
        "SEC PDR stats: ternary reducer timing setup_us=",
        std::chrono::duration_cast<std::chrono::microseconds>(
            reductionStart - setupStart).count(),
        " reduction_us=",
        std::chrono::duration_cast<std::chrono::microseconds>(
            reductionEnd - reductionStart).count(),
        " roots=",
        roots.size(),
        " model_cube=",
        modelCube.size());
  }
  return reduced;
}

StateCube extractSolvedPredecessorCube(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<size_t>& predecessorSymbols,
    const std::vector<PdrTernarySimulationRoot>& roots,
    PdrFormulaSupportCache* supportCache) {
  const StateCube modelCube =
      extractStateCube(solver, variables, predecessorSymbols, 0);
  return reduceSolvedCubeByTernarySimulation(
      solver, variables, roots, modelCube, supportCache);
}

StateCube extractSolvedBadCubeForFormula(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<size_t>& badStateSupport,
    BoolExpr* badFormula,
    size_t level,
    PdrFormulaSupportCache* supportCache) {
  if (isSecDiagEnabled()) {
    emitSecDiag(  // LCOV_EXCL_LINE
        "SEC diag: PDR bad cube uses exact ternary support: ",
        badStateSupport.size(),  // LCOV_EXCL_LINE
        " state symbols at F",
        level);
  }  // LCOV_EXCL_LINE
  const StateCube modelCube =
      extractStateCube(solver, variables, badStateSupport, 0);
  const StateCube cube = reduceSolvedCubeByTernarySimulation(
      solver,
      variables,
      {{badFormula, nullptr, true}},
      modelCube,
      supportCache);
  if (pdrStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: bad cube level=", level,
        " source=ternary_simulation",
        " state_symbols=", badStateSupport.size(),
        " model_cube=", modelCube.size(),
        " cube=", cube.size(),
        " hash=", cubeFingerprint(cube));
  }
  return cube;
}

std::optional<StateCube> findBadCubeForFormula(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    BoolExpr* badFormula,
    const std::optional<std::vector<size_t>>& preciseBadStateSupport,
    size_t level,
    const ComplementPartnerIndex& complementPartners,
    BadCubeAssumptionCache* badCubeAssumptionCache,
    PredecessorAssumptionCache* predecessorAssumptionCache,
    PdrFormulaSupportCache* supportCache) {
  if (!preciseBadStateSupport.has_value()) {
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: bad cube support budget exhausted level=",
          level,
          " node_limit=",
          kMaxPreciseBadCubeSupportNodes);
    }
    markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);
    return std::nullopt;
  }

  // Search the current frame for a concrete state that still satisfies bad
  // after all learned blocking clauses and optional strengthening are applied.
  std::vector<size_t> computedSolverSymbols;
  const std::vector<size_t>* solverSymbolsPtr = nullptr;
  const bool useSharedFrameZeroSolver =
      level == 0 && predecessorAssumptionCache != nullptr &&
      predecessorAssumptionCache->sharedFrameZeroPredecessorSolver != nullptr &&
      predecessorAssumptionCache->sharedFrameZeroPredecessorSymbols !=
          nullptr &&
      !predecessorAssumptionCache->sharedFrameZeroPredecessorSymbols->empty();
  if (useSharedFrameZeroSolver) {
    // prepareSharedExactInitQueries() already built this complete immutable
    // surface. Borrow it rather than reconstructing exact Init for every output.
    solverSymbolsPtr =
        predecessorAssumptionCache->sharedFrameZeroPredecessorSymbols;
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: shared exact F[0] symbols reused for bad cube count=",
          solverSymbolsPtr->size());
    }
  } else {
    computedSolverSymbols = findBadQuerySymbols(
        initFormula,
        frameInvariant,
        frames,
        badFormula,
        level,
        complementPartners,
        supportCache);
    solverSymbolsPtr = &computedSolverSymbols;
  }
  const std::vector<size_t>& solverSymbols = *solverSymbolsPtr;
  const unsigned badCubeConflictLimit =
      // LCOV_EXCL_START
      problem.usesDualRailStateEncoding ? dualRailBadCubeConflictLimit() : 0;
      // LCOV_EXCL_STOP
  const size_t badCubeStatsQueryNumber = nextPdrBadCubeQueryNumber();
  const bool emitStatsForBadCubeQuery =
      shouldEmitPdrStats(badCubeStatsQueryNumber);
  BadCubeAssumptionCache* solverCache =
      !useSharedFrameZeroSolver &&
              shouldUseBadCubeSolverCache(problem, level, solverSymbols.size())
          ? badCubeAssumptionCache
          : nullptr;
  if (problem.usesDualRailStateEncoding && badCubeAssumptionCache != nullptr &&
      !useSharedFrameZeroSolver && solverCache == nullptr &&
      emitStatsForBadCubeQuery) {
    emitSecDiag(
        "SEC PDR stats: bad cube cached solver disabled query_symbols=",
        solverSymbols.size(),
        " query_limit=",
        kMaxDualRailBadCubeSolverCacheStateSymbols,
        " total_state_symbols=",
        problem.totalStateCount,
        " level=",
        level);
  }
  if (problem.usesDualRailStateEncoding && useSharedFrameZeroSolver) {
    PredecessorAssumptionSolver* solvedCache = nullptr;
    const auto badSolveStatus = solveFrameZeroBadCubeWithSharedSolver(
        *predecessorAssumptionCache, problem, solverType, transitionByState,
        complementPartners, initFormula, frameInvariant, frames, badFormula,
        solverSymbols, badCubeConflictLimit, supportCache, &solvedCache);
    if (badSolveStatus == SATSolverWrapper::SolveStatus::Unknown) {
      if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
        emitSecDiag(  // LCOV_EXCL_LINE
            "SEC PDR stats: bad cube query budget exhausted limit=",
            badCubeConflictLimit,
            " symbols=",
            solverSymbols.size(),  // LCOV_EXCL_LINE
            " level=",
            level,
            " cached_assumptions=1");
      }  // LCOV_EXCL_LINE
      markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);
      return std::nullopt;
    }
    if (badSolveStatus == SATSolverWrapper::SolveStatus::Unsat) {
      return std::nullopt;
    }
    return extractSolvedBadCubeForFormula(
        *solvedCache->solver,
        *solvedCache->variables,
        *preciseBadStateSupport,
        badFormula,
        level,
        supportCache);
  }
  if (problem.usesDualRailStateEncoding && solverCache != nullptr) {
    BadCubeAssumptionSolver* solvedCache = nullptr;
    const auto badSolveStatus = solveBadCubeWithCachedAssumption(
        *solverCache, problem, solverType, complementPartners, initFormula,
        frameInvariant, frames, level, badFormula, solverSymbols,
        badCubeConflictLimit, &solvedCache);
    if (badSolveStatus == SATSolverWrapper::SolveStatus::Unknown) {
      if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
        emitSecDiag(  // LCOV_EXCL_LINE
            "SEC PDR stats: bad cube query budget exhausted limit=",
            badCubeConflictLimit,
            " symbols=",
            solverSymbols.size(),  // LCOV_EXCL_LINE
            " level=",
            level,
            " cached_assumptions=1");
      }  // LCOV_EXCL_LINE
      markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);  // LCOV_EXCL_LINE
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    if (badSolveStatus == SATSolverWrapper::SolveStatus::Unsat) {
      return std::nullopt;
    }
    return extractSolvedBadCubeForFormula(
        *solvedCache->solver,
        *solvedCache->variables,
        *preciseBadStateSupport,
        badFormula,
        level,
        supportCache);
  }

  SATSolverWrapper solver(solverType);
  // Bad-state queries are local PDR obligations and are rebuilt repeatedly as
  // frames advance. Keep them on the PDR-local profile: small regressions such
  // as GCD can otherwise spend minutes in Kissat's speculative
  // preprocessing/probing before the actual frame query starts.
  // LCOV_EXCL_START
  solver.configureForSecPdrQuery(solverSymbols.size());
  FrameVariableStore variables(solver, solverSymbols, 1);
  complementPartners.addClauses(solver, variables, solverSymbols, 1);
  // LCOV_EXCL_STOP
  addFrameConstraints(solver, variables, initFormula, frameInvariant, frames,
                      level, 0);
  addPostBootstrapResetInputConstraints(solver, variables, problem, 0);
  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(0));
  solver.addClause({encoder.encode(badFormula)});
  SATSolverWrapper::SolveStatus badSolveStatus =
      SATSolverWrapper::SolveStatus::Sat;
  if (badCubeConflictLimit != 0) {
    // Dual-rail residual repairs can be SAT and decision-heavy even when they
    // do not accumulate many conflicts. Bound both resources so a single
    // LCOV_EXCL_START
    // uncovered output cannot dominate the whole workflow.
    // LCOV_EXCL_STOP
    badSolveStatus = solver.solveWithResourceLimits( // LCOV_EXCL_LINE
        badCubeConflictLimit, // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        /*decisionLimit=*/badCubeConflictLimit);
        // LCOV_EXCL_STOP
  } else { // LCOV_EXCL_LINE
    badSolveStatus = solver.solveStatus();
  }
  if (badSolveStatus == SATSolverWrapper::SolveStatus::Unknown) {
    if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC PDR stats: bad cube query budget exhausted limit=",
          badCubeConflictLimit,
          " symbols=",
          solverSymbols.size(),  // LCOV_EXCL_LINE
          " level=",
          level);
    }  // LCOV_EXCL_LINE
    markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  if (badSolveStatus == SATSolverWrapper::SolveStatus::Unsat) {
    return std::nullopt;
  }

  return extractSolvedBadCubeForFormula(
      solver,
      variables,
      *preciseBadStateSupport,
      badFormula,
      level,
      supportCache);
}

std::optional<StateCube> findBadCube(const KInductionProblem& problem,
                                     KEPLER_FORMAL::Config::SolverType solverType,
                                     const TransitionExprResolver& transitionByState,
                                     BoolExpr* initFormula,
                                     BoolExpr* frameInvariant,
                                     const std::vector<FrameClauses>& frames,
                                     BoolExpr* badFormula,
                                     bool decomposeOutputBad,
                                     const std::optional<std::vector<size_t>>&
                                         preciseBadStateSupport,
                                     const std::unordered_set<size_t>& stateSymbols,
                                     size_t level,
                                     const ComplementPartnerIndex& complementPartners,
                                     BadCubeAssumptionCache* badCubeAssumptionCache,
                                     PredecessorAssumptionCache* predecessorAssumptionCache,
                                     PdrFormulaSupportCache* supportCache) {
  if (!decomposeOutputBad || problem.observedOutputExprs0.size() <= 1 ||
      problem.observedOutputExprs0.size() != problem.observedOutputExprs1.size()) {
    return findBadCubeForFormula(
        problem, solverType, transitionByState, initFormula, frameInvariant,
        frames, badFormula,
        preciseBadStateSupport, level,
        complementPartners, badCubeAssumptionCache, predecessorAssumptionCache,
        supportCache);
  }

  // This is an exact decomposition of R[N] & !P: each query uses the complete
  // state and frame constraints, and the disjunction is UNSAT exactly when all
  // output mismatch terms are UNSAT. It avoids one broad unrelated SAT cone.
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {
    BoolExpr* outputBad = BoolExpr::simplify(BoolExpr::Xor(
        problem.observedOutputExprs0[output],
        problem.observedOutputExprs1[output]));
    const auto outputStateSupport = collectBoundedStateSupportSymbols(
        outputBad, kMaxPreciseBadCubeSupportNodes, 0, stateSymbols);
    if (auto cube = findBadCubeForFormula(
            problem, solverType, transitionByState, initFormula, frameInvariant,
            frames, outputBad, outputStateSupport, level,
            complementPartners, badCubeAssumptionCache,
            predecessorAssumptionCache, supportCache);
        cube.has_value()) {
      return cube;
    }
    if (hasPdrBudgetExhaustion()) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
  }
  return std::nullopt;
}

struct PredecessorQueryOutcome {
  bool hasPredecessor = false;
  std::optional<StateCube> predecessor;
};

PredecessorQueryOutcome findPredecessorCube(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const StateCube& targetCube,
    bool excludeTargetOnCurrentFrame,
    PredecessorQueryPurpose purpose,
    const ComplementPartnerIndex& complementPartners,
    PredecessorAssumptionCache* predecessorAssumptionCache = nullptr,
    size_t* predecessorQueryBudget = nullptr,
    PdrFormulaSupportCache* supportCache = nullptr,
    PredecessorAssumptionCache* narrowGeneralizationProbeCache = nullptr) {
  // This is the one-step predecessor query at the heart of PDR: does some
  // state in F[level] transition into the target cube on the next frame?
  std::optional<PredecessorQueryResultKey> exactCacheKey;
  std::optional<PredecessorQueryResultKey> stableUnsatCacheKey;
  const bool usePredecessorQueryResultCache =
      predecessorAssumptionCache != nullptr;
  if (usePredecessorQueryResultCache) {
    const size_t frameFingerprint = frameClausesFingerprint(frames, level);
    exactCacheKey = makeCachedPredecessorQueryResultKey(
        *predecessorAssumptionCache,
        problem,
        transitionByState,
        initFormula,
        frameInvariant,
        level,
        frameFingerprint,
        excludeTargetOnCurrentFrame,
        targetCube);
    stableUnsatCacheKey = *exactCacheKey;
    stableUnsatCacheKey->frameFingerprint = 0;
    if (const auto cached = cachedPredecessorQueryResult(
            *predecessorAssumptionCache, *exactCacheKey,
            *stableUnsatCacheKey);
        cached.has_value()) {
      if (shouldEmitFrequentPdrStats()) {
        emitSecDiag(
            "SEC PDR stats: predecessor result cache hit level=",
            level,
            " has_predecessor=",
            cached->hasPredecessor ? 1 : 0,
            " has_model=",
            cached->hasPredecessorModel ? 1 : 0,
            " shared_f0=",
            level == 0 &&
                    predecessorAssumptionCache
                            ->sharedFrameZeroQueryResultStore != nullptr
                ? 1
                : 0);
      }
      if (cached->hasPredecessor) {
        if (!predecessorQueryNeedsModel(purpose) ||
            cached->hasPredecessorModel) {
          return {
              true,
              cached->hasPredecessorModel
                  ? std::optional<StateCube>(cached->predecessor)
                  : std::nullopt};
        }
        // A status-only query deliberately did not retain a SAT model. Figure
        // 6 still performs the exact solve when recursive blocking needs one.
      }
      if (!cached->hasPredecessor) {
        return {}; // LCOV_EXCL_LINE
      }
    }
    if (const auto cachedCore = cachedPredecessorUnsatCoreForTarget(
            *predecessorAssumptionCache, *stableUnsatCacheKey, targetCube);
        cachedCore.has_value()) {
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: predecessor unsat-core cache hit level=",
            level,
            " target_cube=",
            targetCube.size(),
            " core_cube=",
            cachedCore->size(),
            " target_hash=",
            cubeFingerprint(targetCube),
            " core_hash=",
            cubeFingerprint(*cachedCore),
            " shared_f0=",
            level == 0 &&
                    predecessorAssumptionCache
                            ->sharedFrameZeroQueryResultStore != nullptr
                ? 1
                : 0);
      }
      rememberPredecessorQueryResult(
          *predecessorAssumptionCache,
          *exactCacheKey,
          *stableUnsatCacheKey,
          std::nullopt,
          &*cachedCore);
      return {};
    }
  }
  if (!consumePdrPredecessorQueryBudget(predecessorQueryBudget)) {
    return {};  // LCOV_EXCL_LINE
  }
  const size_t statsQueryNumber = nextPdrPredecessorQueryNumber();
  const bool emitStatsForQuery = shouldEmitPdrStats(statsQueryNumber);
  PredecessorTargetSurface uncachedTargetSurface;
  const PredecessorTargetSurface* targetSurface = nullptr;
  if (predecessorAssumptionCache != nullptr) {
    targetSurface = &predecessorTargetSurfaceFor(
        *predecessorAssumptionCache, problem, transitionByState,
        complementPartners, targetCube, uncachedTargetSurface);
  } else {
    uncachedTargetSurface =
        buildPredecessorTargetSurface(
            problem, transitionByState, complementPartners, targetCube);
    targetSurface = &uncachedTargetSurface;
  }
  const std::vector<size_t>& encodedTargets = targetSurface->encodedTargets;
  const std::vector<size_t>& transitionSupportSymbols =
      targetSurface->transitionSupportSymbols;
  const size_t transitionEncodingNodes =
      targetSurface->transitionEncodingNodes;
  if (problem.usesDualRailStateEncoding) {
    const size_t encodingNodeLimit = dualRailPredecessorEncodingNodeLimit();
    const size_t nodeHintTargetLimit =
        dualRailPredecessorNodeHintTargetLimit();
    const size_t encodingSupportLimit =
        dualRailPredecessorEncodingSupportLimit();
    const bool unknownNodeCount =
        transitionEncodingNodes == 0 &&
        encodedTargets.size() > nodeHintTargetLimit;
    if (unknownNodeCount || transitionEncodingNodes > encodingNodeLimit ||
        transitionSupportSymbols.size() > encodingSupportLimit) {
      if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
        emitSecDiag(  // LCOV_EXCL_LINE
            "SEC PDR stats: predecessor encoding budget exhausted targets=",
            encodedTargets.size(),
            " nodes=",
            transitionEncodingNodes,
            " node_limit=",
            encodingNodeLimit,
            " node_hint_target_limit=",
            nodeHintTargetLimit,
            " transition_support=",
            transitionSupportSymbols.size(),
            " support_limit=",
            encodingSupportLimit,
            " level=",
            level);
      }  // LCOV_EXCL_LINE
      markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);  // LCOV_EXCL_LINE
      return {};  // LCOV_EXCL_LINE
    }
  }

  // Section V materializes only the transition cone needed by this query.
  // Ternary simulation immediately removes every state outside that cone, so
  // asking SAT to assign those absent variables first is redundant. F[level],
  // the target transition functions, and all domain relations they mention
  // remain exact; this is existential CNF construction, not model reduction.
  const std::vector<size_t>& predecessorSymbols =
      targetSurface->predecessorSymbols;
  PredecessorAssumptionCache* solverCache =
      shouldUsePredecessorSolverCache(
          problem, level, problem.totalStateCount)
          ? predecessorAssumptionCache
          : nullptr;
  std::vector<size_t> computedSolverSymbols;
  const std::vector<size_t>* solverSymbolsPtr = nullptr;
  if (level == 0 && solverCache != nullptr &&
      solverCache->sharedFrameZeroPredecessorSymbols != nullptr &&
      !solverCache->sharedFrameZeroPredecessorSymbols->empty()) {
    // prepareSharedExactInitQueries() includes the complete source symbol
    // surface. Borrow that exact sorted vector instead of rebuilding Init's
    // multi-million-symbol support for every predecessor cube.
    solverSymbolsPtr = solverCache->sharedFrameZeroPredecessorSymbols;
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: shared exact F[0] predecessor symbols reused count=",
          solverSymbolsPtr->size());
    }
  } else {
    computedSolverSymbols = predecessorCurrentFrameQuerySymbols(
        problem,
        initFormula,
        frameInvariant,
        frames,
        level,
        targetCube,
        excludeTargetOnCurrentFrame,
        *targetSurface,
        complementPartners,
        // Exact symbol-surface preparation is bounded independently from SAT
        // solver retention and is needed to measure the final query width.
        predecessorAssumptionCache,
        supportCache);
    solverSymbolsPtr = &computedSolverSymbols;
  }
  const std::vector<size_t>& solverSymbols = *solverSymbolsPtr;
  // A whole-chip design can still produce a tiny exact IC3 obligation. Base
  // higher-frame cache retention on that SAT surface after it is known, not on
  // unrelated state outside the query cone.
  if (solverCache == nullptr && predecessorAssumptionCache != nullptr &&
      shouldUsePredecessorSolverCache(
          problem, level, solverSymbols.size())) {
    solverCache = predecessorAssumptionCache;
  }
  const std::vector<size_t>& cachedSolverSymbols =
      predecessorAssumptionCacheSymbols(
          transitionByState,
          level,
          solverSymbols,
          solverCache);
  const unsigned predecessorConflictLimit =
      problem.usesDualRailStateEncoding
          ? dualRailPredecessorConflictLimit(purpose)
          : 0;
  const unsigned predecessorDecisionLimit =
      problem.usesDualRailStateEncoding
          ? dualRailPredecessorDecisionLimit(purpose)
          : std::numeric_limits<unsigned>::max();
  if (emitStatsForQuery) {
    emitSecDiag(
        "SEC PDR stats: predecessor #", statsQueryNumber,
        " level=", level,
        " target_cube=", targetCube.size(),
        " target_hash=", cubeFingerprint(targetCube),
        " encoded_targets=", encodedTargets.size(),
        " transition_support=", transitionSupportSymbols.size(),
        " predecessor_symbols=", predecessorSymbols.size(),
        " solver_symbols=", solverSymbols.size(),
        " cached_solver_symbols=", cachedSolverSymbols.size(),
        " conflict_limit=", predecessorConflictLimit,
        " decision_limit=", predecessorDecisionLimit,
        " frame_clauses=",
        level < frames.size() ? frames[level].clauses.size() : 0,
        " exclude_target=", excludeTargetOnCurrentFrame ? 1 : 0,
        " purpose=", predecessorQueryPurposeName(purpose));
  }
  if (problem.usesDualRailStateEncoding && predecessorAssumptionCache != nullptr &&
      solverCache == nullptr && emitStatsForQuery) {
    emitSecDiag(
        "SEC PDR stats: predecessor cached solver disabled query_symbols=",
        solverSymbols.size(),
        " query_limit=",
        kMaxDualRailPredecessorSolverCacheStateSymbols,
        " total_state_symbols=",
        problem.totalStateCount,
        " level=",
        level);
  }
  if (problem.usesDualRailStateEncoding &&
      purpose == PredecessorQueryPurpose::GeneralizeBlocker &&
      level > 0 &&
      solverCache != nullptr &&
      narrowGeneralizationProbeCache != nullptr) {
    const std::vector<size_t>& narrowSolverSymbols =
        predecessorAssumptionCacheSymbols(
            transitionByState,
            level,
            solverSymbols,
            narrowGeneralizationProbeCache);
    const unsigned probeConflictLimit =
        boundedNarrowGeneralizationProbeLimit(
            predecessorConflictLimit,
            kNarrowGeneralizationProbeConflictLimit);
    const unsigned probeDecisionLimit =
        boundedNarrowGeneralizationProbeLimit(
            predecessorDecisionLimit,
            kNarrowGeneralizationProbeDecisionLimit);
    StateCube narrowUnsatCore;
    const auto narrowStatus = solvePredecessorCubeWithCachedAssumptions(
        *narrowGeneralizationProbeCache,
        problem,
        solverType,
        transitionByState,
        complementPartners,
        initFormula,
        frameInvariant,
        frames,
        level,
        *targetSurface,
        narrowSolverSymbols,
        excludeTargetOnCurrentFrame,
        purpose,
        probeConflictLimit,
        probeDecisionLimit,
        supportCache,
        nullptr,
        &narrowUnsatCore);
    if (narrowStatus.has_value() &&
        *narrowStatus != SATSolverWrapper::SolveStatus::Unknown) {
      const bool hasPredecessor =
          *narrowStatus == SATSolverWrapper::SolveStatus::Sat;
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: predecessor #",
            statsQueryNumber,
            " narrow_generalization_probe result=",
            hasPredecessor ? "sat" : "unsat",
            " symbols=",
            narrowSolverSymbols.size(),
            " persistent_request_symbols=",
            cachedSolverSymbols.size(),
            " conflict_limit=",
            probeConflictLimit,
            " decision_limit=",
            probeDecisionLimit);
      }
      if (exactCacheKey.has_value() &&
          stableUnsatCacheKey.has_value() &&
          predecessorAssumptionCache != nullptr) {
        if (hasPredecessor) {
          rememberPredecessorQueryResult(
              *predecessorAssumptionCache,
              *exactCacheKey,
              *stableUnsatCacheKey,
              std::nullopt,
              nullptr,
              /*predecessorExistsWithoutModel=*/true);
        } else {
          const StateCube* narrowUnsatCorePtr =
              narrowUnsatCore.empty() ? nullptr : &narrowUnsatCore;
          rememberPredecessorQueryResult(
              *predecessorAssumptionCache,
              *exactCacheKey,
              *stableUnsatCacheKey,
              std::nullopt,
              narrowUnsatCorePtr);
        }
      }
      return {hasPredecessor, std::nullopt};
    }
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: predecessor #",
          statsQueryNumber,
          " narrow_generalization_probe result=unknown "
          "fallback=persistent symbols=",
          narrowSolverSymbols.size(),
          " persistent_request_symbols=",
          cachedSolverSymbols.size(),
          " conflict_limit=",
          probeConflictLimit,
          " decision_limit=",
          probeDecisionLimit);
    }
  }
  if (solverCache != nullptr) {
    PredecessorAssumptionSolver* solvedPredecessorCache = nullptr;
    StateCube cachedUnsatCore;
    const auto solveStart = std::chrono::steady_clock::now();
    const auto cachedStatus = solvePredecessorCubeWithCachedAssumptions(
        *solverCache, problem, solverType, transitionByState,
        complementPartners, initFormula, frameInvariant, frames, level,
        *targetSurface,
        cachedSolverSymbols, excludeTargetOnCurrentFrame,
        purpose,
        predecessorConflictLimit, predecessorDecisionLimit,
        supportCache,
        &solvedPredecessorCache, &cachedUnsatCore);
    if (emitStatsForQuery) {
      const auto solveMicros =
          std::chrono::duration_cast<std::chrono::microseconds>(
              std::chrono::steady_clock::now() - solveStart)
              .count();
      emitSecDiag(
          "SEC PDR stats: predecessor #", statsQueryNumber,
          " cached_query_us=", solveMicros,
          " cached_assumptions=1");
    }
    if (cachedStatus.has_value() &&
        *cachedStatus == SATSolverWrapper::SolveStatus::Unknown) {
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: predecessor query budget exhausted limit=",
            predecessorConflictLimit,
            " decision_limit=",
            predecessorDecisionLimit,
            " symbols=",
            cachedSolverSymbols.size(),
            " target_cube=",
            targetCube.size(),
            " observed_outputs=",
            problem.observedOutputExprs0.size(),
            " purpose=",
            predecessorQueryPurposeName(purpose),
            " level=",
            level,
            " cached_assumptions=1");
      }
      markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);
      return {};
    }
    if (cachedStatus.has_value()) {
      if (*cachedStatus == SATSolverWrapper::SolveStatus::Unsat) {
        if (emitStatsForQuery) {
          emitSecDiag(
              "SEC PDR stats: predecessor #", statsQueryNumber,
              " result=unsat cached_assumptions=1");
        }
        if (exactCacheKey.has_value() && stableUnsatCacheKey.has_value()) {
          const StateCube* cachedUnsatCorePtr =
              cachedUnsatCore.empty() ? nullptr : &cachedUnsatCore;
          rememberPredecessorQueryResult(
              *predecessorAssumptionCache,
              *exactCacheKey,
              *stableUnsatCacheKey,
              std::nullopt,
              cachedUnsatCorePtr);
        }
        return {};
      }
      if (*cachedStatus == SATSolverWrapper::SolveStatus::Sat &&
          solvedPredecessorCache != nullptr) {
        const bool extractModel = predecessorQueryNeedsModel(purpose);
        if (emitStatsForQuery) {
          emitSecDiag(
              "SEC PDR stats: predecessor #", statsQueryNumber,
              " result=sat cached_assumptions=1 model_extracted=",
              extractModel ? 1 : 0,
              " purpose=",
              predecessorQueryPurposeName(purpose));
        }
        if (!extractModel) {
          if (exactCacheKey.has_value() && stableUnsatCacheKey.has_value()) {
            rememberPredecessorQueryResult(
                *predecessorAssumptionCache,
                *exactCacheKey,
                *stableUnsatCacheKey,
                std::nullopt,
                nullptr,
                /*predecessorExistsWithoutModel=*/true);
          }
          return {true, std::nullopt};
        }
        StateCube predecessor = extractSolvedPredecessorCube(
            *solvedPredecessorCache->solver,
            *solvedPredecessorCache->variables,
            predecessorSymbols,
            targetSurface->ternaryRoots,
            supportCache);
        if (emitStatsForQuery) {
          emitSecDiag(
              "SEC PDR stats: predecessor #", statsQueryNumber,
              " predecessor_cube=", predecessor.size(),
              " predecessor_hash=", cubeFingerprint(predecessor));
          if (supportCache != nullptr) {
            supportCache->emitMappedTernarySupportStats();
          }
        }
        if (exactCacheKey.has_value() && stableUnsatCacheKey.has_value()) {
          rememberPredecessorQueryResult(
              *predecessorAssumptionCache,
              *exactCacheKey,
              *stableUnsatCacheKey,
              std::optional<StateCube>(predecessor));
        }
        return {true, std::move(predecessor)};
      }
    }
  }
  SATSolverWrapper solver(solverType);
  solver.configureForSecPdrQuery(solverSymbols.size());
  FrameVariableStore variables(solver, solverSymbols, 1);
  complementPartners.addClauses(solver, variables, solverSymbols, 1);
  addFrameConstraints(solver, variables, initFormula, frameInvariant, frames,
                      level, 0);
  addSafeFramePropertyConstraint(
      solver, variables, problem, level, supportCache, 0);
  addPostBootstrapResetInputConstraints(solver, variables, problem, 0);
  // Encode only the next-state equations needed to decide the requested target
  // cube. This keeps one local PDR obligation from materializing the entire
  // design transition relation.
  std::unordered_map<size_t, int> transitionLeafLits;
  addTransitionConstraintsForTargetGroups(
      solver,
      variables,
      transitionByState,
      0,
      targetSurface->transitionGroups,
      transitionSupportSymbols,
      &transitionLeafLits);
  if (excludeTargetOnCurrentFrame) {
    addNegatedCubeClause(solver, variables, targetCube, 0);
  }
  SATSolverWrapper::SolveStatus predecessorSolveStatus =
      SATSolverWrapper::SolveStatus::Sat;
  if (problem.usesDualRailStateEncoding) {
    // Predecessor queries are local PDR obligations. A limit hit is not a
    // proof of UNSAT, so dual-rail mode turns it into an inconclusive leaf
    // instead of letting one hard residual output dominate the regress run.
    predecessorSolveStatus = solver.solveWithResourceLimits(
        predecessorConflictLimit,
        predecessorDecisionLimit);
  } else {
    predecessorSolveStatus = solver.solveStatus();
  }
  if (predecessorSolveStatus == SATSolverWrapper::SolveStatus::Unknown) {
    if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC PDR stats: predecessor query budget exhausted limit=",
          predecessorConflictLimit,
          " decision_limit=",
          predecessorDecisionLimit,
          " symbols=",
          solverSymbols.size(),
          " target_cube=",
          targetCube.size(),
          " observed_outputs=",
          problem.observedOutputExprs0.size(),
          " purpose=",
          predecessorQueryPurposeName(purpose),
          " level=",
          level);
    }  // LCOV_EXCL_LINE
    markPdrBudgetExhausted(PdrBudgetExhaustion::LocalQuery);  // LCOV_EXCL_LINE
    return {};  // LCOV_EXCL_LINE
  }
  const bool hasPredecessor =
      predecessorSolveStatus == SATSolverWrapper::SolveStatus::Sat;
  if (emitStatsForQuery) {
    emitSecDiag(
        "SEC PDR stats: predecessor #", statsQueryNumber,
        " result=", hasPredecessor ? "sat" : "unsat",
        hasPredecessor ? " model_extracted=" : "",
        hasPredecessor
            ? (predecessorQueryNeedsModel(purpose) ? "1" : "0")
            : "",
        hasPredecessor ? " purpose=" : "",
        hasPredecessor ? predecessorQueryPurposeName(purpose) : "");
  }
  if (!hasPredecessor) {
    if (exactCacheKey.has_value() && stableUnsatCacheKey.has_value() &&
        predecessorAssumptionCache != nullptr) {
      rememberPredecessorQueryResult(
          *predecessorAssumptionCache,
          *exactCacheKey,
          *stableUnsatCacheKey,
          std::nullopt);
    }
    return {};
  }
  if (!predecessorQueryNeedsModel(purpose)) {
    if (exactCacheKey.has_value() && stableUnsatCacheKey.has_value() &&
        predecessorAssumptionCache != nullptr) {
      rememberPredecessorQueryResult(
          *predecessorAssumptionCache,
          *exactCacheKey,
          *stableUnsatCacheKey,
          std::nullopt,
          nullptr,
          /*predecessorExistsWithoutModel=*/true);
    }
    return {true, std::nullopt};
  }
  StateCube predecessor = extractSolvedPredecessorCube(
      solver,
      variables,
      predecessorSymbols,
      targetSurface->ternaryRoots,
      supportCache);
  if (emitStatsForQuery) {
    emitSecDiag(
        "SEC PDR stats: predecessor #", statsQueryNumber,
        " predecessor_cube=", predecessor.size(),
        " predecessor_hash=", cubeFingerprint(predecessor));
    if (supportCache != nullptr) {
      supportCache->emitMappedTernarySupportStats();
    }
  }
  if (exactCacheKey.has_value() && stableUnsatCacheKey.has_value() &&
      predecessorAssumptionCache != nullptr) {
    rememberPredecessorQueryResult(
        *predecessorAssumptionCache,
        *exactCacheKey,
        *stableUnsatCacheKey,
        std::optional<StateCube>(predecessor));
  }
  return {true, std::move(predecessor)};
}

InitIntersectionAssumptionSolver& getInitIntersectionAssumptionSolver(
    PredecessorAssumptionCache& cache,
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    BoolExpr* initFormula) {
  auto& solver = cache.initIntersectionSolver;
  const KInductionProblem* problemIdentity = &problem;
  if (solver != nullptr && solver->problem == problemIdentity &&
      solver->initFormula == initFormula &&
      solver->requestedSolverType == solverType) {
    return *solver;
  }

  auto cached = std::make_unique<InitIntersectionAssumptionSolver>();
  cached->problem = problemIdentity;
  cached->initFormula = initFormula;
  cached->requestedSolverType = solverType;
  const std::vector<size_t> solverSymbols =
      initIntersectionSymbols(problem, initFormula);
  cached->solver = std::make_unique<SATSolverWrapper>(
      SATSolverWrapper::assumptionSolverTypeFor(solverType));
  cached->solver->configureForSecPdrQuery(solverSymbols.size());
  cached->variables =
      std::make_unique<FrameVariableStore>(*cached->solver, solverSymbols, 1);
  std::optional<ComplementPartnerIndex> localStateRelations;
  if (cache.stateRelations == nullptr) {
    localStateRelations.emplace(problem); // LCOV_EXCL_LINE
  }
  const ComplementPartnerIndex& stateRelations = cache.stateRelations != nullptr
                                                     ? *cache.stateRelations
                                                     : *localStateRelations;
  stateRelations.addClauses(*cached->solver, *cached->variables, solverSymbols,
                            1);
  FrameFormulaEncoder encoder(*cached->solver,
                              cached->variables->makeLeafLits(0));
  cached->solver->addClause({encoder.encode(initFormula)});

  solver = std::move(cached);
  return *solver;
}

PredecessorAssumptionSolver* exactFrameZeroSolverForInitIntersection(
    PredecessorAssumptionCache& cache,
    BoolExpr* initFormula,
    const StateCube& cube) {
  PredecessorAssumptionSolver* solver = nullptr;
  if (cache.sharedFrameZeroPredecessorSolver != nullptr) {
    solver = cache.sharedFrameZeroPredecessorSolver->get();
  } else if (const auto local = cache.solversByLevel.find(0);
             local != cache.solversByLevel.end()) {
    solver = local->second.get();
  }
  if (solver == nullptr || solver->key.level != 0 ||
      solver->key.initFormula != initFormula) {
    return nullptr;
  }
  for (const auto& literal : cube) {
    if (!solver->variables->hasSymbol(literal.symbol)) {
      return nullptr;
    }
  }
  return solver;
}

bool solveInitIntersectionWithAssumptions(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    std::unordered_map<StateCube, bool, StateCubeHash>& resultByCube,
    const StateCube& cube) {
  const auto cachedResult = resultByCube.find(cube);
  if (cachedResult != resultByCube.end()) {
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: exact F[0] init intersection cache reused cube=",
          cube.size());
    }
    return cachedResult->second;
  }
  std::vector<int> assumptions;
  assumptions.reserve(cube.size());
  for (const auto& literal : cube) {
    if (!variables.hasSymbol(literal.symbol)) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR init-intersection encoding missing state symbol " +
          std::to_string(literal.symbol));  // LCOV_EXCL_LINE
    }
    const int satLiteral = variables.getLiteral(literal.symbol, 0);
    assumptions.push_back(literal.value ? satLiteral : -satLiteral);
  }
  const bool intersects = solver.solveWithAssumptions(assumptions);
  resultByCube.emplace(cube, intersects);
  return intersects;
}

bool cubeIntersectsInit(const KInductionProblem& problem,
                        KEPLER_FORMAL::Config::SolverType solverType,
                        BoolExpr* initFormula,
                        const StateCube& cube,
                        PredecessorAssumptionCache* cache) {
  // Definite conflicts can avoid a SAT call. Otherwise Figure 6 requires the
  // exact F[0] = I query; absence of a known conflict is not reachability.
  if (cubeContradictsKnownInitFacts(
          problem, cube, cache != nullptr ? cache->initFacts : nullptr)) {
    return false;
  }

  PredecessorAssumptionCache localCache;
  PredecessorAssumptionCache& activeCache =
      cache != nullptr ? *cache : localCache;
  if (auto* frameZeroSolver = exactFrameZeroSolverForInitIntersection(
          activeCache, initFormula, cube);
      frameZeroSolver != nullptr) {
    if (shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: shared exact F[0] solver used for init intersection");
    }
    return solveInitIntersectionWithAssumptions(
        *frameZeroSolver->solver,
        *frameZeroSolver->variables,
        frameZeroSolver->initIntersectionResultByCube,
        cube);
  }
  auto& cached = getInitIntersectionAssumptionSolver(
      activeCache, problem, solverType, initFormula);
  return solveInitIntersectionWithAssumptions(
      *cached.solver, *cached.variables, cached.resultByCube, cube);
}

std::optional<StateCube> growCoreOutsideInit(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    BoolExpr* initFormula,
    const StateCube& core,
    const StateCube& targetCube,
    PredecessorAssumptionCache* cache) {
  StateCube candidate = core;
  if (!cubeIntersectsInit(
          problem, solverType, initFormula, candidate, cache)) {
    return candidate;
  }

  // Pdr_ManReduceClause in the authors' ABC implementation strengthens a
  // failed-assumption core with literals from the original cube until it no
  // longer overlaps Init. Strengthening preserves the solved Q2 implication.
  for (const auto& literal : targetCube) {
    if (findCubeLiteralValue(candidate, literal.symbol).has_value()) {
      continue;
    }
    candidate.push_back(literal);
    normalizeCube(candidate);
    if (!cubeIntersectsInit(
            problem, solverType, initFormula, candidate, cache)) {
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: predecessor core kept outside init core=",
            core.size(),
            "->",
            candidate.size(),
            " target=",
            targetCube.size());
      }
      return candidate;
    }
  }
  return std::nullopt;
}

class BlockedCubeReductionChecker {
 public:
  BlockedCubeReductionChecker(
      const KInductionProblem& problem,
      KEPLER_FORMAL::Config::SolverType solverType,
      const TransitionExprResolver& transitionByState,
      BoolExpr* initFormula,
      BoolExpr* frameInvariant,
      const std::vector<FrameClauses>& frames,
      size_t level,
      PredecessorAssumptionCache* predecessorAssumptionCache,
      PredecessorAssumptionCache* narrowGeneralizationProbeCache,
      const ComplementPartnerIndex& complementPartners,
      size_t* predecessorQueryBudget,
      PdrFormulaSupportCache* supportCache)
      : problem_(problem),
        solverType_(solverType),
        transitionByState_(transitionByState),
        initFormula_(initFormula),
        frameInvariant_(frameInvariant),
        frames_(frames),
        level_(level),
        predecessorAssumptionCache_(predecessorAssumptionCache),
        narrowGeneralizationProbeCache_(narrowGeneralizationProbeCache),
        complementPartners_(complementPartners),
        predecessorQueryBudget_(predecessorQueryBudget),
        supportCache_(supportCache) {}

  std::optional<StateCube> cachedCore(const StateCube& cube) const {
    if (predecessorAssumptionCache_ == nullptr) {
      return std::nullopt;
    }
    const auto core = cachedPredecessorUnsatCoreForCube(
        *predecessorAssumptionCache_,
        problem_,
        transitionByState_,
        initFormula_,
        frameInvariant_,
        frames_,
        level_ - 1,
        cube,
        /*excludeTargetOnCurrentFrame=*/true);
    if (!core.has_value() || core->size() >= cube.size()) {
      return std::nullopt;
    }
    return growCoreOutsideInit(
        problem_,
        solverType_,
        initFormula_,
        *core,
        cube,
        predecessorAssumptionCache_);
  }

  std::optional<StateCube> generalize(const StateCube& reduced) const {
    if (reduced.empty() ||
        cubeIntersectsInit(
            problem_,
            solverType_,
            initFormula_,
            reduced,
            predecessorAssumptionCache_)) {
      return std::nullopt;
    }
    // Figure 7 generalizes with Q2: F[k-1] & !s & T & s'.
    const auto predecessor = findPredecessorCube(
        problem_,
        solverType_,
        transitionByState_,
        initFormula_,
        frameInvariant_,
        frames_,
        level_ - 1,
        reduced,
        /*excludeTargetOnCurrentFrame=*/true,
        PredecessorQueryPurpose::GeneralizeBlocker,
        complementPartners_,
        predecessorAssumptionCache_,
        predecessorQueryBudget_,
        supportCache_,
        narrowGeneralizationProbeCache_);
    if (hasPdrBudgetExhaustion() || predecessor.hasPredecessor) {
      return std::nullopt;
    }
    if (const auto core = cachedCore(reduced); core.has_value()) {
      return core;
    }
    return reduced;
  }

 private:
  const KInductionProblem& problem_;
  KEPLER_FORMAL::Config::SolverType solverType_;
  const TransitionExprResolver& transitionByState_;
  BoolExpr* initFormula_ = nullptr;
  BoolExpr* frameInvariant_ = nullptr;
  const std::vector<FrameClauses>& frames_;
  size_t level_ = 0;
  PredecessorAssumptionCache* predecessorAssumptionCache_ = nullptr;
  PredecessorAssumptionCache* narrowGeneralizationProbeCache_ = nullptr;
  const ComplementPartnerIndex& complementPartners_;
  size_t* predecessorQueryBudget_ = nullptr;
  PdrFormulaSupportCache* supportCache_ = nullptr;
};

StateCube generalizeBlockedCube(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const StateCube& cube,
    PredecessorAssumptionCache* predecessorAssumptionCache,
    const ComplementPartnerIndex& complementPartners,
    size_t* predecessorQueryBudget,
    PdrFormulaSupportCache* supportCache) {
  // Figure 7's literal-removal checks share one exact local SAT context. Its
  // lifetime ends with this cube, so unrelated output cones cannot widen it.
  PredecessorAssumptionCache narrowGeneralizationProbeCache;
  narrowGeneralizationProbeCache.stateRelations = &complementPartners;
  BlockedCubeReductionChecker reductionChecker(
      problem,
      solverType,
      transitionByState,
      initFormula,
      frameInvariant,
      frames,
      level,
      predecessorAssumptionCache,
      &narrowGeneralizationProbeCache,
      complementPartners,
      predecessorQueryBudget,
      supportCache);

  // Figure 7 first uses the failed assumptions from solveRelative, then visits
  // each remaining literal once with the same Q2 query. Keep that original
  // order while the cube shrinks; retrying an earlier SAT removal would be the
  // stronger non-monotone generalization that Section VI-B found unhelpful.
  StateCube generalized = cube;
  if (const auto core = reductionChecker.cachedCore(generalized);
      core.has_value()) {
    generalized = *core;
  }

  const StateCube literalsToTry = generalized;
  size_t checks = 0;
  for (const auto& literal : literalsToTry) {
    if (generalized.size() <= 1) {
      break;
    }
    const auto current = std::lower_bound(
        generalized.begin(), generalized.end(), literal, cubeLiteralLess);
    if (current == generalized.end() || !(*current == literal)) {
      continue;
    }
    StateCube reduced = generalized;
    reduced.erase(
        reduced.begin() + std::distance(generalized.begin(), current));
    ++checks;
    const auto result = reductionChecker.generalize(reduced);
    if (hasPdrBudgetExhaustion()) {
      retireGeneralizationStatusQ2Selectors(predecessorAssumptionCache);
      return generalized;
    }
    if (!result.has_value()) {
      continue;
    }
    generalized = *result;
  }

  if (generalized.size() != cube.size() &&
      shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: generalized blocked cube level=",
        level,
        " size=",
        cube.size(),
        "->",
        generalized.size(),
        " checks=",
        checks);
  }
  // Figure 7 has finished consuming this local family of Q2 assumptions.
  // Retire only status selectors; recursively blocked cubes remain cached.
  retireGeneralizationStatusQ2Selectors(predecessorAssumptionCache);
  return generalized;
}

bool framesConverged(const FrameClauses& lhs, const FrameClauses& rhs) {
  if (lhs.clauses.size() != rhs.clauses.size()) {
    return false;
  }
  for (const auto& clause : lhs.clauses) {
    if (!frameHasSubsumingClause(rhs, clause)) {
      return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
  }
  for (const auto& clause : rhs.clauses) {
    if (!frameHasSubsumingClause(lhs, clause)) {
      return false;  // LCOV_EXCL_LINE
    }
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  return true;
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
bool obligationAlreadyBlocked(const std::vector<FrameClauses>& frames,
// LCOV_EXCL_STOP
                              const ProofObligation& obligation) {
  return frameHasSubsumingClause(frames[obligation.level], clauseFromCube(obligation.cube));
}  // LCOV_EXCL_LINE

StateCube generalizeInitExcludedCube(const KInductionProblem& problem,  // LCOV_EXCL_LINE
                                     KEPLER_FORMAL::Config::SolverType solverType,
                                     BoolExpr* initFormula,
                                     const StateCube& cube,
                                     PredecessorAssumptionCache* cache) {
  // Ordinary Init can also be a relational frontier made of equality facts.
  // When a predecessor violates that frontier, learn a generalized
  // LCOV_EXCL_STOP
  // F[0] clause immediately instead of relying on many small seed clauses to
  // LCOV_EXCL_START
  // rediscover adjacent impossible cubes one at a time.
  StateCube candidate = cube;  // LCOV_EXCL_LINE
  size_t index = 0;  // LCOV_EXCL_LINE
  size_t attempts = 0;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  while (index < candidate.size() &&  // LCOV_EXCL_LINE
         attempts < kMaxInitExcludedCubeGeneralizationAttempts) {  // LCOV_EXCL_LINE
    ++attempts;  // LCOV_EXCL_LINE
    StateCube reduced = candidate;  // LCOV_EXCL_LINE
    reduced.erase(reduced.begin() + static_cast<std::ptrdiff_t>(index));  // LCOV_EXCL_LINE
    if (!cubeIntersectsInit(  // LCOV_EXCL_LINE
            problem, solverType, initFormula, reduced, cache)) {  // LCOV_EXCL_LINE
      candidate = std::move(reduced);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    ++index;  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return candidate;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

ProofObligationKey proofObligationKey(const ProofObligation& obligation) {
  ProofObligationKey key;
  key.level = obligation.level;
  key.badFrame = obligation.badFrame;
  key.cube = obligation.cube;
  return key;
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

class ProofObligationLowerPriority {
 public:
  bool operator()(const ProofObligation& lhs,
                  const ProofObligation& rhs) const {
    // priority_queue places the element for which this relation is false on
    // top. Reverse the existing Figure 6 priority predicate exactly.
    return detail::pdrProofObligationPriorityLess(rhs.level, rhs.sequence,
                                                  lhs.level, lhs.sequence);
  }
};

class ProofObligationQueue {
 public:
  bool empty() const { return queue_.empty(); }

  bool enqueue(ProofObligation obligation) {
    if (!queuedKeys_.insert(proofObligationKey(obligation)).second) {
      return false; // LCOV_EXCL_LINE
    }
    obligation.sequence = nextSequence_++;
    queue_.push(std::move(obligation));
    return true;
  }

  ProofObligation pop() {
    ProofObligation obligation = queue_.top();
    queue_.pop();
    queuedKeys_.erase(proofObligationKey(obligation));
    return obligation;
  }

  void enqueueNext(const ProofObligation& obligation, size_t rootLevel) {
    if (obligation.level >= rootLevel) {
      return;
    }
    // Figure 6 requeues next(s); advance badFrame too so the path suffix length
    // remains unchanged for counterexample reporting.
    ProofObligation next = obligation;
    ++next.level;
    ++next.badFrame;
    if (enqueue(std::move(next)) && shouldEmitFrequentPdrStats()) {
      emitSecDiag(
          "SEC PDR stats: proof obligation requeued level=", obligation.level,
          "->", obligation.level + 1, " bad_frame=", obligation.badFrame + 1);
    }
  }

 private:
  std::priority_queue<ProofObligation, std::vector<ProofObligation>,
                      ProofObligationLowerPriority>
      queue_;
  std::unordered_set<ProofObligationKey, ProofObligationKeyHash> queuedKeys_;
  size_t nextSequence_ = 0;
};

void learnBlockedObligation(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    std::vector<FrameClauses>& frames,
    size_t rootLevel,
    const ComplementPartnerIndex& complementPartners,
    PredecessorAssumptionCache& predecessorAssumptionCache,
    size_t* predecessorQueryBudget,
    PdrFormulaSupportCache* supportCache,
    const ProofObligation& obligation) {
  StateCube cube = generalizeBlockedCube(
      problem, solverType, transitionByState, initFormula, frameInvariant,
      frames, obligation.level, obligation.cube, &predecessorAssumptionCache,
      complementPartners, predecessorQueryBudget, supportCache);
  size_t learnedLevel = obligation.level;
  // Figure 6 repeatedly applies solveRelative(next(z)) and keeps the highest
  // frame where the blocker is relatively inductive.
  while (learnedLevel + 1 < rootLevel) {
    const auto predecessor = findPredecessorCube(
        problem, solverType, transitionByState, initFormula, frameInvariant,
        frames, learnedLevel, cube,
        /*excludeTargetOnCurrentFrame=*/true,
        PredecessorQueryPurpose::LiftBlocker,
        complementPartners,
        &predecessorAssumptionCache, predecessorQueryBudget,
        supportCache);
    if (hasPdrBudgetExhaustion() || predecessor.hasPredecessor) {
      break;
    }
    ++learnedLevel;
  }
  addClauseToFrames(frames, clauseFromCube(cube), learnedLevel);
  if (learnedLevel != obligation.level &&
      shouldEmitFrequentPdrStats()) {
    emitSecDiag(
        "SEC PDR stats: blocked cube lifted level=",
        obligation.level,
        "->",
        learnedLevel,
        " cube=",
        cube.size());
  }
}

bool blockProofObligations(const KInductionProblem& problem,
                           KEPLER_FORMAL::Config::SolverType solverType,
                           const TransitionExprResolver& transitionByState,
                           BoolExpr* initFormula,
                           BoolExpr* frameInvariant,
                           std::vector<FrameClauses>& frames,
                           const InitFactIndex& initFacts,
                           const StateCube& badCube,
                           size_t rootLevel,
                           size_t& badFrame,
                           const ComplementPartnerIndex& complementPartners,
                           PredecessorAssumptionCache& predecessorAssumptionCache,
                           size_t* predecessorQueryBudget,
                           PdrFormulaSupportCache* supportCache) {
  // This is the paper's recursive blocking idea expressed as an explicit queue
  // so we do not depend on deep recursion for large obligation stacks.
  ProofObligationQueue queue;
  (void)queue.enqueue(ProofObligation{badCube, rootLevel, rootLevel});

  while (!queue.empty()) {
    const ProofObligation obligation = queue.pop();

    if (obligationAlreadyBlocked(frames, obligation)) {
      continue;  // LCOV_EXCL_LINE
    }

    if (obligation.level == 0) {
      if (const auto conflictCube =
              knownInitConflictCube(initFacts, obligation.cube);
          conflictCube.has_value()) {
        // When the cube visibly contradicts a structured exact Init fact,
        // learn only that conflict instead of a wide SAT-model cube;
        // LCOV_EXCL_STOP
        // this keeps large ASIC output slices from rediscovering the same
        // LCOV_EXCL_START
        // state equality violation thousands of times.
        if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
          emitSecDiag(  // LCOV_EXCL_LINE
              "SEC PDR stats: known init conflict ",
              "cube=", obligation.cube.size(),  // LCOV_EXCL_LINE
              " core=", conflictCube->size(),  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
              " bad_frame=", obligation.badFrame,  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              " hash=", cubeFingerprint(*conflictCube));  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
        }  // LCOV_EXCL_LINE
        addClauseToFrame(frames[0], clauseFromCube(*conflictCube));  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      if (!cubeIntersectsInit(
              problem,
              solverType,
              initFormula,
              obligation.cube,
              &predecessorAssumptionCache)) {
        const StateCube generalizedCube = generalizeInitExcludedCube(  // LCOV_EXCL_LINE
            problem,  // LCOV_EXCL_LINE
            solverType,  // LCOV_EXCL_LINE
            initFormula,  // LCOV_EXCL_LINE
            obligation.cube,  // LCOV_EXCL_LINE
            &predecessorAssumptionCache);  // LCOV_EXCL_LINE
        addClauseToFrame(frames[0], clauseFromCube(generalizedCube));  // LCOV_EXCL_LINE
        continue;
      }  // LCOV_EXCL_LINE
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: counterexample candidate reached init ",
            "bad_frame=", obligation.badFrame,
            // LCOV_EXCL_START
            " cube=", obligation.cube.size());
      }
      // LCOV_EXCL_STOP
      badFrame = obligation.badFrame;
      return false;
    }

    // Q2 is sound only for cubes outside Init. Figure 6 makes this invariant
    // explicit before every relative-induction query.
    if (cubeIntersectsInit(
            problem,
            solverType,
            initFormula,
            obligation.cube,
            &predecessorAssumptionCache)) {
      if (pdrStatsEnabled()) {
        emitSecDiag(
            "SEC PDR stats: counterexample candidate intersects init ",
            "level=", obligation.level,
            " bad_frame=", obligation.badFrame,
            " cube=", obligation.cube.size());
      }
      badFrame = obligation.badFrame;
      return false;
    }

    const auto predecessor = findPredecessorCube(
        problem,
        solverType,
        transitionByState,
        initFormula,
        frameInvariant,
        frames,
        obligation.level - 1,
        obligation.cube,
        /*excludeTargetOnCurrentFrame=*/true,
        PredecessorQueryPurpose::BlockObligation,
        complementPartners,
        &predecessorAssumptionCache,
        predecessorQueryBudget,
        supportCache);
    if (hasPdrBudgetExhaustion()) {
      return true;  // LCOV_EXCL_LINE
    }
    if (!predecessor.hasPredecessor) {
      learnBlockedObligation(
          problem, solverType, transitionByState, initFormula, frameInvariant,
          frames, rootLevel, complementPartners, predecessorAssumptionCache,
          predecessorQueryBudget, supportCache, obligation);
      if (hasPdrBudgetExhaustion()) {
        return true;  // LCOV_EXCL_LINE
      }
      queue.enqueueNext(obligation, rootLevel);
      continue;
    }
    if (!predecessor.predecessor.has_value()) {  // LCOV_EXCL_LINE
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "PDR blocking predecessor model was not extracted");  // LCOV_EXCL_LINE
    }
    ProofObligation predecessorObligation{
        *predecessor.predecessor, obligation.level - 1,
        obligation.badFrame};
    (void)queue.enqueue(obligation);
    (void)queue.enqueue(std::move(predecessorObligation));
  }

  return true;
}

std::vector<StateClause> buildSeedClauses(const KInductionProblem& problem,
                                          const InitFactIndex& initFacts) {
  (void)problem;
  (void)initFacts;
  std::vector<StateClause> seedClauses;
  // Cross-design internal state equality seeds are forbidden.
  return seedClauses;
}

BoolExpr* selectPdrFrameInvariant(const KInductionProblem& problem,
                                  BoolExpr* initFormula,
                                  KEPLER_FORMAL::Config::SolverType solverType) {
  FormulaSupportCache invariantSupportCache;
  auto validateCandidate =
      [&](const char* label, BoolExpr* candidate) -> BoolExpr* {
    if (candidate == nullptr) {
      if (pdrStatsEnabled()) {
        emitSecDiag("SEC PDR stats: frame invariant ", label, " unavailable");
      }
      return nullptr;
    }

    const bool initImpliesCandidate =
        initialFrontierImplies(initFormula, candidate, solverType);
    const bool inductive =
        initImpliesCandidate &&
        isInductiveInvariant(
            problem, candidate, solverType, invariantSupportCache);
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: frame invariant ", label,
          " support=", candidate->getSupportVars().size(),
          " init=", initImpliesCandidate ? "pass" : "fail",
          " inductive=", inductive ? "pass" : "fail");
    }
    if (!initImpliesCandidate || !inductive) {
      return nullptr;
    }
    return candidate;
  };

  // PDR may reuse a validated public SEC strengthening lemma as a frame fact,
  // but it must not build a frame invariant from cross-design internal state
  // equalities.
  BoolExpr* sharedStrengthening =
      selectValidatedStrengtheningInvariant(problem, initFormula, solverType);
  if (BoolExpr* strengthenedInvariant =
          validateCandidate("shared_strengthening", sharedStrengthening)) {
    if (isSecDiagEnabled()) {
      emitSecDiag(
          "SEC diag: PDR using validated SEC strengthening frame invariant");
    }
    return strengthenedInvariant;
  }
  return nullptr;
}

void propagateClauses(const KInductionProblem& problem,
                      KEPLER_FORMAL::Config::SolverType solverType,
                      const TransitionExprResolver& transitionByState,
                      BoolExpr* initFormula,
                      BoolExpr* frameInvariant,
                      std::vector<FrameClauses>& frames,
                      size_t maxLevel,
                      const ComplementPartnerIndex& complementPartners,
                      PredecessorAssumptionCache* predecessorAssumptionCache,
                      size_t* predecessorQueryBudget,
                      PdrFormulaSupportCache* supportCache) {
  // Standard PDR propagation: if F[i] /\ T implies a clause on the next frame,
  // move that clause forward into F[i+1].
  for (size_t level = 1; level <= maxLevel; ++level) {
    const auto snapshot = frames[level].clauses;
    for (const auto& clause : snapshot) {
      // Only propagate clauses that are not already known to hold on the next frame,
      // otherwise we would be doing redundant work and risking over-blocking by
      // adding the same clause again after generalization.
      if (frameHasSubsumingClause(frames[level + 1], clause)) {
        continue;
      }
      const StateCube violatingCube = cubeFromClauseNegation(clause);
      // A clause is only safe to propagate if it does not block a real bad path, so check
      // whether any predecessor of the negated cube survives in the current frame. If not, the
      // clause can be added to the next frame without risking over-blocking.
      const auto predecessor = findPredecessorCube(
          problem,
          solverType,
          transitionByState,
          initFormula,
          frameInvariant,
          frames,
          level,
          violatingCube,
          false,
          PredecessorQueryPurpose::PropagateClause,
          complementPartners,
          predecessorAssumptionCache,
          predecessorQueryBudget,
          supportCache);
      if (hasPdrBudgetExhaustion()) {
        // Figure 9 propagation is opportunistic: only a proved-UNSAT query
        // moves the clause. UNKNOWN leaves this clause in its current frame;
        // it must not abort otherwise exact blocking work for the whole output.
        if (pdrStatsEnabled()) {
          emitSecDiag(
              "SEC PDR stats: propagation left clause in frame level=",
              level,
              " clause_literals=",
              clause.size());
        }
        resetPdrBudgetExhaustion();
        continue;
      }
      if (!predecessor.hasPredecessor) {
        addClauseToFrame(frames[level + 1], clause);
      }
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
  }
}

bool isSecPdrTraceEnabled() {
  return std::getenv("KEPLER_SEC_PDR_TRACE") != nullptr;
}

std::string formatSymbolForPdrTrace(size_t symbol) {
  if (symbol == 0) {
    return "FALSE";  // LCOV_EXCL_LINE
  }
  if (symbol == 1) {
    return "TRUE";  // LCOV_EXCL_LINE
  }
  return "x" + std::to_string(symbol);
}

std::string formatCubeForPdrTrace(const StateCube& cube) {
  std::ostringstream oss;
  oss << "{";
  for (size_t i = 0; i < cube.size(); ++i) {
    if (i != 0) {
      // LCOV_EXCL_START
      oss << ", ";
    }
    // LCOV_EXCL_STOP
    oss << formatSymbolForPdrTrace(cube[i].symbol) << "=" << (cube[i].value ? "1" : "0");
  }
  oss << "}";
  return oss.str();
}

std::string formatClauseForPdrTrace(const StateClause& clause) {
  std::ostringstream oss;
  oss << "(";
  for (size_t i = 0; i < clause.size(); ++i) {
    if (i != 0) {
      oss << " OR ";  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    if (!clause[i].positive) {
      oss << "!";
    }
    oss << formatSymbolForPdrTrace(clause[i].symbol);
  }
  oss << ")";
  return oss.str();
}

std::string formatFramesForPdrTrace(const std::vector<FrameClauses>& frames) {
  std::ostringstream oss;
  for (size_t level = 0; level < frames.size(); ++level) {
    oss << "  F[" << level << "]: ";
    if (level == 0) {
      oss << "Init";
    }
    oss << "\n";
    if (frames[level].clauses.empty()) {
      oss << "    <empty>\n";
      continue;
    }
    for (const auto& clause : frames[level].clauses) {
      oss << "    " << formatClauseForPdrTrace(clause) << "\n";
    }
  }
  return oss.str();
}

void emitPdrTrace(std::string_view label, const std::string& body) {
  if (!isSecPdrTraceEnabled()) {
    return;
  }
  emitSecDiag("SEC PDR trace: ", label, "\n", body);
}

void emitPdrTraceProblem(const KInductionProblem& problem) {
  if (!isSecPdrTraceEnabled()) {
    return;
  }
  // Full formula formatting recursively walks every transition/property DAG.
  // That is useful for small debug tests, but on ASIC-size SEC problems it can
  // allocate gigabytes before PDR starts.  Keep the expensive string build
  // strictly behind the explicit PDR trace flag.
  emitSecDiag("SEC PDR trace: problem\n", formatKInductionProblemForDebug(problem));
}

void emitPdrTraceFrames(std::string_view label,
                        const std::vector<FrameClauses>& frames) {
  if (!isSecPdrTraceEnabled()) {
    return;
  }
  emitSecDiag("SEC PDR trace: ", label, "\n", formatFramesForPdrTrace(frames));
}

BoolExpr* appendAssignmentFormula(
    BoolExpr* formula,
    const std::vector<std::pair<size_t, bool>>& assignments) {
  for (const auto& [symbol, value] : assignments) {
    BoolExpr* variable = BoolExpr::Var(symbol);
    formula = BoolExpr::And(
        formula, value ? variable : BoolExpr::Not(variable));
  }
  return formula;
}

BoolExpr* makeBalancedConjunction(std::vector<BoolExpr*> terms) {
  if (terms.empty()) {
    return BoolExpr::createTrue();
  }
  while (terms.size() > 1) {
    std::vector<BoolExpr*> next;
    next.reserve((terms.size() + 1) / 2);
    for (size_t index = 0; index < terms.size(); index += 2) {
      next.push_back(index + 1 < terms.size()
                         ? BoolExpr::And(terms[index], terms[index + 1])
                         : terms[index]);
    }
    terms = std::move(next);
  }
  return terms.front();
}

class ExactPdrBootstrapInitBuilder {
 public:
  explicit ExactPdrBootstrapInitBuilder(const KInductionProblem& problem)
      : problem_(problem),
        transitionByState_(problem),
        stateSymbols_(problem.combinedStateSymbols()),
        symbolAtFrame_(problem.resetBootstrapCycles + 1),
        remapMemoByFrame_(problem.resetBootstrapCycles) {
    std::sort(stateSymbols_.begin(), stateSymbols_.end());
    stateSymbols_.erase(
        std::unique(stateSymbols_.begin(), stateSymbols_.end()),
        stateSymbols_.end());
  }

  BoolExpr* build() {
    collectReservedSymbolsAndTransitions();
    initializeStateSymbols();

    std::vector<BoolExpr*> terms;
    terms.reserve(
        transitions_.size() * problem_.resetBootstrapCycles +
        problem_.initialStateAssignments.size() + 32);
    // The reset prefix starts from the design's actual ternary initialization:
    // known registers use 01/10 and resetless registers use X=11.
    addInitialStateAssignments(terms);
    for (size_t frame = 0; frame < problem_.resetBootstrapCycles; ++frame) {
      addResetInputAssignments(frame, /*asserted=*/true, terms);
      addStateDomainRelations(frame, terms);
      addTransitionRelation(frame, terms);
    }
    addStateDomainRelations(problem_.resetBootstrapCycles, terms);

    // The paper requires F[0] = I.  The final reset frame is therefore part of
    // I as well: reset is deasserted and the observed frontier property is the
    // same one used by the exact bounded SEC base query.
    addResetInputAssignments(
        problem_.resetBootstrapCycles, /*asserted=*/false, terms);
    if (problem_.usesResetBootstrapObservationFrontier()) {
      terms.push_back(problem_.property);
    }
    return makeBalancedConjunction(std::move(terms));
  }

 private:
  void reserveFormulaSymbols(BoolExpr* formula) {
    if (formula == nullptr) {
      return;
    }
    const auto support = formula->getSupportVars();
    reservedSymbols_.insert(support.begin(), support.end());
  }

  void reserveSupportSymbols(const std::set<size_t>& support) {
    reservedSymbols_.insert(support.begin(), support.end());
  }

  void reserveAssignments(
      const std::vector<std::pair<size_t, bool>>& assignments) {
    for (const auto& [symbol, /*value*/ _] : assignments) {
      reservedSymbols_.insert(symbol);
    }
  }

  void collectReservedSymbolsAndTransitions() {
    reservedSymbols_.insert(problem_.allSymbols.begin(), problem_.allSymbols.end());
    reservedSymbols_.insert(stateSymbols_.begin(), stateSymbols_.end());
    reserveAssignments(problem_.resetBootstrapInputs);
    reserveAssignments(problem_.initialStateAssignments);
    reserveFormulaSymbols(problem_.property);
    reserveFormulaSymbols(problem_.bad);

    transitions_.reserve(stateSymbols_.size());
    for (const size_t stateSymbol : stateSymbols_) {
      if (!transitionByState_.contains(stateSymbol)) {
        continue;
      }
      BoolExpr* transition = transitionByState_.at(stateSymbol);
      transitions_.emplace_back(stateSymbol, transition);
      // The resolver already caches exact support for eager and lazy
      // transitions. Reuse it here so reset-prefix construction does not walk
      // each large materialized dual-rail DAG once for every use.
      reserveSupportSymbols(transitionByState_.support(stateSymbol));
    }
  }

  size_t allocateFreshSymbol() {
    while (reservedSymbols_.find(nextFreshSymbol_) != reservedSymbols_.end()) {
      ++nextFreshSymbol_;
    }
    const size_t symbol = nextFreshSymbol_++;
    reservedSymbols_.insert(symbol);
    return symbol;
  }

  void initializeStateSymbols() {
    const size_t finalFrame = problem_.resetBootstrapCycles;
    for (size_t frame = 0; frame <= finalFrame; ++frame) {
      auto& frameSymbols = symbolAtFrame_[frame];
      frameSymbols.reserve(stateSymbols_.size());
      for (const size_t stateSymbol : stateSymbols_) {
        frameSymbols.emplace(
            stateSymbol,
            frame == finalFrame ? stateSymbol : allocateFreshSymbol());
      }
    }
  }

  size_t mappedSymbol(size_t frame, size_t symbol) {
    auto& frameSymbols = symbolAtFrame_.at(frame);
    if (const auto found = frameSymbols.find(symbol);
        found != frameSymbols.end()) {
      return found->second;
    }
    const size_t mapped = allocateFreshSymbol();
    frameSymbols.emplace(symbol, mapped);
    return mapped;
  }

  BoolExpr* remapAtFrame(BoolExpr* formula,
                         size_t frame,
                         const std::set<size_t>& support) {
    for (const size_t symbol : support) {
      if (symbol >= 2) {
        static_cast<void>(mappedSymbol(frame, symbol));
      }
    }
    return remapBoolExprVariables(
        formula, symbolAtFrame_.at(frame), remapMemoByFrame_.at(frame));
  }

  void addInitialStateAssignments(std::vector<BoolExpr*>& terms) {
    for (const auto& [symbol, value] : problem_.initialStateAssignments) {
      BoolExpr* variable = BoolExpr::Var(mappedSymbol(0, symbol));
      terms.push_back(value ? variable : BoolExpr::Not(variable));
    }
  }

  void addResetInputAssignments(size_t frame,
                                bool asserted,
                                std::vector<BoolExpr*>& terms) {
    for (const auto& [symbol, assertedValue] : problem_.resetBootstrapInputs) {
      const bool value = asserted ? assertedValue : !assertedValue;
      BoolExpr* variable =
          frame == problem_.resetBootstrapCycles
              ? BoolExpr::Var(symbol)
              : BoolExpr::Var(mappedSymbol(frame, symbol));
      terms.push_back(value ? variable : BoolExpr::Not(variable));
    }
  }

  void addComplementRelations(
      size_t frame,
      const std::vector<std::pair<size_t, size_t>>& pairs,
      std::vector<BoolExpr*>& terms) {
    for (const auto& [primary, complement] : pairs) {
      terms.push_back(makeEqualityExpr(
          BoolExpr::Var(mappedSymbol(frame, complement)),
          BoolExpr::Not(BoolExpr::Var(mappedSymbol(frame, primary)))));
    }
  }

  void addEqualityRelations(
      size_t frame,
      const std::vector<std::pair<size_t, size_t>>& pairs,
      std::vector<BoolExpr*>& terms) {
    for (const auto& [lhs, rhs] : pairs) {
      terms.push_back(makeEqualityExpr(
          BoolExpr::Var(mappedSymbol(frame, lhs)),
          BoolExpr::Var(mappedSymbol(frame, rhs))));
    }
  }

  void addDualRailValidity(size_t frame,
                           std::vector<BoolExpr*>& terms) {
    // Historical reset states belong to the same exact dual-rail domain as
    // F[0]; (may-be-one, may-be-zero) = (0, 0) is not a valid state.
    for (const auto& rails : problem_.dualRailStatePairs) {
      terms.push_back(BoolExpr::Or(
          BoolExpr::Var(mappedSymbol(frame, rails.mayBeOne)),
          BoolExpr::Var(mappedSymbol(frame, rails.mayBeZero))));
    }
  }

  void addStateDomainRelations(size_t frame,
                               std::vector<BoolExpr*>& terms) {
    addComplementRelations(frame, problem_.complementedStatePairs0, terms);
    addComplementRelations(frame, problem_.complementedStatePairs1, terms);
    addEqualityRelations(frame, problem_.sameFrameStateEqualityPairs0, terms);
    addEqualityRelations(frame, problem_.sameFrameStateEqualityPairs1, terms);
    addDualRailValidity(frame, terms);
  }

  void addTransitionRelation(size_t frame,
                             std::vector<BoolExpr*>& terms) {
    for (const auto& [stateSymbol, transition] : transitions_) {
      terms.push_back(makeEqualityExpr(
          BoolExpr::Var(mappedSymbol(frame + 1, stateSymbol)),
          remapAtFrame(
              transition, frame, transitionByState_.support(stateSymbol))));
    }
  }

  const KInductionProblem& problem_;
  TransitionExprResolver transitionByState_;
  std::vector<size_t> stateSymbols_;
  std::vector<std::pair<size_t, BoolExpr*>> transitions_;
  std::unordered_set<size_t> reservedSymbols_;
  size_t nextFreshSymbol_ = 2;
  std::vector<std::unordered_map<size_t, size_t>> symbolAtFrame_;
  std::vector<std::unordered_map<BoolExpr*, BoolExpr*>> remapMemoByFrame_;
};

BoolExpr* buildExactPdrInitFormula(
    const KInductionProblem& problem,
    PDRExactInitCache::Impl* sharedExactInit) {
  if (sharedExactInit != nullptr && sharedExactInit->initFormula != nullptr) {
    if (pdrStatsEnabled()) {
      emitSecDiag("SEC PDR stats: shared exact F[0] cache reused");
    }
    return sharedExactInit->initFormula;
  }

  const KInductionProblem& source =
      sharedExactInit != nullptr ? *sharedExactInit->sourceProblem : problem;
  BoolExpr* initFormula = nullptr;
  if (source.resetBootstrapCycles != 0) {
    // PDR requires F[0] to be the initial-state predicate. For SEC that starts
    // after reset, this predicate is the exact reset transition image.
    initFormula = ExactPdrBootstrapInitBuilder(source).build();
  } else {
    BoolExpr* init = source.initialCondition != nullptr
                         ? source.initialCondition
                         : BoolExpr::createTrue();
    initFormula = BoolExpr::simplify(
        appendAssignmentFormula(init, source.initialStateAssignments));
  }
  if (sharedExactInit != nullptr) {
    sharedExactInit->initFormula = initFormula;
    if (pdrStatsEnabled()) {
      emitSecDiag("SEC PDR stats: shared exact F[0] cache built");
    }
  }
  return initFormula;
}

}  // namespace

PDREngine::PDREngine(const KInductionProblem& problem,
                     KEPLER_FORMAL::Config::SolverType solverType,
                     size_t maxPredecessorQueries,
                     std::shared_ptr<PDRExactInitCache> exactInitCache)
    : problem_(problem),
      solverType_(solverType),
      maxPredecessorQueries_(maxPredecessorQueries),
      exactInitCache_(std::move(exactInitCache)) {}

PDRResult PDREngine::run(size_t maxFrames) const {
  return runWithQueryLimits(maxFrames, problem_.property, nullptr);
}

PDRResult PDREngine::run(size_t maxFrames, BoolExpr* property) const {
  return runWithQueryLimits(maxFrames, property, nullptr);
}

PDRResult PDREngine::run(size_t maxFrames,
                         BoolExpr* property,
                         const PDRQueryLimits& queryLimits) const {
  return runWithQueryLimits(maxFrames, property, &queryLimits);
}

PDRResult PDREngine::runWithQueryLimits(
    size_t maxFrames,
    BoolExpr* property,
    const PDRQueryLimits* queryLimits) const {
  if (property == nullptr) {
    return {PDRStatus::Inconclusive, 0};  // LCOV_EXCL_LINE
  }
  // Batch probes use smaller SAT limits only to decide whether to split. A
  // limit hit remains UNKNOWN, and singleton leaves run with the normal limits.
  const ScopedPdrQueryLimits scopedQueryLimits(queryLimits);
  const bool usesDefaultProperty = property == problem_.property;
  BoolExpr* normalizedProperty = BoolExpr::simplify(property);
  BoolExpr* normalizedBad =
      BoolExpr::simplify(BoolExpr::Not(normalizedProperty));
  std::optional<KInductionProblem> alternateProblem;
  const KInductionProblem* runProblem = &problem_;
  const bool canUseOriginalProblem =
      usesDefaultProperty && problem_.property == normalizedProperty &&
      problem_.bad == normalizedBad;
  if (!canUseOriginalProblem) {
    // Normal SEC output batches already contain their selected property.  Copy
    // the large immutable model only for the alternate-property API or an
    // unusual caller whose stored bad root is not the normalized complement.
    alternateProblem.emplace(problem_);
    alternateProblem->property = normalizedProperty;
    alternateProblem->bad = normalizedBad;
    if (!usesDefaultProperty) {
      // Alternate targets are independent PDR safety properties. Do not
      // inherit a target-specific induction hypothesis from normal SEC.
      alternateProblem->inductionProperty = nullptr;
      alternateProblem->inductionBad = nullptr;
    }
    runProblem = &*alternateProblem;
  }

  // Build the SEC startup frontier once so every frame query shares the same
  // interpretation of reset/bootstrap and frame-0 equality constraints.
  resetPdrBudgetExhaustion();
  setPdrPredecessorQueryLimit(maxPredecessorQueries_);
  emitPdrTraceProblem(*runProblem);
  PDRExactInitCache::Impl* sharedExactInit = nullptr;
  if (exactInitCache_ != nullptr &&
      exactInitCache_->impl_->matches(problem_, solverType_)) {
    sharedExactInit = exactInitCache_->impl_.get();
  }
  BoolExpr* initFormula =
      buildExactPdrInitFormula(problem_, sharedExactInit);
  if (initFormula == nullptr) {
    if (pdrStatsEnabled()) {
      emitSecDiag(
          "SEC PDR stats: inconclusive reason=exact_f0_unavailable ",
          "reset_bootstrap_cycles=", problem_.resetBootstrapCycles,
          " bootstrap_assignments=",
          problem_.bootstrapStateAssignments.size(),
          " state_symbols=", problem_.combinedStateSymbols().size());
    }
    return {PDRStatus::Inconclusive, 0};  // LCOV_EXCL_LINE
  }

  // PDR still establishes convergence through its own frame/blocking loop, but
  // it may reuse a validated public SEC strengthening lemma as a frame
  // invariant after checking init coverage and transition preservation.
  BoolExpr* frameInvariant =
      selectPdrFrameInvariant(problem_, initFormula, solverType_);

  std::unique_ptr<TransitionExprResolver> localTransitionByState;
  std::unique_ptr<ComplementPartnerIndex> localStateRelations;
  std::shared_ptr<PdrFormulaSupportCache> localFormulaSupportCache;
  TransitionExprResolver* transitionByStatePtr = nullptr;
  ComplementPartnerIndex* complementPartnersPtr = nullptr;
  PdrFormulaSupportCache* formulaSupportCachePtr = nullptr;
  if (sharedExactInit != nullptr) {
    transitionByStatePtr = &sharedExactInit->transitionByState;
    complementPartnersPtr = &sharedExactInit->stateRelations;
    if (sharedExactInit->formulaSupportCache == nullptr) {
      sharedExactInit->formulaSupportCache =
          std::make_shared<PdrFormulaSupportCache>();
    }
    formulaSupportCachePtr = sharedExactInit->formulaSupportCache.get();
    ++sharedExactInit->immutableMetadataUses;
    if (pdrStatsEnabled()) {
      emitSecDiag("SEC PDR stats: immutable model metadata ",
                  sharedExactInit->immutableMetadataUses == 1 ? "built"
                                                              : "reused",
                  " use=", sharedExactInit->immutableMetadataUses);
    }
  } else {
    localTransitionByState =
        std::make_unique<TransitionExprResolver>(*runProblem);
    localStateRelations = std::make_unique<ComplementPartnerIndex>(*runProblem);
    localFormulaSupportCache = std::make_shared<PdrFormulaSupportCache>();
    transitionByStatePtr = localTransitionByState.get();
    complementPartnersPtr = localStateRelations.get();
    formulaSupportCachePtr = localFormulaSupportCache.get();
  }
  TransitionExprResolver& transitionByState = *transitionByStatePtr;
  ComplementPartnerIndex& complementPartners = *complementPartnersPtr;
  PdrFormulaSupportCache& formulaSupportCache = *formulaSupportCachePtr;
  if (sharedExactInit != nullptr) {
    prepareSharedExactInitQueries(
        *sharedExactInit,
        initFormula,
        complementPartners,
        &formulaSupportCache);
  }
  // The bad predicate is the same for every frame query. Cache its support too
  // so repeated checks do not walk the combined mismatch formula again.
  const auto preciseBadStateSupport = collectBoundedStateSupportSymbols(
      runProblem->bad, std::numeric_limits<size_t>::max(), 0,
      transitionByState.stateSymbols());
  BadCubeAssumptionCache badCubeAssumptionCache;
  PredecessorAssumptionCache predecessorAssumptionCache;
  predecessorAssumptionCache.stateRelations = &complementPartners;
  if (sharedExactInit != nullptr) {
    const size_t sharedRunId = sharedExactInit->nextHigherFrameRunId++;
    if (sharedExactInit->nextHigherFrameRunId == 0) {  // LCOV_EXCL_LINE
      sharedExactInit->nextHigherFrameRunId = 1;  // LCOV_EXCL_LINE
    }
    predecessorAssumptionCache.sharedTargetSurfaces =
        &sharedExactInit->targetSurfaces;
    predecessorAssumptionCache.sharedFrameZeroPredecessorSolver =
        &sharedExactInit->frameZeroPredecessorSolver;
    predecessorAssumptionCache.sharedFrameZeroPredecessorSymbols =
        &sharedExactInit->frameZeroPredecessorSymbols;
    predecessorAssumptionCache.sharedFrameZeroPredecessorProblem =
        sharedExactInit->sourceProblem;
    predecessorAssumptionCache.sharedFrameZeroTransitionModel =
        sharedExactInit->sourceProblem;
    predecessorAssumptionCache.sharedFrameZeroQueryResultStore =
        &sharedExactInit->frameZeroPredecessorResults;
    predecessorAssumptionCache.sharedFrameZeroQueryProblem =
        sharedExactInit->sourceProblem;
    predecessorAssumptionCache.sharedFrameZeroQueryTransition =
        &sharedExactInit->transitionByState;
    predecessorAssumptionCache.sharedHigherFrameSolverPools =
        &sharedExactInit->higherFramePredecessorSolverPools;
    predecessorAssumptionCache.sharedHigherFrameProblem =
        sharedExactInit->sourceProblem;
    predecessorAssumptionCache.sharedHigherFrameTransitionModel =
        sharedExactInit->sourceProblem;
    predecessorAssumptionCache.sharedHigherFrameRunId = sharedRunId;
    predecessorAssumptionCache.sharedHigherFrameFamilySymbols =
        &formulaSupportCache.relationClosedSupport(
            runProblem->property, complementPartners);
    predecessorAssumptionCache.usePathLocalHigherFrameSolverReuse =
        runProblem->usesStrictDualRailEqualityProperty;
  }
  size_t remainingPredecessorQueries = maxPredecessorQueries_;
  size_t* predecessorQueryBudget =
      maxPredecessorQueries_ == 0 ? nullptr : &remainingPredecessorQueries;
  std::vector<FrameClauses> frames(1);
  if (sharedExactInit != nullptr) {
    refreshReusableInvariant(*sharedExactInit, initFormula);
  }
  ReusableInvariantCandidateRecorder reusableInvariantRecorder(
      sharedExactInit, frames);
  emitPdrTraceFrames("initial_frames", frames);

  // Before growing any frame sequence, check whether exact Init itself already
  // contains a bad state.
  if (auto badCube = findBadCube(
          *runProblem, solverType_, transitionByState, initFormula,
          frameInvariant, frames, runProblem->bad, usesDefaultProperty,
          preciseBadStateSupport, transitionByState.stateSymbols(), 0,
          complementPartners, &badCubeAssumptionCache,
          &predecessorAssumptionCache, &formulaSupportCache);
      badCube.has_value()) {
    emitPdrTrace("bad_cube@F0", formatCubeForPdrTrace(*badCube));
    return {PDRStatus::Different, 0};
  }
  if (hasPdrBudgetExhaustion()) {
    return {PDRStatus::Inconclusive, 0};  // LCOV_EXCL_LINE
  }

  if (maxFrames == 0) {
    return {PDRStatus::Inconclusive, 0};
  }

  // Init/bootstrap facts are static for a PDR run. Wide dual-rail SEC problems
  // can carry tens of thousands of boot assignments, so build the lookup index
  // once instead of rebuilding it for every blocked obligation.
  std::optional<InitFactIndex> localInitFacts;
  const InitFactIndex* initFactsPtr = nullptr;
  if (sharedExactInit != nullptr) {
    if (!sharedExactInit->initFacts.has_value()) {
      sharedExactInit->initFacts.emplace(
          buildInitFactIndex(*sharedExactInit->sourceProblem));
    }
    initFactsPtr = &*sharedExactInit->initFacts;
  } else {
    localInitFacts.emplace(buildInitFactIndex(*runProblem));
    initFactsPtr = &*localInitFacts;
  }
  const InitFactIndex& initFacts = *initFactsPtr;
  predecessorAssumptionCache.initFacts = initFactsPtr;
  const auto seedClauses = buildSeedClauses(*runProblem, initFacts);
  frames.emplace_back(FrameClauses{seedClauses});
  if (sharedExactInit != nullptr) {
    injectReusableInvariantClauses(*sharedExactInit, frames[1]);
  }
  emitPdrTraceFrames("seeded_frames", frames);
  for (size_t level = 1; level <= maxFrames; ++level) {
    // Phase 1: exhaust the proof obligations created by bad states that still
    // survive in the current frontier.
    while (true) {
      const auto badCube = findBadCube(
          *runProblem, solverType_, transitionByState, initFormula,
          frameInvariant, frames, runProblem->bad, usesDefaultProperty,
          preciseBadStateSupport, transitionByState.stateSymbols(), level,
          complementPartners, &badCubeAssumptionCache,
          &predecessorAssumptionCache, &formulaSupportCache);
      if (hasPdrBudgetExhaustion()) {
        return {PDRStatus::Inconclusive, level};  // LCOV_EXCL_LINE
      }
      if (!badCube.has_value()) {
        break;
      }
      emitPdrTrace(("bad_cube@F" + std::to_string(level)).c_str(),
                   formatCubeForPdrTrace(*badCube));
      size_t badFrame = level;
      if (!blockProofObligations(
              *runProblem, solverType_, transitionByState, initFormula,
              frameInvariant, frames, initFacts, *badCube, level, badFrame,
              complementPartners, predecessorAssumptionCache,
              predecessorQueryBudget, &formulaSupportCache)) {
        if (hasPdrBudgetExhaustion()) {
          return {PDRStatus::Inconclusive, level};  // LCOV_EXCL_LINE
        }
        emitPdrTraceFrames("frames_before_counterexample", frames);
        return {PDRStatus::Different, badFrame};
      }
      if (hasPdrBudgetExhaustion()) {
        return {PDRStatus::Inconclusive, level};  // LCOV_EXCL_LINE
      }
      emitPdrTraceFrames("frames_after_blocking", frames);
    }

    // Phase 2: create the next frame, seed it with already-known startup
    // facts
    frames.emplace_back(FrameClauses{seedClauses});
    if (sharedExactInit != nullptr) {
      injectReusableInvariantClauses(
          *sharedExactInit, frames[level + 1]);
    }
    // and then push learned clauses forward.
    // We push in order to reach covergence and the condition is that that 
    // the clause is not preventing an actual bad path
    propagateClauses(*runProblem, solverType_, transitionByState, initFormula,
                     frameInvariant, frames, level, complementPartners,
                     &predecessorAssumptionCache, predecessorQueryBudget,
                     &formulaSupportCache);
    if (hasPdrBudgetExhaustion()) {
      return {PDRStatus::Inconclusive, level};  // LCOV_EXCL_LINE
    }
    emitPdrTraceFrames(("frames_after_propagation@F" + std::to_string(level)).c_str(),
                       frames);

    // Phase 3: convergence means F[i] == F[i+1], so the frame has become an
    // inductive invariant and the SEC property is proved.
    for (size_t convergenceLevel = 1; convergenceLevel <= level; ++convergenceLevel) {
      if (framesConverged(frames[convergenceLevel], frames[convergenceLevel + 1])) {
        emitPdrTraceFrames(
            ("frames_converged@F" + std::to_string(convergenceLevel)).c_str(), frames);
        return {PDRStatus::Equivalent, convergenceLevel};
      }
    }
  }
  if (pdrStatsEnabled()) {  // LCOV_EXCL_LINE
    emitSecDiag(  // LCOV_EXCL_LINE
        "SEC PDR stats: max frame budget exhausted max_frames=",
        maxFrames);
  }  // LCOV_EXCL_LINE
  return {PDRStatus::Inconclusive, maxFrames};  // LCOV_EXCL_LINE
}

}  // namespace KEPLER_FORMAL::SEC
