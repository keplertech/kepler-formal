// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "proof/TransitionExprResolver.h"

#include <algorithm>
#include <memory_resource>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "common/BoolExprUtils.h"
#include "common/PrivateProofSymbol.h"

namespace KEPLER_FORMAL::SEC {

namespace {

constexpr size_t kCachedUnionSupportTargetThreshold = 16;

size_t mapLazyTransitionSymbol(
    size_t designIndex,
    size_t localSymbol,
    std::unordered_map<size_t, size_t>& symbolMap);

class LazyDualRailVariableMapper final : public DualRailVariableMapper {  // LCOV_EXCL_LINE
 public:
  LazyDualRailVariableMapper(
      size_t designIndex,
      std::unordered_map<size_t, size_t>& binaryMap,
      const std::unordered_map<size_t, DualRailSymbolPair>& stateRails)
      : designIndex_(designIndex),
        binaryMap_(binaryMap),
        stateRails_(stateRails) {}

  DualRailBoolExpr mapVariable(size_t symbol) override {
    if (const auto stateIt = stateRails_.find(symbol);
        stateIt != stateRails_.end()) {
      return DualRailBoolExpr{
          BoolExpr::Var(stateIt->second.mayBeOne),
          BoolExpr::Var(stateIt->second.mayBeZero)};
    }

    if (symbol < 2) {
      // LCOV_EXCL_START
      return symbol == 1  // LCOV_EXCL_LINE
                 ? DualRailBoolExpr{  // LCOV_EXCL_LINE
                       BoolExpr::createTrue(), BoolExpr::createFalse()}  // LCOV_EXCL_LINE
                 : DualRailBoolExpr{  // LCOV_EXCL_LINE
                       BoolExpr::createFalse(), BoolExpr::createTrue()};  // LCOV_EXCL_LINE
                       // LCOV_EXCL_STOP
    }

    const size_t mapped = mapLazyTransitionSymbol(designIndex_, symbol, binaryMap_);
    return DualRailBoolExpr{
        BoolExpr::Var(mapped),
        BoolExpr::Not(BoolExpr::Var(mapped))};
  }

 private:
  size_t designIndex_ = 0;
  std::unordered_map<size_t, size_t>& binaryMap_;
  const std::unordered_map<size_t, DualRailSymbolPair>& stateRails_;
};

size_t countBoolExprNodes(BoolExpr* formula) {
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return 0;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  std::pmr::monotonic_buffer_resource visitedResource;
  std::pmr::unordered_set<BoolExpr*> visited{&visitedResource};
  visited.reserve(4096);
  std::vector<BoolExpr*> stack;
  stack.reserve(1024);
  stack.push_back(formula);
  while (!stack.empty()) {
    BoolExpr* node = stack.back();
    stack.pop_back();
    if (!visited.insert(node).second) {
      continue;
    }
    if (node->getRight() != nullptr) {
      stack.push_back(node->getRight());
    }
    if (node->getLeft() != nullptr) {
      stack.push_back(node->getLeft());
    }
  }
  return visited.size();
}

template <typename SymbolMapper>
std::set<size_t> collectBoolExprSupport(BoolExpr* formula,
                                        SymbolMapper&& mapSymbol) {
  std::set<size_t> support;
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return support;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // BoolExpr::getSupportVars() is intentionally stateless, but SEC/PDR asks
  // for many transition supports while validating projected candidates. Use an
  // arena-backed visited set here so each support walk avoids thousands of
  // tiny malloc/free operations before the result is stored in the resolver's
  // shared per-transition cache.
  std::pmr::monotonic_buffer_resource visitedResource;
  std::pmr::unordered_set<BoolExpr*> visited{&visitedResource};
  visited.reserve(4096);
  std::vector<BoolExpr*> stack;
  stack.reserve(1024);
  stack.push_back(formula);
  while (!stack.empty()) {
    BoolExpr* node = stack.back();
    stack.pop_back();
    if (node == nullptr || !visited.insert(node).second) {
      continue;
    }
    if (node->getOp() == Op::VAR) {
      support.insert(mapSymbol(node->getId()));
      continue;
    }
    if (node->getRight() != nullptr) {
      stack.push_back(node->getRight());
    }
    if (node->getLeft() != nullptr) {
      stack.push_back(node->getLeft());
    }
  }
  return support;
}

struct SupportVisitKey {
  BoolExpr* node = nullptr;
  std::unordered_map<size_t, size_t>* symbolMap = nullptr;
  const std::unordered_map<size_t, DualRailSymbolPair>* stateRails = nullptr;
  size_t designIndex = 0;

  bool operator==(const SupportVisitKey& other) const {
    return node == other.node && symbolMap == other.symbolMap &&
           stateRails == other.stateRails && designIndex == other.designIndex;
  }
};

struct SupportVisitKeyHash {
  size_t operator()(const SupportVisitKey& key) const {
    return std::hash<BoolExpr*>{}(key.node) ^
           (std::hash<void*>{}(key.symbolMap) << 1) ^
           (std::hash<const void*>{}(key.stateRails) << 2) ^
           (std::hash<size_t>{}(key.designIndex) << 3);
  }
};

std::set<size_t> identitySupport(BoolExpr* formula) {
  return collectBoolExprSupport(
      formula, [](size_t symbol) { return symbol; });
}

struct EncodingPostorderVisit {
  BoolExpr* node = nullptr;
  bool childrenVisited = false;
};

std::vector<BoolExpr*> buildEncodingPostorder(BoolExpr* formula) {
  std::vector<BoolExpr*> postorder;
  if (formula == nullptr) {
    return postorder; // LCOV_EXCL_LINE
  }

  // Match FrameFormulaEncoder's left-before-right iterative DFS exactly. The
  // resulting recipe stores node order only; solver literals and clauses stay
  // private to each fresh SAT query.
  std::unordered_set<BoolExpr*> encoded;
  std::vector<EncodingPostorderVisit> stack;
  stack.push_back({formula, false});
  while (!stack.empty()) {
    const EncodingPostorderVisit visit = stack.back();
    stack.pop_back();
    if (encoded.find(visit.node) != encoded.end()) {
      continue;
    }
    if (!visit.childrenVisited && visit.node->getOp() != Op::VAR) {
      stack.push_back({visit.node, true});
      if (visit.node->getRight() != nullptr) {
        stack.push_back({visit.node->getRight(), false});
      }
      if (visit.node->getLeft() != nullptr) {
        stack.push_back({visit.node->getLeft(), false});
      }
      continue;
    }
    encoded.insert(visit.node);
    postorder.push_back(visit.node);
  }
  return postorder;
}

size_t mapLazyTransitionSymbol(
    size_t designIndex,
    size_t localSymbol,
    std::unordered_map<size_t, size_t>& symbolMap) {
  if (localSymbol < 2) {
    // LCOV_EXCL_START
    return localSymbol;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (const auto mappedIt = symbolMap.find(localSymbol);
      mappedIt != symbolMap.end()) {
    return mappedIt->second;
  }

  const size_t privateSymbol =
      makePrivateProofLeafSymbol(designIndex, localSymbol);
  symbolMap.emplace(localSymbol, privateSymbol);
  return privateSymbol;
}

std::set<size_t> remappedSupport(
    BoolExpr* formula,
    size_t designIndex,
    std::unordered_map<size_t, size_t>& symbolMap) {
  return collectBoolExprSupport(formula, [&](size_t localSymbol) {
    if (localSymbol < 2) {
      // LCOV_EXCL_START
      return localSymbol;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    return mapLazyTransitionSymbol(designIndex, localSymbol, symbolMap);
  });
}

std::set<size_t> remappedDualRailSupport(
    BoolExpr* formula,
    size_t designIndex,
    std::unordered_map<size_t, size_t>& binaryMap,
    const std::unordered_map<size_t, DualRailSymbolPair>& stateRails) {
  std::vector<DualRailSymbolPair> touchedStateRails;
  std::set<size_t> support = collectBoolExprSupport(formula, [&](size_t localSymbol) {
    if (localSymbol < 2) {
      return localSymbol;  // LCOV_EXCL_LINE
    }
    if (const auto stateIt = stateRails.find(localSymbol);
        stateIt != stateRails.end()) {
      // A local state bit becomes two possible-value rails in the shared
      // dual-rail SEC problem.  Support queries need both rails even when only
      // one of the lifted transition rails is being requested.
      touchedStateRails.push_back(stateIt->second);
      return stateIt->second.mayBeOne;
    }
    return mapLazyTransitionSymbol(designIndex, localSymbol, binaryMap);
  });
  for (const auto& rails : touchedStateRails) {
    support.insert(rails.mayBeZero);
  }
  support.erase(0);
  support.erase(1);
  return support;
}

void addCollectedSupportSymbol(
    size_t symbol,
    const std::unordered_set<size_t>& knownStateSymbols,
    std::unordered_set<size_t>& stateSupport,
    std::unordered_set<size_t>& allSupport) {
  if (symbol < 2) {
    return;
  }
  allSupport.insert(symbol);
  if (knownStateSymbols.find(symbol) != knownStateSymbols.end()) {
    stateSupport.insert(symbol);
  }
}

BoolExpr* materializeLazyDualRailTransition(
    LazyTransitionSource source,
    LazyTransitionStore& store) {
  if (source.designIndex >= store.localToCombinedByDesign.size()) {
    // LCOV_EXCL_START
    throw std::runtime_error("Invalid lazy transition design index");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  LazyDualRailVariableMapper mapper(
      source.designIndex,
      store.localToCombinedByDesign[source.designIndex],
      store.dualRailStateByLocalSymbolByDesign[source.designIndex]);
  const auto lifted = buildDualRailBoolExpr(
      source.localExpr,
      mapper,
      store.dualRailRemapMemoByDesign[source.designIndex]);
  if (source.rail == LazyTransitionRail::DualRailOne) {
    return lifted.mayBeOne;
  }
  if (source.rail == LazyTransitionRail::DualRailZero) {
    return lifted.mayBeZero;
  }
  // LCOV_EXCL_START
  throw std::runtime_error("Lazy dual-rail transition requested for binary rail");  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

}  // namespace

TransitionExprResolver::TransitionExprResolver(const KInductionProblem& problem)
    : problem_(problem) {
  eagerByStateSymbol_.reserve(
      problem.transitions0.size() + problem.transitions1.size() +
      problem.auxiliaryTransitions.size());
  for (const auto& [stateSymbol, expr] : problem.transitions0) {
    eagerByStateSymbol_.emplace(stateSymbol, expr);
  }
  for (const auto& [stateSymbol, expr] : problem.transitions1) {
    eagerByStateSymbol_.emplace(stateSymbol, expr);
  }
  for (const auto& [stateSymbol, expr] : problem.auxiliaryTransitions) {
    eagerByStateSymbol_.emplace(stateSymbol, expr);
  }
}

bool TransitionExprResolver::contains(size_t stateSymbol) const {
  if (eagerByStateSymbol_.find(stateSymbol) != eagerByStateSymbol_.end()) {
    return true;
  }
  return problem_.lazyTransitions != nullptr &&
         problem_.lazyTransitions->sourceByStateSymbol.find(stateSymbol) !=
             problem_.lazyTransitions->sourceByStateSymbol.end();
}

BoolExpr* TransitionExprResolver::at(size_t stateSymbol) const {
  if (const auto eagerIt = eagerByStateSymbol_.find(stateSymbol);
      eagerIt != eagerByStateSymbol_.end()) {
    return eagerIt->second;
  }

  if (problem_.lazyTransitions == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing transition expression for state symbol " +  // LCOV_EXCL_LINE
        std::to_string(stateSymbol));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }

  auto& store = *problem_.lazyTransitions;
  if (const auto cachedIt = store.remappedByStateSymbol.find(stateSymbol);
      cachedIt != store.remappedByStateSymbol.end()) {
    return cachedIt->second;
  }

  const auto sourceIt = store.sourceByStateSymbol.find(stateSymbol);
  if (sourceIt == store.sourceByStateSymbol.end()) {
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing lazy transition expression for state symbol " +  // LCOV_EXCL_LINE
        std::to_string(stateSymbol));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }
  const LazyTransitionSource& source = sourceIt->second;
  if (source.designIndex >= store.localToCombinedByDesign.size()) {
    // LCOV_EXCL_START
    throw std::runtime_error("Invalid lazy transition design index");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (source.rail != LazyTransitionRail::Binary) {
    BoolExpr* remapped = materializeLazyDualRailTransition(source, store);
    store.remappedByStateSymbol.emplace(stateSymbol, remapped);
    return remapped;
  }

  // Populate design-private mappings for transition-only local support before
  // materializing the lazy BoolExpr.  This keeps unmodeled internal leaves
  // design-local without forcing a full transition remap during COI discovery.
  (void)remappedSupport(
      source.localExpr,
      source.designIndex,
      store.localToCombinedByDesign[source.designIndex]);
  BoolExpr* remapped = remapBoolExprVariables(
      source.localExpr,
      store.localToCombinedByDesign[source.designIndex],
      store.remapMemoByDesign[source.designIndex]);
  store.remappedByStateSymbol.emplace(stateSymbol, remapped);
  return remapped;
}

TransitionExprView TransitionExprResolver::expressionView(size_t stateSymbol) const {
  if (const auto eagerIt = eagerByStateSymbol_.find(stateSymbol);
      eagerIt != eagerByStateSymbol_.end()) {
    return TransitionExprView{eagerIt->second, nullptr};
  }

  if (problem_.lazyTransitions == nullptr) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing transition expression for state symbol " +  // LCOV_EXCL_LINE
        std::to_string(stateSymbol));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }

  const auto& store = *problem_.lazyTransitions;  // LCOV_EXCL_LINE
  const auto sourceIt = store.sourceByStateSymbol.find(stateSymbol);  // LCOV_EXCL_LINE
  if (sourceIt == store.sourceByStateSymbol.end()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing lazy transition expression for state symbol " +  // LCOV_EXCL_LINE
        std::to_string(stateSymbol));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }
  const LazyTransitionSource& source = sourceIt->second;  // LCOV_EXCL_LINE
  if (source.designIndex >= store.localToCombinedByDesign.size()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    throw std::runtime_error("Invalid lazy transition design index");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (source.rail != LazyTransitionRail::Binary) {
    return TransitionExprView{at(stateSymbol), nullptr};
  }
  // LCOV_EXCL_START
  return TransitionExprView{  // LCOV_EXCL_LINE
      source.localExpr, &store.localToCombinedByDesign[source.designIndex]};  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
}

const std::set<size_t>& TransitionExprResolver::support(size_t stateSymbol) const {
  if (const auto cachedIt = supportByStateSymbol_.find(stateSymbol);
      cachedIt != supportByStateSymbol_.end()) {
    return cachedIt->second;
  }

  // PDR asks for the same transition cone many times while blocking and
  // generalizing related cubes. BoolExpr::getSupportVars() intentionally walks
  // the DAG without owning a global cache, so the transition resolver keeps a
  // local per-proof cache and avoids rebuilding identical support sets.  For
  // lazy SEC transitions, compute support in the source model's local symbol
  // space and remap only the support IDs; remapping the full Boolean DAG just
  // to know its support dominated BlackParrot batch setup.
  if (const auto eagerIt = eagerByStateSymbol_.find(stateSymbol);
      eagerIt != eagerByStateSymbol_.end()) {
    auto [insertedIt, _] =
        supportByStateSymbol_.emplace(stateSymbol, identitySupport(eagerIt->second));
    return insertedIt->second;
  }

  if (problem_.lazyTransitions == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing transition expression for state symbol " +  // LCOV_EXCL_LINE
        std::to_string(stateSymbol));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }
  auto& store = *problem_.lazyTransitions;
  if (const auto cachedIt = store.supportByStateSymbol.find(stateSymbol);
      cachedIt != store.supportByStateSymbol.end()) {
    return cachedIt->second;
  }
  const auto sourceIt = store.sourceByStateSymbol.find(stateSymbol);
  if (sourceIt == store.sourceByStateSymbol.end()) {
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing lazy transition expression for state symbol " +  // LCOV_EXCL_LINE
        std::to_string(stateSymbol));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }
  if (sourceIt->second.designIndex >= store.localToCombinedByDesign.size()) {
    // LCOV_EXCL_START
    throw std::runtime_error("Invalid lazy transition design index");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (sourceIt->second.rail != LazyTransitionRail::Binary) {
    auto [insertedIt, _] = store.supportByStateSymbol.emplace(
        stateSymbol,
        remappedDualRailSupport(
            sourceIt->second.localExpr,
            sourceIt->second.designIndex,
            store.localToCombinedByDesign[sourceIt->second.designIndex],
            store.dualRailStateByLocalSymbolByDesign[sourceIt->second.designIndex]));
    return insertedIt->second;
  }
  auto [insertedIt, _] = store.supportByStateSymbol.emplace(
      stateSymbol,
      remappedSupport(
          sourceIt->second.localExpr,
          sourceIt->second.designIndex,
          store.localToCombinedByDesign[sourceIt->second.designIndex]));
  return insertedIt->second;
}

const std::vector<BoolExpr*>&
TransitionExprResolver::encodingPostorder(size_t stateSymbol) const {
  const TransitionExprView view = expressionView(stateSymbol);
  if (const auto cachedIt = encodingPostorderByExpr_.find(view.expr);
      cachedIt != encodingPostorderByExpr_.end()) {
    return cachedIt->second;
  }
  auto [insertedIt, _] = encodingPostorderByExpr_.emplace(
      view.expr, buildEncodingPostorder(view.expr));
  return insertedIt->second;
}

void TransitionExprResolver::collectSupportForTargets(
    const std::vector<size_t>& stateSymbols,
    const std::unordered_set<size_t>& knownStateSymbols,
    std::unordered_set<size_t>& stateSupport,
    std::unordered_set<size_t>& allSupport) const {
  if (stateSymbols.empty()) {
    return;
  }

  if (!problem_.usesDualRailStateEncoding &&
      stateSymbols.size() >= kCachedUnionSupportTargetThreshold) {
    // Large reset-frontier COI builders repeatedly ask for neighboring target
    // sets. Reuse the resolver's per-transition support cache instead of
    // rebuilding a large per-query visited table and throwing it away.
    for (const auto stateSymbol : stateSymbols) {
      for (const auto symbol : support(stateSymbol)) {
        addCollectedSupportSymbol(
            symbol, knownStateSymbols, stateSupport, allSupport);
      }
    }
    return;
  }

  // COI builders often need the union support for many transition targets in
  // the same frame. Walking each target separately repeats shared BoolExpr DAG
  // regions on ASIC designs. Traverse all requested roots together and key the
  // visited set by the source symbol map, so lazy design-local expressions are
  // still remapped exactly without forcing full DAG remapping.
  std::pmr::monotonic_buffer_resource visitedResource;
  std::pmr::unordered_set<SupportVisitKey, SupportVisitKeyHash> visited{
      &visitedResource};
  visited.reserve(std::min(
      static_cast<size_t>(1 << 20),
      std::max(static_cast<size_t>(4096), stateSymbols.size() * 1024)));
  std::vector<SupportVisitKey> stack;
  stack.reserve(stateSymbols.size());

  for (const auto stateSymbol : stateSymbols) {
    if (const auto eagerIt = eagerByStateSymbol_.find(stateSymbol);
        eagerIt != eagerByStateSymbol_.end()) {
      stack.push_back({eagerIt->second, nullptr, nullptr});
      continue;
    }

    if (problem_.lazyTransitions == nullptr) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "Missing transition expression for state symbol " +  // LCOV_EXCL_LINE
          std::to_string(stateSymbol));  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
    }
    auto& store = *problem_.lazyTransitions;  // LCOV_EXCL_LINE
    const auto sourceIt = store.sourceByStateSymbol.find(stateSymbol);  // LCOV_EXCL_LINE
    if (sourceIt == store.sourceByStateSymbol.end()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "Missing lazy transition expression for state symbol " +  // LCOV_EXCL_LINE
          std::to_string(stateSymbol));  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
    }
    const LazyTransitionSource& source = sourceIt->second;  // LCOV_EXCL_LINE
    if (source.designIndex >= store.localToCombinedByDesign.size()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      throw std::runtime_error("Invalid lazy transition design index");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (source.rail != LazyTransitionRail::Binary) {
      // Dual-rail lazy transitions for a wide bus often share most of their
      // source DAG. Walk those roots together instead of computing one full
      // support set per rail/state bit while building the same COI.
      stack.push_back(
          {source.localExpr,
           &store.localToCombinedByDesign[source.designIndex],
           &store.dualRailStateByLocalSymbolByDesign[source.designIndex],
           source.designIndex});
      continue;  // LCOV_EXCL_LINE
    }
    stack.push_back(  // LCOV_EXCL_LINE
        {source.localExpr,
         &store.localToCombinedByDesign[source.designIndex],
         nullptr,
         source.designIndex});  // LCOV_EXCL_LINE
  }

  while (!stack.empty()) {
    const SupportVisitKey current = stack.back();
    stack.pop_back();
    BoolExpr* node = current.node;
    if (node == nullptr || !visited.insert(current).second) {
      continue;
    }
    if (node->getOp() == Op::VAR) {
      size_t symbol = node->getId();
      if (current.stateRails != nullptr && symbol >= 2) {
        if (const auto stateIt = current.stateRails->find(symbol);
            stateIt != current.stateRails->end()) {
          addCollectedSupportSymbol(
              stateIt->second.mayBeOne,
              knownStateSymbols,
              stateSupport,
              allSupport);
          addCollectedSupportSymbol(
              stateIt->second.mayBeZero,
              knownStateSymbols,
              stateSupport,
              allSupport);
          continue;
        }
      }
      if (current.symbolMap != nullptr && symbol >= 2) {
        symbol = mapLazyTransitionSymbol(
            current.designIndex, symbol, *current.symbolMap);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      addCollectedSupportSymbol(
          symbol, knownStateSymbols, stateSupport, allSupport);
      continue;
    }
    if (node->getRight() != nullptr) {
      stack.push_back(
          {node->getRight(),
           current.symbolMap,
           current.stateRails,
           current.designIndex});
    }
    if (node->getLeft() != nullptr) {
      stack.push_back(
          {node->getLeft(),
           current.symbolMap,
           current.stateRails,
           current.designIndex});
    }
  }
}

size_t TransitionExprResolver::nodeCount(size_t stateSymbol) const {
  if (const auto cachedIt = nodeCountByStateSymbol_.find(stateSymbol);
      cachedIt != nodeCountByStateSymbol_.end()) {
    return cachedIt->second;
  }

  // The SAT encoder needs a rough size hint before it starts streaming clauses.
  // Computing this once per transition expression is cheaper than repeatedly
  // letting large PDR predecessor queries grow and rehash their per-query DAG
  // map while the same state transition is encoded over and over.  Lazy
  // transitions can use the source expression's node count because variable
  // remapping preserves the DAG shape.
  BoolExpr* expr = nullptr;
  if (const auto eagerIt = eagerByStateSymbol_.find(stateSymbol);
      eagerIt != eagerByStateSymbol_.end()) {
    expr = eagerIt->second;
  } else if (problem_.lazyTransitions != nullptr) {
    const auto& store = *problem_.lazyTransitions;
    if (const auto cachedIt = store.nodeCountByStateSymbol.find(stateSymbol);
        cachedIt != store.nodeCountByStateSymbol.end()) {
      nodeCountByStateSymbol_.emplace(stateSymbol, cachedIt->second);
      return cachedIt->second;
    }
    const auto sourceIt = store.sourceByStateSymbol.find(stateSymbol);
    if (sourceIt != store.sourceByStateSymbol.end()) {
      if (sourceIt->second.rail == LazyTransitionRail::Binary) {
        expr = sourceIt->second.localExpr;
      } else {
        expr = at(stateSymbol);
      }
    }
  }
  if (expr == nullptr) {
    // LCOV_EXCL_START
    expr = at(stateSymbol);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  const size_t nodeCount = countBoolExprNodes(expr);
  if (problem_.lazyTransitions != nullptr &&
      eagerByStateSymbol_.find(stateSymbol) == eagerByStateSymbol_.end()) {
    problem_.lazyTransitions->nodeCountByStateSymbol.emplace(stateSymbol, nodeCount);
  }
  auto [insertedIt, _] = nodeCountByStateSymbol_.emplace(stateSymbol, nodeCount);
  return insertedIt->second;
}

const std::unordered_set<size_t>& TransitionExprResolver::stateSymbols() const {
  if (stateSymbolsInitialized_) {
    return stateSymbols_;
  }

  // The PDR predecessor loop repeatedly asks whether a symbol belongs to the
  // combined state space. Build that lookup once per proof instead of
  // allocating the same set for every obligation.
  stateSymbols_.reserve(
      problem_.state0Symbols.size() + problem_.state1Symbols.size() +
      problem_.auxiliaryStateSymbols.size());
  stateSymbols_.insert(problem_.state0Symbols.begin(), problem_.state0Symbols.end());
  stateSymbols_.insert(problem_.state1Symbols.begin(), problem_.state1Symbols.end());
  stateSymbols_.insert(
      problem_.auxiliaryStateSymbols.begin(),
      problem_.auxiliaryStateSymbols.end());
  stateSymbolsInitialized_ = true;
  return stateSymbols_;
}

const std::unordered_map<size_t, size_t>&
TransitionExprResolver::primaryByComplement() const {
  if (primaryByComplementInitialized_) {
    return primaryByComplement_;
  }

  // Complemented flop outputs do not have independent transition equations;
  // their next value is tied to the primary flop. Cache the reverse lookup so
  // each PDR target expansion does not rescan all complemented pairs.
  primaryByComplement_.reserve(
      problem_.complementedStatePairs0.size() +
      problem_.complementedStatePairs1.size());
  for (const auto& [primarySymbol, complementedSymbol] :
       problem_.complementedStatePairs0) {
    primaryByComplement_.emplace(complementedSymbol, primarySymbol);
  }
  for (const auto& [primarySymbol, complementedSymbol] :
       problem_.complementedStatePairs1) {
    primaryByComplement_.emplace(complementedSymbol, primarySymbol);
  }
  primaryByComplementInitialized_ = true;
  return primaryByComplement_;
}

}  // namespace KEPLER_FORMAL::SEC
