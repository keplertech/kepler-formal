// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "latch/LatchResetAdapter.h"

#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using Expressions = std::unordered_map<size_t, BoolExpr*>;

// A cumulative reachable-DAG work budget, not a process peak-memory limit.
class CompositionBudget {
 public:
  explicit CompositionBudget(size_t limit) : limit_(limit) {
    if (!limit) throw std::invalid_argument("reset composition node budget must be positive");
  }
  void account(BoolExpr* root) {
    std::vector<BoolExpr*> pending{root};
    while (!pending.empty()) {
      auto* node = pending.back();
      pending.pop_back();
      if (!node || !seen_.insert(node).second) continue;
      if (seen_.size() > limit_)
        throw std::runtime_error("reset composition node budget exceeded");
      if (node->getLeft()) pending.push_back(node->getLeft());
      if (node->getRight()) pending.push_back(node->getRight());
    }
  }
 private:
  size_t limit_;
  std::unordered_set<BoolExpr*> seen_;
};

BoolExpr* constant(bool value) { return value ? BoolExpr::createTrue() : BoolExpr::createFalse(); }
BoolExpr* mux(BoolExpr* select, BoolExpr* yes, BoolExpr* no) {
  if (yes == no) return yes;
  return BoolExpr::Or(BoolExpr::And(select, yes), BoolExpr::And(BoolExpr::Not(select), no));
}

// Simultaneous substitution: replacement expressions are not themselves
// substituted. One memo is shared across every output and state in a stage.
class Substitute {
 public:
  Substitute(const Expressions& replacements, CompositionBudget& budget)
      : replacements_(replacements), budget_(budget) {}
  BoolExpr* operator()(BoolExpr* root) {
    std::vector<std::pair<BoolExpr*, bool>> pending{{root, false}};
    while (!pending.empty()) {
      const auto [node, ready] = pending.back();
      pending.pop_back();
      if (!node || !node->isValid()) throw std::invalid_argument("invalid event expression");
      if (memo_.count(node)) continue;
      if (node->getOp() == Op::VAR) {
        const auto found = replacements_.find(node->getId());
        memo_[node] = found == replacements_.end() ? node : found->second;
      } else if (!ready) {
        pending.emplace_back(node, true);
        if (node->getRight()) pending.emplace_back(node->getRight(), false);
        pending.emplace_back(node->getLeft(), false);
      } else {
        auto* left = memo_.at(node->getLeft());
        switch (node->getOp()) {
          case Op::NOT: memo_[node] = BoolExpr::Not(left); break;
          case Op::AND: memo_[node] = BoolExpr::And(left, memo_.at(node->getRight())); break;
          case Op::OR: memo_[node] = BoolExpr::Or(left, memo_.at(node->getRight())); break;
          case Op::XOR: memo_[node] = BoolExpr::Xor(left, memo_.at(node->getRight())); break;
          default: throw std::invalid_argument("unsupported event expression operation");
        }
      }
      if (memo_.count(node)) budget_.account(memo_.at(node));
    }
    return memo_.at(root);
  }
 private:
  const Expressions& replacements_;
  CompositionBudget& budget_;
  std::unordered_map<BoolExpr*, BoolExpr*> memo_;
};

class Composer {
 public:
  Composer(const SequentialDesignModel& model, const EventResetInterface& interface,
           CompositionBudget& budget)
      : model_(model), interface_(interface), budget_(budget) {}

  Expressions initial() const {
    Expressions state;
    for (const auto& key : model_.stateBits) {
      const auto id = model_.inputVarByKey.at(key);
      state.emplace(id, BoolExpr::Var(id));
    }
    return state;
  }

  Expressions force(const Expressions& state, size_t index, bool value) const {
    Expressions environment;
    if (interface_.singleInputChange) {
      select(environment, index, constant(value));
    } else {
      Substitute previous(state, budget_);
      for (size_t i = 0; i < interface_.inputKeys.size(); ++i)
        environment.emplace(model_.inputVarByKey.at(interface_.inputKeys[i]),
            i == index ? constant(value) : previous(interface_.currentInputs[i]));
    }
    return transition(state, environment);
  }

  Expressions external(const Expressions& state, size_t resetIndex, bool active,
                       std::optional<size_t> heldClock) const {
    return transition(state, externalEnvironment(resetIndex, active, heldClock));
  }

  Expressions sampleAll(Expressions state, const std::vector<size_t>& pins,
                        const std::vector<std::vector<BoolExpr*>>& ordering) const {
    std::vector<BoolExpr*> used(pins.size(), BoolExpr::createFalse());
    for (size_t stage = 0; stage < pins.size(); ++stage) {
      std::vector<BoolExpr*> requested;
      auto* valid = BoolExpr::createFalse();
      for (size_t candidate = 0; candidate < pins.size(); ++candidate) {
        auto* match = BoolExpr::Not(used[candidate]);
        for (size_t bit = 0; bit < ordering[stage].size(); ++bit)
          match = BoolExpr::And(match, ((candidate >> bit) & 1) ? ordering[stage][bit]
              : BoolExpr::Not(ordering[stage][bit]));
        requested.push_back(match);
        valid = BoolExpr::Or(valid, match);
      }
      std::vector<BoolExpr*> selected;
      auto* earlierUsed = BoolExpr::createTrue();
      for (size_t candidate = 0; candidate < pins.size(); ++candidate) {
        auto* fallback = BoolExpr::And(earlierUsed, BoolExpr::Not(used[candidate]));
        selected.push_back(BoolExpr::Or(requested[candidate],
            BoolExpr::And(BoolExpr::Not(valid), fallback)));
        earlierUsed = BoolExpr::And(earlierUsed, used[candidate]);
      }
      Expressions environment;
      for (size_t bit = 0; bit < interface_.selectorSymbols.size(); ++bit) {
        auto* value = BoolExpr::createFalse();
        for (size_t candidate = 0; candidate < pins.size(); ++candidate)
          if ((pins[candidate] >> bit) & 1) value = BoolExpr::Or(value, selected[candidate]);
        environment.emplace(interface_.selectorSymbols[bit], value);
      }
      auto* value = BoolExpr::createFalse();
      for (size_t candidate = 0; candidate < pins.size(); ++candidate) {
        value = BoolExpr::Or(value, BoolExpr::And(selected[candidate],
            BoolExpr::Var(model_.inputVarByKey.at(interface_.inputKeys[pins[candidate]]))));
        used[candidate] = BoolExpr::Or(used[candidate], selected[candidate]);
      }
      environment.emplace(*interface_.valueSymbol, value);
      state = transition(state, environment);
    }
    return state;
  }

  Expressions externalEnvironment(size_t resetIndex, bool active,
                                  std::optional<size_t> heldClock) const {
    Expressions environment;
    if (interface_.singleInputChange) {
      auto* blocked = selected(resetIndex);
      if (heldClock) blocked = BoolExpr::Or(blocked, selected(*heldClock));
      const size_t stutter = interface_.inputKeys.size();
      for (size_t bit = 0; bit < interface_.selectorSymbols.size(); ++bit) {
        const auto id = interface_.selectorSymbols[bit];
        environment.emplace(id, mux(blocked, constant((stutter >> bit) & 1), BoolExpr::Var(id)));
      }
      environment.emplace(*interface_.valueSymbol, BoolExpr::Var(*interface_.valueSymbol));
    } else {
      for (size_t index = 0; index < interface_.inputKeys.size(); ++index) {
        const auto id = model_.inputVarByKey.at(interface_.inputKeys[index]);
        environment.emplace(id, index == resetIndex ? constant(active) :
            heldClock && index == *heldClock ? constant(false) : BoolExpr::Var(id));
      }
    }
    return environment;
  }

 private:
  BoolExpr* selected(size_t index) const {
    auto* result = BoolExpr::createTrue();
    for (size_t bit = 0; bit < interface_.selectorSymbols.size(); ++bit) {
      auto* variable = BoolExpr::Var(interface_.selectorSymbols[bit]);
      result = BoolExpr::And(result, ((index >> bit) & 1) ? variable : BoolExpr::Not(variable));
    }
    return result;
  }
  void select(Expressions& environment, size_t index, BoolExpr* value) const {
    for (size_t bit = 0; bit < interface_.selectorSymbols.size(); ++bit)
      environment.emplace(interface_.selectorSymbols[bit], constant((index >> bit) & 1));
    environment.emplace(*interface_.valueSymbol, value);
  }
  Expressions transition(const Expressions& state, const Expressions& environment) const {
    Expressions replacements = state;
    for (const auto& [id, expression] : environment) {
      budget_.account(expression);
      replacements.insert_or_assign(id, expression);
    }
    Substitute substitute(replacements, budget_);
    Expressions next;
    for (const auto& key : model_.stateBits)
      next.emplace(model_.inputVarByKey.at(key), substitute(model_.nextStateExprByStateKey.at(key)));
    return next;
  }
  const SequentialDesignModel& model_;
  const EventResetInterface& interface_;
  CompositionBudget& budget_;
};

std::string validate(const SequentialDesignModel& model, const SecResetSpec& reset,
                     size_t& resetIndex) {
  if (model.eventContract.empty() || !model.eventResetInterface)
    return "Latch reset cycles require retained event reset interface metadata";
  if (reset.cycles == 0) return "Latch reset cycles must be positive";
  if (reset.cycles == std::numeric_limits<size_t>::max())
    return "Latch reset cycle count is too large";
  if (reset.ports.size() != 1)
    return "Latch reset cycles currently require exactly one reset port";
  const auto& interface = *model.eventResetInterface;
  if (!interface.clockError.empty()) return interface.clockError;
  const size_t inputs = interface.inputKeys.size();
  if (!interface.clockInputIndex || *interface.clockInputIndex >= inputs)
    return "Latch reset cycles require one unambiguous external clock carrier";
  if (interface.inputNames.size() != inputs || interface.currentInputs.size() != inputs)
    return "Latch reset interface has inconsistent input metadata";
  const auto& requested = reset.ports.front().name;
  std::vector<size_t> exact, bits;
  for (size_t i = 0; i < inputs; ++i) {
    const auto& name = interface.inputNames[i];
    if (name == requested) exact.push_back(i);
    else if (requested.find('[') == std::string::npos && name.size() > requested.size() + 2 &&
        name.compare(0, requested.size(), requested) == 0 && name[requested.size()] == '[' && name.back() == ']')
      bits.push_back(i);
  }
  if (exact.size() == 1) resetIndex = exact.front();
  else if (exact.empty() && bits.size() == 1 && interface.inputNames[bits.front()] == requested + "[0]")
    resetIndex = bits.front();
  else return "Latch reset port `" + requested + "` is missing or ambiguous; name a bus bit explicitly";
  if (resetIndex == *interface.clockInputIndex)
    return "Latch reset port cannot also be the cycle clock carrier";
  std::set<size_t> stateSymbols, environmentSymbols;
  for (const auto& key : model.stateBits) {
    const auto found = model.inputVarByKey.find(key);
    if (found == model.inputVarByKey.end() || found->second < 2 ||
        !stateSymbols.insert(found->second).second ||
        !model.nextStateExprByStateKey.count(key) || !model.initialStateValueByKey.count(key))
      return "Latch reset cycles require fully initialized event state";
  }
  for (const auto& key : model.environmentInputs) {
    const auto found = model.inputVarByKey.find(key);
    if (found == model.inputVarByKey.end() || found->second < 2 ||
        stateSymbols.count(found->second) || !environmentSymbols.insert(found->second).second)
      return "Latch reset interface has invalid environment symbols";
  }
  std::set<size_t> inputSymbols;
  for (size_t i = 0; i < inputs; ++i) {
    if (!interface.currentInputs[i] || !interface.currentInputs[i]->isValid())
      return "Latch reset interface lacks remembered external input levels";
    for (auto id : interface.currentInputs[i]->getSupportVars())
      if (id >= 2 && !stateSymbols.count(id))
        return "Latch reset remembered input levels must depend only on state";
    if (!model.inputVarByKey.count(interface.inputKeys[i]) ||
        !environmentSymbols.count(model.inputVarByKey.at(interface.inputKeys[i])) ||
        !inputSymbols.insert(model.inputVarByKey.at(interface.inputKeys[i])).second)
      return "Latch reset interface input is not an environment variable";
  }
  if (interface.singleInputChange) {
    if (!interface.valueSymbol || !environmentSymbols.count(*interface.valueSymbol) ||
        interface.selectorSymbols.empty() ||
        interface.selectorSymbols.size() >= std::numeric_limits<size_t>::digits ||
        (size_t(1) << interface.selectorSymbols.size()) <= inputs)
      return "Latch reset interface has invalid selector encoding";
    std::set<size_t> selectors{*interface.valueSymbol};
    for (auto id : interface.selectorSymbols)
      if (!environmentSymbols.count(id) || !selectors.insert(id).second)
        return "Latch reset interface has invalid selector symbols";
  }
  return {};
}
}  // namespace

ResetCycleAdaptation adaptResetCycles(const SequentialDesignModel& model, const SecResetSpec& reset) {
  size_t resetIndex = 0;
  if (auto error = validate(model, reset, resetIndex); !error.empty()) return {{}, std::move(error)};
  try {
    const auto& interface = *model.eventResetInterface;
    const size_t clock = *interface.clockInputIndex;
    const bool asserted = reset.ports.front().activeValue;
    CompositionBudget budget(interface.maxCompositionNodes);
    for (const auto& [key, expression] : model.nextStateExprByStateKey) budget.account(expression);
    for (const auto& [key, expression] : model.observedOutputExprByKey) budget.account(expression);
    for (auto* expression : interface.currentInputs) budget.account(expression);
    SequentialDesignModel adapted = model;
    size_t nextId = 2;
    for (const auto& [key, id] : model.inputVarByKey) {
      if (id >= std::numeric_limits<size_t>::max() - std::numeric_limits<size_t>::digits - 1)
        return {{}, "Latch reset adapter has no free state symbol range"};
      nextId = std::max(nextId, id + 1);
    }
    std::vector<size_t> sampledPins;
    std::vector<std::vector<BoolExpr*>> ordering;
    if (interface.singleInputChange) {
      for (size_t i = 0; i < interface.inputKeys.size(); ++i)
        if (i != resetIndex && i != clock) sampledPins.push_back(i);
      if (!sampledPins.empty() && sampledPins.size() > interface.maxCompositionNodes / sampledPins.size())
        return {{}, "Latch reset ordering exceeds the composition node budget"};
      size_t width = 0;
      for (size_t count = sampledPins.empty() ? 0 : sampledPins.size() - 1; count; count >>= 1) ++width;
      ordering.resize(sampledPins.size());
      for (size_t stage = 0; stage < sampledPins.size(); ++stage) {
        for (size_t bit = 0; bit < width; ++bit) {
          if (nextId >= std::numeric_limits<size_t>::max() - std::numeric_limits<size_t>::digits - 1)
            return {{}, "Latch reset adapter has no free ordering symbol range"};
          SignalKey key{{uint64_t(1) << 61, 5, stage, bit}, {0}};
          if (adapted.inputVarByKey.count(key)) return {{}, "Latch reset ordering key collides with existing input"};
          ordering[stage].push_back(BoolExpr::Var(nextId));
          budget.account(ordering[stage].back());
          adapted.environmentInputs.push_back(key);
          adapted.inputVarByKey.emplace(key, nextId++);
          adapted.displayNameByKey.emplace(key, "$event.reset.order[" + std::to_string(stage) +
              "][" + std::to_string(bit) + "]");
        }
      }
    }
    Composer composer(model, interface, budget);
    const auto original = composer.initial();
    auto boot = composer.force(original, resetIndex, asserted);
    boot = composer.force(boot, clock, false);
    boot = interface.singleInputChange ? composer.sampleAll(std::move(boot), sampledPins, ordering)
        : composer.external(boot, resetIndex, asserted, clock);
    boot = composer.force(boot, clock, true);
    boot = composer.force(boot, clock, false);
    const auto released = composer.force(boot, resetIndex, !asserted);
    const auto normal = composer.external(original, resetIndex, !asserted, {});
    std::vector<BoolExpr*> count;
    std::vector<SignalKey> countKeys;
    for (size_t value = reset.cycles, bit = 0; value != 0; value >>= 1, ++bit) {
      SignalKey key{{uint64_t(1) << 61, 4, bit}, {0}};
      if (adapted.inputVarByKey.count(key)) return {{}, "Latch reset counter key collides with existing state"};
      countKeys.push_back(key);
      count.push_back(BoolExpr::Var(nextId));
      adapted.inputVarByKey.emplace(key, nextId++);
      adapted.stateBits.push_back(key);
      adapted.initialStateValueByKey.emplace(key, (reset.cycles >> bit) & 1);
      adapted.displayNameByKey.emplace(key, "$event.reset.remaining[" + std::to_string(bit) + "]");
    }
    auto* active = BoolExpr::createFalse();
    auto* last = count.front();
    for (size_t bit = 0; bit < count.size(); ++bit) {
      active = BoolExpr::Or(active, count[bit]);
      if (bit) last = BoolExpr::And(last, BoolExpr::Not(count[bit]));
    }
    auto* borrow = BoolExpr::createTrue();
    for (size_t bit = 0; bit < count.size(); ++bit) {
      adapted.nextStateExprByStateKey.emplace(countKeys[bit],
          BoolExpr::And(active, BoolExpr::Xor(count[bit], borrow)));
      budget.account(adapted.nextStateExprByStateKey.at(countKeys[bit]));
      borrow = BoolExpr::And(borrow, BoolExpr::Not(count[bit]));
    }
    for (const auto& key : model.stateBits) {
      const auto id = model.inputVarByKey.at(key);
      adapted.nextStateExprByStateKey[key] = mux(active,
          mux(last, released.at(id), boot.at(id)), normal.at(id));
      budget.account(adapted.nextStateExprByStateKey.at(key));
    }
    const auto environment = composer.externalEnvironment(resetIndex, !asserted, {});
    Substitute observations(environment, budget);
    for (const auto& key : model.observedOutputs) {
      adapted.observedOutputExprByKey[key] = BoolExpr::And(BoolExpr::Not(active),
          observations(model.observedOutputExprByKey.at(key)));
      budget.account(adapted.observedOutputExprByKey.at(key));
    }
    adapted.eventContract += ";reset_cycles=" + std::to_string(reset.cycles) +
        ";reset_port=" + interface.inputNames[resetIndex] +
        ";reset_active=" + std::to_string(asserted) +
        ";reset_clock=" + interface.inputNames[clock] + ";reset_protocol=low-sample-high-low-release" +
        (interface.singleInputChange ? ";reset_sampling=all-levels-all-arrival-orders" : ";reset_sampling=atomic-level-vector");
    adapted.eventResetCycles = reset.cycles;
    // A second adaptation must not accidentally apply another prefix.
    adapted.eventResetInterface.reset();
    return {std::move(adapted), {}};
  } catch (const std::exception& error) {
    return {{}, std::string("Cannot compose latch reset cycles: ") + error.what()};
  }
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
