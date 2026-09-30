// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "LatchEventModel.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

bool binary(const Bits& bits) {
  return std::all_of(bits.begin(), bits.end(), [](uint8_t bit) { return bit <= 1; });
}

Bits gather(const Bits& bits, const std::vector<size_t>& indices) {
  Bits result;
  result.reserve(indices.size());
  for (size_t index : indices) result.push_back(bits[index]);
  return result;
}

State failure(State state, std::string reason) {
  state.error = true;
  state.errorReason = std::move(reason);
  return state;
}

void appendSize(std::string& key, size_t value) {
  const auto wide = static_cast<uint64_t>(value);
  for (size_t i = 0; i < sizeof(wide); ++i) {
    key.push_back(static_cast<char>((wide >> (i * 8)) & 0xff));
  }
}

void appendBits(std::string& key, const Bits& bits) {
  appendSize(key, bits.size());
  for (auto bit : bits) key.push_back(static_cast<char>(bit));
}

auto reactionTuple(const Reaction& reaction) {
  return std::tie(reaction.storage, reaction.outputs, reaction.error, reaction.reason);
}

}  // namespace

bool State::operator==(const State& other) const {
  return std::tie(current, previous, storage, active, bootstrap, error, errorReason) ==
         std::tie(other.current, other.previous, other.storage, other.active,
                  other.bootstrap, other.error, other.errorReason);
}

bool State::operator<(const State& other) const {
  return std::tie(current, previous, storage, active, bootstrap, error, errorReason) <
         std::tie(other.current, other.previous, other.storage, other.active,
                  other.bootstrap, other.error, other.errorReason);
}

std::string State::key() const {
  std::string result;
  appendBits(result, current);
  appendBits(result, previous);
  appendSize(result, storage.size());
  for (const auto& bits : storage) appendBits(result, bits);
  appendBits(result, active);
  result.push_back(static_cast<char>(bootstrap));
  result.push_back(static_cast<char>(error));
  appendSize(result, errorReason.size());
  result += errorReason;
  return result;
}

EventModel::EventModel(Network network, Limits limits, size_t workerCount)
    : network_(std::move(network)), limits_(limits), workerCount_(workerCount),
      consumers_(network_.netCount) {
  if (!limits_.maxPinOrderings || !limits_.maxSuccessors) {
    throw std::invalid_argument("Latch event resource limits must be positive");
  }
  if (workerCount_ > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("Latch event worker count is too large");
  }
  if (network_.constantByNet.empty()) network_.constantByNet.resize(network_.netCount);
  if (network_.constantByNet.size() != network_.netCount) {
    throw std::invalid_argument("Latch event constant vector has the wrong width");
  }
  Bits driven(network_.netCount, 0);
  const auto writer = [&](size_t net) {
    if (net >= network_.netCount) throw std::invalid_argument("Latch event net out of range");
    if (driven[net]) throw std::invalid_argument("Latch event net has multiple drivers");
    driven[net] = 1;
  };
  for (size_t net : network_.externalInputs) writer(net);
  for (size_t net = 0; net < network_.netCount; ++net) {
    if (network_.constantByNet[net].has_value()) writer(net);
  }
  for (size_t i = 0; i < network_.primitives.size(); ++i) {
    const auto& primitive = network_.primitives[i];
    for (size_t net : primitive.outputs) writer(net);
    for (size_t net : primitive.inputs) {
      if (net >= network_.netCount) throw std::invalid_argument("Latch event input out of range");
      consumers_[net].push_back(i);
    }
  }
  if (std::find(driven.begin(), driven.end(), 0) != driven.end()) {
    throw std::invalid_argument("Latch event net has no driver");
  }
  for (auto& fanout : consumers_) {
    std::sort(fanout.begin(), fanout.end());
    fanout.erase(std::unique(fanout.begin(), fanout.end()), fanout.end());
  }
}

void EventModel::validateState(const State& state) const {
  if (state.current.size() != network_.netCount ||
      state.previous.size() != network_.netCount || !binary(state.current) ||
      !binary(state.previous) || state.storage.size() != network_.primitives.size() ||
      state.active.size() != network_.primitives.size() || !binary(state.active)) {
    throw std::invalid_argument("Malformed Boolean latch event state");
  }
  for (size_t i = 0; i < state.storage.size(); ++i) {
    if (state.storage[i].size() != network_.primitives[i].storageBits ||
        !binary(state.storage[i])) {
      throw std::invalid_argument("Malformed latch event storage");
    }
  }
  for (size_t i = 0; i < network_.netCount; ++i) {
    if (network_.constantByNet[i] &&
        (state.current[i] != *network_.constantByNet[i] ||
         state.previous[i] != *network_.constantByNet[i])) {
      throw std::invalid_argument("Latch event state changed a constant net");
    }
  }
}

void EventModel::normalizeBoundary(State& state) const {
  if (!state.bootstrap && !state.error &&
      std::none_of(state.active.begin(), state.active.end(), [](auto bit) { return bit; })) {
    state.previous = state.current;
  }
}

bool EventModel::stable(const State& state) const {
  validateState(state);
  return !state.bootstrap && !state.error && state.previous == state.current &&
         std::none_of(state.active.begin(), state.active.end(), [](auto bit) { return bit; });
}

State EventModel::bootstrap(const Bits& inputs, const std::vector<Bits>& initialStorage,
                            const Bits& internalSeeds) const {
  if (inputs.size() != network_.externalInputs.size() || !binary(inputs) ||
      internalSeeds.size() != network_.netCount || !binary(internalSeeds) ||
      initialStorage.size() != network_.primitives.size()) {
    throw std::invalid_argument("Malformed latch event bootstrap parameters");
  }
  State state;
  state.current = internalSeeds;
  state.storage = initialStorage;
  state.active.assign(network_.primitives.size(), 1);
  state.bootstrap = true;
  for (size_t i = 0; i < inputs.size(); ++i) state.current[network_.externalInputs[i]] = inputs[i];
  for (size_t i = 0; i < network_.netCount; ++i) {
    if (network_.constantByNet[i]) state.current[i] = *network_.constantByNet[i];
  }
  const Bits projectionSnapshot = state.current;
  for (size_t i = 0; i < network_.primitives.size(); ++i) {
    const auto& primitive = network_.primitives[i];
    const auto& storage = initialStorage[i];
    if (storage.size() != primitive.storageBits || !binary(storage)) {
      throw std::invalid_argument("Malformed latch event bootstrap storage");
    }
    if (!primitive.storageBits) continue;
    Bits outputs;
    try {
      if (primitive.initialOutputValues) {
        outputs = primitive.initialOutputValues(storage, gather(projectionSnapshot, primitive.inputs));
      } else if (primitive.initialOutputs) outputs = primitive.initialOutputs(storage);
      else if (primitive.outputs.size() == storage.size()) outputs = storage;
      else {
        state = failure(std::move(state), primitive.name + ": missing initial output mapping");
        continue;
      }
    } catch (const Limit&) {
      throw;
    } catch (const std::bad_alloc&) {
      throw;
    } catch (const std::exception& exception) {
      state = failure(std::move(state), primitive.name + ": initial output mapping: " + exception.what());
      continue;
    }
    if (outputs.size() != primitive.outputs.size() || !binary(outputs)) {
      state = failure(std::move(state), primitive.name + ": invalid initial output mapping");
      continue;
    }
    for (size_t j = 0; j < outputs.size(); ++j) state.current[primitive.outputs[j]] = outputs[j];
  }
  state.previous = state.current;  // In particular, no invented external clock edge.
  return state;
}

State EventModel::admit(const State& quiescent, const Bits& inputs) const {
  validateState(quiescent);
  if (quiescent.error) return quiescent;
  if (!stable(quiescent)) return failure(quiescent, "External transaction before quiescence");
  if (inputs.size() != network_.externalInputs.size() || !binary(inputs)) {
    return failure(quiescent, "Invalid external Boolean transaction");
  }
  State state = quiescent;
  state.previous = quiescent.current;
  state.active.assign(network_.primitives.size(), 0);
  for (size_t i = 0; i < inputs.size(); ++i) {
    const size_t net = network_.externalInputs[i];
    state.current[net] = inputs[i];
    if (inputs[i] != state.previous[net]) {
      for (size_t consumer : consumers_[net]) state.active[consumer] = 1;
    }
  }
  normalizeBoundary(state);
  return state;
}

std::vector<Reaction> EventModel::evaluate(size_t index, const State& state) const {
  const auto& primitive = network_.primitives[index];
  const auto oldOutputs = gather(state.current, primitive.outputs);
  const auto oldStorage = state.storage[index];
  if (!state.active[index]) return {{oldStorage, oldOutputs}};
  const auto current = gather(state.current, primitive.inputs);
  const auto previous = gather(state.previous, primitive.inputs);
  const auto invoke = [&](const Bits& storage, const Bits& before, const Bits& pins,
                          std::optional<size_t> changedPin) -> Reaction {
    auto invalid = [&](const std::string& reason) {
      return Reaction{storage, oldOutputs, true, primitive.name + ": " + reason};
    };
    if (!primitive.react) return invalid("missing primitive reaction");
    Reaction result;
    try {
      result = primitive.react(storage, before, pins, changedPin, state.bootstrap);
    } catch (const Limit&) {
      throw;
    } catch (const std::bad_alloc&) {
      throw;
    } catch (const std::exception& exception) {
      return invalid(std::string("primitive reaction: ") + exception.what());
    }
    if (result.error) return invalid(result.reason.empty() ? "unsupported primitive reaction" : result.reason);
    if (result.storage.size() != primitive.storageBits || !binary(result.storage) ||
        result.outputs.size() != primitive.outputs.size() || !binary(result.outputs)) {
      return invalid("invalid primitive reaction shape or Boolean value");
    }
    return result;
  };
  if (state.bootstrap || primitive.storageBits == 0) {
    return {invoke(oldStorage, previous, current, std::nullopt)};
  }
  std::vector<size_t> changed;
  for (size_t pin = 0; pin < current.size(); ++pin) {
    if (current[pin] != previous[pin]) changed.push_back(pin);
  }
  size_t orderingCount = 1;
  for (size_t i = 2; i <= changed.size(); ++i) {
    if (orderingCount > limits_.maxPinOrderings / i) {
      throw Limit(primitive.name + ": changed-pin ordering limit exceeded");
    }
    orderingCount *= i;
  }
  if (changed.empty()) return {invoke(oldStorage, previous, current, std::nullopt)};
  std::vector<Reaction> results;
  do {
    Bits pins = previous;
    Reaction result{oldStorage, oldOutputs};
    for (size_t pin : changed) {
      const Bits before = pins;
      pins[pin] = current[pin];
      result = invoke(result.storage, before, pins, pin);
      if (result.error) break;
    }
    results.push_back(std::move(result));
  } while (std::next_permutation(changed.begin(), changed.end()));
  std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) {
    return reactionTuple(a) < reactionTuple(b);
  });
  results.erase(std::unique(results.begin(), results.end(), [](const auto& a, const auto& b) {
    return reactionTuple(a) == reactionTuple(b);
  }), results.end());
  return results;
}

std::vector<State> EventModel::successors(const State& state) const {
  validateState(state);
  if (state.error || stable(state)) return {state};
  std::vector<std::vector<Reaction>> choices(network_.primitives.size());
  const auto evaluateAll = [&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, choices.size()), [&](const auto& range) {
      for (size_t i = range.begin(); i != range.end(); ++i) choices[i] = evaluate(i, state);
    });
  };
  tbb::task_arena arena(workerCount_ ? static_cast<int>(workerCount_) : tbb::task_arena::automatic);
  arena.execute(evaluateAll);

  size_t count = 1;
  for (const auto& alternatives : choices) {
    if (count > limits_.maxSuccessors / alternatives.size()) {
      throw Limit("Latch event successor limit exceeded");
    }
    count *= alternatives.size();
  }
  State base = state;
  base.bootstrap = false;
  base.previous = state.current;
  base.active.assign(network_.primitives.size(), 0);
  std::vector<State> results{std::move(base)};
  for (size_t i = 0; i < choices.size(); ++i) {
    std::vector<State> expanded;
    expanded.reserve(results.size() * choices[i].size());
    for (const auto& partial : results) {
      for (const auto& reaction : choices[i]) {
        State next = partial;
        next.storage[i] = reaction.storage;
        for (size_t j = 0; j < reaction.outputs.size(); ++j) {
          next.current[network_.primitives[i].outputs[j]] = reaction.outputs[j];
        }
        if (reaction.error && !next.error) {
          next.error = true;
          next.errorReason = reaction.reason;
        }
        expanded.push_back(std::move(next));
      }
    }
    results = std::move(expanded);
  }
  for (auto& next : results) {
    for (size_t net = 0; net < network_.netCount; ++net) {
      if (next.current[net] != state.current[net]) {
        for (size_t consumer : consumers_[net]) next.active[consumer] = 1;
      }
    }
    normalizeBoundary(next);
  }
  std::sort(results.begin(), results.end());
  results.erase(std::unique(results.begin(), results.end()), results.end());
  return results;
}

Primitive combinational(std::string name, std::vector<size_t> inputs,
                        std::vector<size_t> outputs,
                        std::function<Bits(const Bits&)> function) {
  Primitive primitive;
  primitive.name = std::move(name);
  primitive.inputs = std::move(inputs);
  primitive.outputs = std::move(outputs);
  if (function) {
    primitive.react = [function = std::move(function)](const Bits&, const Bits&, const Bits& pins,
                                                      std::optional<size_t>, bool) {
      return Reaction{{}, function(pins)};
    };
  }
  return primitive;
}

Primitive latch(std::string name, size_t data, size_t enable, size_t output, bool activeHigh) {
  Primitive primitive;
  primitive.name = std::move(name);
  primitive.inputs = {data, enable};
  primitive.outputs = {output};
  primitive.storageBits = 1;
  primitive.react = [activeHigh](const Bits& storage, const Bits&, const Bits& pins,
                                std::optional<size_t>, bool) {
    const Bits value{pins[1] == activeHigh ? pins[0] : storage[0]};
    return Reaction{value, value};
  };
  return primitive;
}

Primitive flipFlop(std::string name, size_t data, size_t clock, size_t output, bool risingEdge) {
  Primitive primitive;
  primitive.name = std::move(name);
  primitive.inputs = {data, clock};
  primitive.outputs = {output};
  primitive.storageBits = 1;
  primitive.react = [risingEdge](const Bits& storage, const Bits& before, const Bits& pins,
                                std::optional<size_t> changedPin, bool bootstrap) {
    const bool edge = !bootstrap && changedPin == std::optional<size_t>{1} &&
                      before[1] != pins[1] && pins[1] == risingEdge;
    const Bits value{edge ? pins[0] : storage[0]};
    return Reaction{value, value};
  };
  return primitive;
}

}  // namespace KEPLER_FORMAL::SEC::LATCH
