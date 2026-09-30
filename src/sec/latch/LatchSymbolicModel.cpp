// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchSymbolicModel.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
BoolExpr* zero() { return BoolExpr::createFalse(); }
BoolExpr* one() { return BoolExpr::createTrue(); }
BoolExpr* negate(BoolExpr* value) { return BoolExpr::Not(value); }
BoolExpr* conjunction(BoolExpr* a, BoolExpr* b) { return BoolExpr::And(a, b); }
BoolExpr* disjunction(BoolExpr* a, BoolExpr* b) { return BoolExpr::Or(a, b); }
BoolExpr* different(BoolExpr* a, BoolExpr* b) { return BoolExpr::Xor(a, b); }

void validate(const SymbolicBits& bits, size_t size) {
  if (bits.size() != size || std::any_of(bits.begin(), bits.end(), [](auto* bit) {
        return !bit || !bit->isValid();
      })) throw std::invalid_argument("Malformed symbolic Boolean vector");
}
void validate(const Network& network, const SymbolicState& state) {
  validate(state.current, network.netCount);
  validate(state.previous, network.netCount);
  validate(state.active, network.primitives.size());
  validate({state.bootstrap, state.error}, 2);
  if (state.storage.size() != network.primitives.size())
    throw std::invalid_argument("Malformed symbolic primitive storage");
  for (size_t i = 0; i < state.storage.size(); ++i)
    validate(state.storage[i], network.primitives[i].storageBits);
}
SymbolicBits gather(const SymbolicBits& bits, const std::vector<size_t>& indices) {
  SymbolicBits result;
  for (auto index : indices) result.push_back(bits.at(index));
  return result;
}
SymbolicBits choose(BoolExpr* condition, const SymbolicBits& yes, const SymbolicBits& no) {
  if (yes.size() != no.size()) throw std::invalid_argument("Symbolic mux width mismatch");
  SymbolicBits result;
  for (size_t i = 0; i < yes.size(); ++i) result.push_back(symbolicMux(condition, yes[i], no[i]));
  return result;
}
BoolExpr* equal(const SymbolicBits& a, const SymbolicBits& b) {
  if (a.size() != b.size()) throw std::invalid_argument("Symbolic equality width mismatch");
  auto* result = one();
  for (size_t i = 0; i < a.size(); ++i) result = conjunction(result, negate(different(a[i], b[i])));
  return result;
}
BoolExpr* inactive(const SymbolicState& state) {
  auto* result = conjunction(negate(state.bootstrap), negate(state.error));
  for (auto* active : state.active) result = conjunction(result, negate(active));
  return result;
}
void normalizeCompletedActivation(SymbolicState& state,
    const std::vector<std::vector<size_t>>& consumers) {
  // Only call after admission/wave computed the COMPLETE fanout of changes
  // against previous. For a consumed net, inactivity already implies that net
  // did not change, so its history needs no mux. Sink nets still need global
  // normalization. This keeps a global guard out of every primitive's inputs.
  auto* settled = inactive(state);
  for (size_t net = 0; net < state.current.size(); ++net)
    if (consumers[net].empty())
      state.previous[net] = symbolicMux(settled, state.current[net], state.previous[net]);
}
SymbolicState choose(BoolExpr* condition, const SymbolicState& yes, const SymbolicState& no) {
  SymbolicState result;
  result.current = choose(condition, yes.current, no.current);
  result.previous = choose(condition, yes.previous, no.previous);
  result.active = choose(condition, yes.active, no.active);
  for (size_t i = 0; i < yes.storage.size(); ++i)
    result.storage.push_back(choose(condition, yes.storage[i], no.storage[i]));
  result.bootstrap = symbolicMux(condition, yes.bootstrap, no.bootstrap);
  result.error = symbolicMux(condition, yes.error, no.error);
  return result;
}
SymbolicReaction choose(BoolExpr* condition, const SymbolicReaction& yes, const SymbolicReaction& no) {
  return {choose(condition, yes.storage, no.storage), choose(condition, yes.outputs, no.outputs),
          symbolicMux(condition, yes.error, no.error)};
}

SymbolicReaction invoke(const Primitive& primitive, const SymbolicPrimitive& callback,
    const SymbolicBits& storage, const SymbolicBits& before, const SymbolicBits& current,
    const SymbolicBits& oldOutputs, std::optional<size_t> changed, bool bootstrap) {
  SymbolicReaction result{storage, oldOutputs, one()};
  if (!callback.react) return result;
  try {
    result = callback.react(storage, before, current, changed, bootstrap);
    validate(result.storage, primitive.storageBits);
    validate(result.outputs, primitive.outputs.size());
    validate({result.error}, 1);
  } catch (const Limit&) { throw;
  } catch (const std::bad_alloc&) { throw;
  } catch (const std::exception&) {
    return {storage, oldOutputs, one()};
  }
  // Concrete invalid/error reactions keep the incoming storage and original
  // physical outputs, even after earlier successful visits in a permutation.
  result.storage = choose(result.error, storage, result.storage);
  result.outputs = choose(result.error, oldOutputs, result.outputs);
  return result;
}

size_t permutationCount(size_t pins, size_t limit) {
  size_t count = 1;
  for (size_t i = 2; i <= pins; ++i) {
    if (count > limit / i) throw Limit("symbolic changed-pin ordering limit exceeded");
    count *= i;
  }
  return count;
}
BoolExpr* code(const SymbolicBits& bits, size_t value) {
  auto* result = one();
  for (size_t i = 0; i < bits.size(); ++i)
    result = conjunction(result, (value >> i) & 1 ? bits[i] : negate(bits[i]));
  return result;
}

bool binary(const Bits& bits) {
  return std::all_of(bits.begin(), bits.end(), [](auto bit) { return bit <= 1; });
}
size_t domainSize(size_t storage, size_t pins, size_t copies, size_t maximum) {
  if (storage > maximum || pins > (maximum - storage) / copies)
    throw Limit("local primitive truth-table lifting limit exceeded");
  const size_t bits = storage + copies * pins;
  if (bits >= std::numeric_limits<size_t>::digits)
    throw Limit("local primitive truth-table lifting index overflow");
  return bits;
}
BoolExpr* minterm(const SymbolicBits& expressions, size_t value) {
  return code(expressions, value);
}
Bits unpack(size_t value, size_t offset, size_t size) {
  Bits result(size);
  for (size_t i = 0; i < size; ++i) result[i] = (value >> (offset + i)) & 1;
  return result;
}
}  // namespace

BoolExpr* symbolicMux(BoolExpr* condition, BoolExpr* yes, BoolExpr* no) {
  if (yes == no) return yes;
  if (condition == one()) return yes;
  if (condition == zero()) return no;
  return disjunction(conjunction(condition, yes), conjunction(negate(condition), no));
}

SymbolicBits flattenSymbolicBoundary(const SymbolicState& state) {
  auto result = state.current;
  for (const auto& storage : state.storage) result.insert(result.end(), storage.begin(), storage.end());
  return result;
}

SymbolicEventModel::SymbolicEventModel(SymbolicNetwork network, Limits limits, size_t workers)
    : network_(std::move(network)), limits_(limits), workers_(workers) {
  // Reuse all concrete structural/single-writer validation, not its callbacks.
  EventModel validated(network_.reference, limits, workers);
  network_.reference = validated.network();
  if (network_.primitives.size() != network_.reference.primitives.size())
    throw std::invalid_argument("Symbolic primitive count mismatch");
  consumers_.resize(network_.reference.netCount);
  for (size_t i = 0; i < network_.reference.primitives.size(); ++i)
    for (auto net : network_.reference.primitives[i].inputs) consumers_[net].push_back(i);
  for (auto& consumers : consumers_) {
    std::sort(consumers.begin(), consumers.end());
    consumers.erase(std::unique(consumers.begin(), consumers.end()), consumers.end());
  }
}

SymbolicState SymbolicEventModel::boundary(const SymbolicBits& current,
    const std::vector<SymbolicBits>& storage) const {
  SymbolicState result;
  result.current = result.previous = current;
  result.storage = storage;
  result.active.assign(network_.primitives.size(), zero());
  validate(network_.reference, result);
  return result;
}

SymbolicState SymbolicEventModel::bootstrap(const SymbolicBits& inputs,
    const std::vector<SymbolicBits>& storage, const SymbolicBits& seeds) const {
  const auto& reference = network_.reference;
  validate(inputs, reference.externalInputs.size());
  auto result = boundary(seeds, storage);
  result.bootstrap = one();
  result.active.assign(network_.primitives.size(), one());
  for (size_t i = 0; i < inputs.size(); ++i) result.current[reference.externalInputs[i]] = inputs[i];
  for (size_t i = 0; i < reference.netCount; ++i)
    if (reference.constantByNet[i]) result.current[i] = *reference.constantByNet[i] ? one() : zero();
  const auto snapshot = result.current;
  for (size_t i = 0; i < reference.primitives.size(); ++i) {
    const auto& primitive = reference.primitives[i];
    if (!primitive.storageBits) continue;
    SymbolicBits outputs;
    try {
      if (network_.primitives[i].initialOutputs)
        outputs = network_.primitives[i].initialOutputs(storage[i], gather(snapshot, primitive.inputs));
      else if (!primitive.initialOutputs && !primitive.initialOutputValues &&
               primitive.outputs.size() == primitive.storageBits) outputs = storage[i];
      else throw std::invalid_argument("Missing symbolic initial output projection");
      validate(outputs, primitive.outputs.size());
    } catch (const Limit&) { throw;
    } catch (const std::bad_alloc&) { throw;
    } catch (const std::exception&) {
      result.error = one();
      continue;
    }
    for (size_t bit = 0; bit < outputs.size(); ++bit) result.current[primitive.outputs[bit]] = outputs[bit];
  }
  result.previous = result.current;
  return result;
}

BoolExpr* SymbolicEventModel::stable(const SymbolicState& state) const {
  validate(network_.reference, state);
  return conjunction(inactive(state), equal(state.current, state.previous));
}

SymbolicState SymbolicEventModel::admit(const SymbolicState& state, const SymbolicBits& inputs) const {
  validate(network_.reference, state);
  if (inputs.size() != network_.reference.externalInputs.size()) {
    auto result = state;
    result.error = one();
    return result;
  }
  validate(inputs, inputs.size());
  auto result = state;
  result.previous = state.current;
  result.active.assign(network_.primitives.size(), zero());
  for (size_t i = 0; i < inputs.size(); ++i) {
    const auto net = network_.reference.externalInputs[i];
    result.current[net] = inputs[i];
    auto* changed = different(inputs[i], state.current[net]);
    for (auto consumer : consumers_[net]) result.active[consumer] = disjunction(result.active[consumer], changed);
  }
  normalizeCompletedActivation(result, consumers_);
  auto invalid = state;
  invalid.error = one();
  return choose(state.error, state, choose(stable(state), result, invalid));
}

SymbolicState SymbolicEventModel::wave(const SymbolicState& state, const FreshSymbol& fresh) const {
  const auto& reference = network_.reference;
  validate(reference, state);
  if (state.error == one() || stable(state) == one()) return state;
  // Allocate every nondeterministic symbol before entering parallel workers.
  std::vector<SymbolicBits> selections(reference.primitives.size());
  std::set<size_t> allocated;
  for (size_t i = 0; i < selections.size(); ++i) {
    const auto& primitive = reference.primitives[i];
    if (!primitive.storageBits || state.bootstrap == one() || state.active[i] == zero()) continue;
    const auto count = permutationCount(primitive.inputs.size(), limits_.maxPinOrderings);
    for (size_t remaining = count - 1; remaining; remaining >>= 1) {
      if (!fresh) throw std::invalid_argument("Missing symbolic choice allocator");
      auto* symbol = fresh();
      validate({symbol}, 1);
      if (symbol->getOp() != Op::VAR ||
          (symbol->getId() >= 2 && !allocated.insert(symbol->getId()).second))
        throw std::invalid_argument("Symbolic choice allocator must return distinct variables");
      selections[i].push_back(symbol);
    }
  }
  std::vector<SymbolicReaction> results(reference.primitives.size());
  const auto evaluate = [&](size_t index) {
    const auto& primitive = reference.primitives[index];
    const auto& callback = network_.primitives[index];
    const auto& storage = state.storage[index];
    const auto previous = gather(state.previous, primitive.inputs);
    const auto current = gather(state.current, primitive.inputs);
    const auto outputs = gather(state.current, primitive.outputs);
    SymbolicReaction held{storage, outputs, zero()};
    if (state.active[index] == zero()) { results[index] = held; return; }
    auto boot = invoke(primitive, callback, storage, previous, current, outputs, {}, true);
    if (state.bootstrap == one()) {
      results[index] = choose(state.active[index], boot, held);
      return;
    }
    SymbolicReaction ordinary;
    if (!primitive.storageBits) {
      ordinary = invoke(primitive, callback, storage, previous, current, outputs, {}, false);
    } else {
      auto* noChanges = one();
      SymbolicBits changed;
      for (size_t pin = 0; pin < current.size(); ++pin) {
        changed.push_back(different(current[pin], previous[pin]));
        noChanges = conjunction(noChanges, negate(changed.back()));
      }
      const auto fallback = invoke(primitive, callback, storage, previous, current, outputs, {}, false);
      std::vector<size_t> permutation(current.size());
      std::iota(permutation.begin(), permutation.end(), 0);
      size_t ordinal = 0;
      do {
        auto pins = previous;
        auto result = held;
        for (auto pin : permutation) {
          const auto before = pins;
          pins[pin] = current[pin];
          auto reaction = invoke(primitive, callback, result.storage, before, pins, outputs, pin, false);
          auto* visit = conjunction(changed[pin], negate(result.error));
          result = choose(visit, reaction, result);
        }
        result = choose(noChanges, fallback, result);
        // All unused binary rank encodings select the canonical first order.
        ordinary = ordinal == 0 ? result : choose(code(selections[index], ordinal), result, ordinary);
        ++ordinal;
      } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    results[index] = choose(state.active[index], choose(state.bootstrap, boot, ordinary), held);
  };
  tbb::task_arena arena(workers_ ? static_cast<int>(workers_) : tbb::task_arena::automatic);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, results.size()), [&](const auto& range) {
      for (size_t i = range.begin(); i != range.end(); ++i) evaluate(i);
    });
  });
  auto result = state;
  result.bootstrap = zero();
  result.previous = state.current;
  result.active.assign(results.size(), zero());
  for (size_t i = 0; i < results.size(); ++i) {
    result.storage[i] = results[i].storage;
    result.error = disjunction(result.error, results[i].error);
    for (size_t bit = 0; bit < results[i].outputs.size(); ++bit)
      result.current[reference.primitives[i].outputs[bit]] = results[i].outputs[bit];
  }
  for (size_t net = 0; net < reference.netCount; ++net) {
    auto* changed = different(result.current[net], state.current[net]);
    for (auto consumer : consumers_[net]) result.active[consumer] = disjunction(result.active[consumer], changed);
  }
  normalizeCompletedActivation(result, consumers_);
  // Stable identity needs no global mux: every inactive primitive already
  // holds, BOOT is false and previous=current, so every field above is unchanged
  // on a stable valuation. Only errors require explicit absorbing padding.
  return choose(state.error, state, result);
}

BoolExpr* SymbolicEventModel::boundaryInvariant(const SymbolicState& state) const {
  auto* result = stable(state);
  const auto& reference = network_.reference;
  for (size_t net = 0; net < reference.netCount; ++net)
    if (reference.constantByNet[net])
      result = conjunction(result, *reference.constantByNet[net] ? state.current[net] : negate(state.current[net]));
  for (size_t i = 0; i < reference.primitives.size(); ++i) {
    const auto& primitive = reference.primitives[i];
    const auto pins = gather(state.current, primitive.inputs);
    const auto outputs = gather(state.current, primitive.outputs);
    const auto reaction = invoke(primitive, network_.primitives[i], state.storage[i], pins, pins, outputs, {}, true);
    result = conjunction(result, negate(reaction.error));
    result = conjunction(result, equal(state.storage[i], reaction.storage));
    result = conjunction(result, equal(outputs, reaction.outputs));
  }
  return result;
}

SymbolicPrimitive liftPrimitive(const Primitive& primitive, size_t maxLocalBits) {
  const auto storageBits = primitive.storageBits;
  const auto pins = primitive.inputs.size();
  const auto total = domainSize(storageBits, pins, 2, maxLocalBits);
  SymbolicPrimitive result;
  result.react = [primitive, storageBits, pins, total](const SymbolicBits& storage,
      const SymbolicBits& before, const SymbolicBits& current,
      std::optional<size_t> changed, bool bootstrap) {
    validate(storage, storageBits); validate(before, pins); validate(current, pins);
    SymbolicBits arguments = storage;
    arguments.insert(arguments.end(), before.begin(), before.end());
    arguments.insert(arguments.end(), current.begin(), current.end());
    SymbolicReaction lifted{SymbolicBits(storageBits, zero()), SymbolicBits(primitive.outputs.size(), zero()), zero()};
    for (size_t assignment = 0; assignment < (size_t{1} << total); ++assignment) {
      auto* condition = minterm(arguments, assignment);
      if (condition == zero()) continue;
      Reaction concrete;
      try {
        if (!primitive.react) throw std::invalid_argument("missing reaction");
        concrete = primitive.react(unpack(assignment, 0, storageBits),
            unpack(assignment, storageBits, pins), unpack(assignment, storageBits + pins, pins), changed, bootstrap);
        if (concrete.storage.size() != storageBits || concrete.outputs.size() != primitive.outputs.size() ||
            !binary(concrete.storage) || !binary(concrete.outputs)) concrete.error = true;
      } catch (const Limit&) { throw;
      } catch (const std::bad_alloc&) { throw;
      } catch (const std::exception&) { concrete.error = true; }
      if (concrete.error) { lifted.error = disjunction(lifted.error, condition); continue; }
      for (size_t i = 0; i < storageBits; ++i)
        if (concrete.storage[i]) lifted.storage[i] = disjunction(lifted.storage[i], condition);
      for (size_t i = 0; i < concrete.outputs.size(); ++i)
        if (concrete.outputs[i]) lifted.outputs[i] = disjunction(lifted.outputs[i], condition);
    }
    return lifted;
  };
  const auto initialBits = domainSize(storageBits, pins, 1, maxLocalBits);
  result.initialOutputs = [primitive, storageBits, pins, initialBits](const SymbolicBits& storage,
                                                                   const SymbolicBits& inputs) {
    validate(storage, storageBits); validate(inputs, pins);
    SymbolicBits arguments = storage;
    arguments.insert(arguments.end(), inputs.begin(), inputs.end());
    SymbolicBits outputs(primitive.outputs.size(), zero());
    for (size_t assignment = 0; assignment < (size_t{1} << initialBits); ++assignment) {
      auto* condition = minterm(arguments, assignment);
      if (condition == zero()) continue;
      const auto bits = unpack(assignment, 0, storageBits);
      Bits values;
      if (primitive.initialOutputValues) values = primitive.initialOutputValues(bits, unpack(assignment, storageBits, pins));
      else if (primitive.initialOutputs) values = primitive.initialOutputs(bits);
      else if (primitive.outputs.size() == storageBits) values = bits;
      else throw std::invalid_argument("missing initial output mapping");
      if (values.size() != outputs.size() || !binary(values))
        throw std::invalid_argument("invalid initial output mapping");
      for (size_t i = 0; i < values.size(); ++i)
        if (values[i]) outputs[i] = disjunction(outputs[i], condition);
    }
    return outputs;
  };
  return result;
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
