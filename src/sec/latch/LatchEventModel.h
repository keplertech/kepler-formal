// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace KEPLER_FORMAL::SEC::LATCH {

using Bits = std::vector<uint8_t>;

struct Reaction {
  Bits storage;
  Bits outputs;
  bool error = false;
  std::string reason{};
};

struct Primitive {
  std::string name;
  std::vector<size_t> inputs;
  std::vector<size_t> outputs;
  size_t storageBits = 0;
  // Callbacks must be pure and thread-safe. A sequential reaction consumes only
  // changedPin's event; bootstrap has no edge. Zero-storage cells evaluate once
  // per activation, with changedPin unset and the complete current pin vector.
  std::function<Reaction(const Bits&, const Bits&, const Bits&,
                        std::optional<size_t> changedPin, bool bootstrap)> react;
  // Maps remembered storage to physical outputs before forced bootstrap. If
  // omitted, identity is permitted only when the widths match. Gates use seeds.
  std::function<Bits(const Bits&)> initialOutputs;
  // Optional input-dependent physical projection (for example, an integrated
  // clock gate). All such projections read one frozen seed/input snapshot and
  // commit together. Auxiliary seeds still require bootstrap certification.
  std::function<Bits(const Bits& storage, const Bits& pins)> initialOutputValues;
};

struct Network {
  size_t netCount = 0;
  std::vector<size_t> externalInputs;
  std::vector<Primitive> primitives;
  // Empty means no constants; otherwise exactly netCount entries.
  std::vector<std::optional<bool>> constantByNet{};
};

struct State {
  Bits current;
  Bits previous;
  std::vector<Bits> storage;
  Bits active;
  bool bootstrap = false;
  bool error = false;
  std::string errorReason;

  bool operator==(const State& other) const;
  bool operator!=(const State& other) const { return !(*this == other); }
  bool operator<(const State& other) const;
  // Complete, unambiguous serialization suitable for hash-map keys.
  std::string key() const;
};

struct Limits {
  size_t maxPinOrderings = 100000;
  size_t maxSuccessors = 100000;
};

// Resource exhaustion never silently truncates the reference relation.
class Limit : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Finite Boolean evaluate/update semantics. Every wave is an epoch barrier:
// primitive evaluations may run concurrently but all read the same snapshot.
// Simultaneous changed pin positions retain ALL orders and producer choices are
// committed once, shared by every fanout. Worker count cannot change behavior.
class EventModel {
 public:
  explicit EventModel(Network network, Limits limits = {}, size_t workerCount = 0);

  const Network& network() const { return network_; }
  // inputs follows externalInputs order. initialStorage has one entry per cell,
  // including empty entries for gates. Seeds has netCount bits. External and
  // constant positions are overwritten. Physical storage outputs are projected
  // before BOOT; with input-dependent projections their original seeds can
  // affect other projections and must also be quantified by the certifier.
  State bootstrap(const Bits& inputs, const std::vector<Bits>& initialStorage,
                  const Bits& internalSeeds) const;
  // Invalid admission is an explicit absorbing, nonstable error state.
  State admit(const State& quiescent, const Bits& inputs) const;
  std::vector<State> successors(const State& state) const;
  bool stable(const State& state) const;

 private:
  void validateState(const State& state) const;
  std::vector<Reaction> evaluate(size_t primitive, const State& state) const;
  void normalizeBoundary(State& state) const;

  Network network_;
  Limits limits_;
  size_t workerCount_;
  std::vector<std::vector<size_t>> consumers_;
};

// Convenience primitives for tests and explicitly modeled Boolean cells. More
// elaborate reset priorities/output mappings use the Primitive reaction API.
Primitive combinational(std::string name, std::vector<size_t> inputs,
                        std::vector<size_t> outputs,
                        std::function<Bits(const Bits&)> function);
Primitive latch(std::string name, size_t data, size_t enable, size_t output,
                bool activeHigh = true);
Primitive flipFlop(std::string name, size_t data, size_t clock, size_t output,
                   bool risingEdge = true);

}  // namespace KEPLER_FORMAL::SEC::LATCH
