// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "LatchEventModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {

enum class CertificationStatus {
  Certified,
  NonSettling,
  OrderDependent,
  Invalid,
  ResourceLimit,
  UnprovedBound,
};

const char* certificationStatusName(CertificationStatus status);

// Every limit is a rejection threshold, never an assumption removing executions.
// maxWaves limits compilation depth; cycle detection explores complete states
// independently of that candidate depth, within the graph resource limits.
struct CompilerLimits {
  size_t maxStateBits = 512;
  size_t maxEpisodeStates = 65536;
  size_t maxEpisodeTransitions = 262144;
  size_t maxBoundaryStates = 4096;
  size_t maxTransactions = 65536;
  size_t maxInitialConfigurations = 4096;
  size_t maxExternalBits = 12;
  size_t maxInitialStorageBits = 12;
  size_t maxBootstrapSeedBits = 12;
  size_t maxWaves = 256;
};

// The generic interface also permits testing totality, absorbing stability, and
// complete-state uniqueness independently of primitive-specific construction.
// Callbacks describe the entire relation, not one sampled execution.
struct EpisodeRelation {
  std::function<bool(const State&)> stable;
  std::function<std::vector<State>(const State&)> successors;
};

struct EpisodeResult {
  CertificationStatus status = CertificationStatus::Invalid;
  std::string detail;
  // Complete stable states, not merely equal present outputs. Certified implies
  // exactly one result. On failure this vector is diagnostic, not a certificate.
  std::vector<State> stableStates;
  size_t maxWaves = 0;
  size_t exploredStates = 0;
  size_t exploredTransitions = 0;

  bool certified() const { return status == CertificationStatus::Certified; }
};

// Explore the complete reachable episode graph. A nonstable reachable cycle is
// NonSettling, a missing/error transition Invalid, and distinct full boundaries
// OrderDependent. Acyclic graphs give an exact longest wave count. Stable states
// must have identity successors; they are not permitted to vanish from a bound.
EpisodeResult certifyEpisode(const std::vector<State>& entries,
                             const EpisodeRelation& relation,
                             const CompilerLimits& limits = {});
EpisodeResult certifyEpisode(const EventModel& model,
                             const std::vector<State>& entries,
                             const CompilerLimits& limits = {});

// Fix intended storage and external inputs, independently enumerate ALL Boolean
// seeds of combinational outputs, and certify their joint set of BOOT entries.
// If a storage-output BOOT projection reads pins, enumerate every internal output
// seed: another storage output's overwritten seed can affect that projection.
EpisodeResult certifyBootstrap(const EventModel& model, const Bits& inputs,
                               const std::vector<Bits>& initialStorage,
                               const CompilerLimits& limits = {});

struct CompileOptions {
  Bits initialInputs;
  // One vector per primitive, including empty vectors for zero-storage gates.
  // nullopt means a fixed-once, universally included initial Boolean choice.
  // An empty outer vector means every primitive's storage bit is unspecified.
  std::vector<std::vector<std::optional<uint8_t>>> initialStorage;
  // An explicit environment restriction, OFF by default. It constrains only
  // external transactions, never simultaneous internally generated pin changes.
  bool singleExternalInputChange = false;
  CompilerLimits limits;
};

struct BoundaryTransition {
  size_t from = 0;
  Bits input;
  size_t to = 0;
};

struct InitialBoundary {
  // Retain every prescribed origin, even if different storage choices settle
  // to the same boundary. Never cherry-pick a convenient initialization pair.
  std::vector<Bits> storage;
  size_t boundary = 0;
};

struct TransitionTable {
  std::vector<State> boundaries;
  std::vector<BoundaryTransition> rows;
  std::vector<InitialBoundary> initials;
  Bits initialInputs;
  bool singleExternalInputChange = false;
  size_t maxWaves = 0;
};

struct CompileResult {
  CertificationStatus status = CertificationStatus::Invalid;
  std::string detail;
  // Present ONLY after all initial origins and every permitted transaction from
  // every reachable complete boundary have passed progress and uniqueness.
  std::optional<TransitionTable> table;

  bool certified() const { return status == CertificationStatus::Certified; }
};

CompileResult compileTransitionTable(const EventModel& model,
                                     const CompileOptions& options);

}  // namespace KEPLER_FORMAL::SEC::LATCH
