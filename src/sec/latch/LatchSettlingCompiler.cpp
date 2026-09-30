// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "LatchSettlingCompiler.h"

#include <algorithm>
#include <limits>
#include <new>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

struct Failure {
  CertificationStatus status;
  std::string detail;
};

[[noreturn]] void invalid(const std::string& detail) {
  throw Failure{CertificationStatus::Invalid, detail};
}

[[noreturn]] void resource(const std::string& detail) {
  throw Failure{CertificationStatus::ResourceLimit, detail};
}

void addWidth(size_t& total, size_t width, size_t maximum) {
  if (total > maximum || width > maximum - total) {
    resource("Complete reference state exceeds the state-bit limit");
  }
  total += width;
}

void checkBits(const Bits& bits) {
  if (std::any_of(bits.begin(), bits.end(), [](uint8_t bit) { return bit > 1; })) {
    invalid("A reference value is not Boolean");
  }
}

void checkState(const State& state, const State& shape,
                const CompilerLimits& limits) {
  if (state.current.size() != shape.current.size() ||
      state.previous.size() != state.current.size() ||
      state.storage.size() != shape.storage.size() ||
      state.active.size() != shape.active.size() ||
      state.active.size() != state.storage.size()) {
    invalid("Reference transitions changed the complete-state layout");
  }
  size_t bits = 0;
  addWidth(bits, 2, limits.maxStateBits);  // BOOT and sticky error.
  addWidth(bits, state.current.size(), limits.maxStateBits);
  addWidth(bits, state.previous.size(), limits.maxStateBits);
  addWidth(bits, state.active.size(), limits.maxStateBits);
  checkBits(state.current);
  checkBits(state.previous);
  checkBits(state.active);
  for (size_t primitive = 0; primitive < state.storage.size(); ++primitive) {
    if (state.storage[primitive].size() != shape.storage[primitive].size()) {
      invalid("Reference transitions changed the primitive storage layout");
    }
    addWidth(bits, state.storage[primitive].size(), limits.maxStateBits);
    checkBits(state.storage[primitive]);
  }
}

void checkNetworkSize(const Network& network, const CompilerLimits& limits) {
  size_t bits = 0;
  addWidth(bits, 2, limits.maxStateBits);
  addWidth(bits, network.netCount, limits.maxStateBits);
  addWidth(bits, network.netCount, limits.maxStateBits);
  addWidth(bits, network.primitives.size(), limits.maxStateBits);
  for (const auto& primitive : network.primitives) {
    addWidth(bits, primitive.storageBits, limits.maxStateBits);
  }
}

size_t valuationCount(size_t bits, size_t bitLimit, const char* description) {
  if (bits > bitLimit || bits >= std::numeric_limits<size_t>::digits) {
    resource(std::string(description) + " exceeds the Boolean-enumeration limit");
  }
  return size_t{1} << bits;
}

Bits valuation(size_t code, size_t width) {
  Bits result(width);
  for (size_t bit = 0; bit < width; ++bit) {
    result[bit] = static_cast<uint8_t>((code >> bit) & size_t{1});
  }
  return result;
}

std::vector<size_t> bootstrapSeedNets(const Network& network) {
  std::vector<size_t> result;
  const bool inputDependentProjection = std::any_of(
      network.primitives.begin(), network.primitives.end(),
      [](const Primitive& primitive) { return bool(primitive.initialOutputValues); });
  for (const auto& primitive : network.primitives) {
    // Input-dependent BOOT projections read one frozen preprojection snapshot.
    // A storage output's seed can therefore influence another cell even though
    // that output is overwritten by its own storage projection afterward.
    if (primitive.storageBits == 0 || inputDependentProjection) {
      result.insert(result.end(), primitive.outputs.begin(), primitive.outputs.end());
    }
  }
  // The EventModel constructor enforces unique ownership; sorting also makes
  // enumeration independent of how the caller listed the primitives.
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

struct GraphNode {
  State state;
  std::vector<size_t> next;
  bool stable = false;
};

void exploreEpisode(const std::vector<State>& entries,
                    const EpisodeRelation& relation,
                    const CompilerLimits& limits, EpisodeResult& result) {
  if (entries.empty()) {
    invalid("An empty episode-entry relation cannot establish a certificate");
  }
  if (!relation.stable || !relation.successors) {
    invalid("The episode reference relation is missing a required callback");
  }
  std::vector<GraphNode> nodes;
  std::unordered_map<std::string, size_t> byState;
  std::vector<size_t> roots;
  const auto insert = [&](const State& state) {
    checkState(state, entries.front(), limits);
    auto key = state.key();
    const auto found = byState.find(key);
    if (found != byState.end()) {
      return found->second;
    }
    if (nodes.size() >= limits.maxEpisodeStates) {
      resource("Episode graph exceeds the complete-state limit");
    }
    const size_t index = nodes.size();
    byState.emplace(std::move(key), index);
    nodes.push_back({state, {}, false});
    result.exploredStates = nodes.size();
    return index;
  };
  for (const auto& entry : entries) {
    roots.push_back(insert(entry));
  }

  for (size_t index = 0; index < nodes.size(); ++index) {
    // Insertion below can reallocate nodes; never keep a reference into it.
    const State state = nodes[index].state;
    if (state.error) {
      invalid("Reachable reference error: " + state.errorReason);
    }
    const bool stable = relation.stable(state);
    nodes[index].stable = stable;
    const auto successors = relation.successors(state);
    if (successors.empty()) {
      invalid("A reachable reference state has no successor (missing totality)");
    }
    if (successors.size() >
        limits.maxEpisodeTransitions - result.exploredTransitions) {
      resource("Episode graph exceeds the transition limit");
    }
    result.exploredTransitions += successors.size();
    if (stable) {
      for (const auto& successor : successors) {
        if (successor != state) {
          invalid("Stable reference states must have only identity successors");
        }
      }
      result.stableStates.push_back(state);
      continue;  // Identity padding does not count toward the settling depth.
    }
    std::vector<size_t> next;
    for (const auto& successor : successors) {
      next.push_back(insert(successor));
    }
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
    nodes[index].next = std::move(next);
  }

  // Remove stable identity edges and topologically sort the remaining graph.
  // A cycle is an allowed infinite execution, even if other branches settle.
  std::vector<size_t> incoming(nodes.size(), 0);
  for (const auto& node : nodes) {
    for (const size_t next : node.next) {
      ++incoming[next];
    }
  }
  std::vector<size_t> order;
  for (size_t index = 0; index < nodes.size(); ++index) {
    if (incoming[index] == 0) {
      order.push_back(index);
    }
  }
  for (size_t cursor = 0; cursor < order.size(); ++cursor) {
    for (const size_t next : nodes[order[cursor]].next) {
      if (--incoming[next] == 0) {
        order.push_back(next);
      }
    }
  }
  if (order.size() != nodes.size()) {
    throw Failure{CertificationStatus::NonSettling,
                  "A reachable nonstable cycle admits an infinite episode"};
  }

  std::vector<size_t> depth(nodes.size(), 0);
  for (auto position = order.rbegin(); position != order.rend(); ++position) {
    for (const size_t next : nodes[*position].next) {
      depth[*position] = std::max(depth[*position], depth[next] + 1);
    }
  }
  for (const size_t root : roots) {
    result.maxWaves = std::max(result.maxWaves, depth[root]);
  }
  if (result.stableStates.empty()) {
    invalid("The reference graph has no stable result");
  }
  std::sort(result.stableStates.begin(), result.stableStates.end());
  if (result.stableStates.size() != 1) {
    throw Failure{CertificationStatus::OrderDependent,
                  "Allowed executions reach different complete boundary states"};
  }
  if (result.maxWaves > limits.maxWaves) {
    throw Failure{CertificationStatus::UnprovedBound,
                  "The exact settling depth exceeds the permitted unfolding bound"};
  }
  result.status = CertificationStatus::Certified;
}

EpisodeResult bootstrapImpl(const EventModel& model, const Bits& inputs,
                            const std::vector<Bits>& initialStorage,
                            const CompilerLimits& limits) {
  checkNetworkSize(model.network(), limits);
  const auto seedNets = bootstrapSeedNets(model.network());
  const size_t count = valuationCount(seedNets.size(), limits.maxBootstrapSeedBits,
                                     "Auxiliary bootstrap seed space");
  if (count > limits.maxInitialConfigurations || count > limits.maxEpisodeStates) {
    resource("Bootstrap seed space exceeds the initial-configuration limit");
  }
  std::vector<State> entries;
  entries.reserve(count);
  for (size_t code = 0; code < count; ++code) {
    Bits seeds(model.network().netCount, 0);
    for (size_t bit = 0; bit < seedNets.size(); ++bit) {
      seeds[seedNets[bit]] = static_cast<uint8_t>((code >> bit) & size_t{1});
    }
    entries.push_back(model.bootstrap(inputs, initialStorage, seeds));
  }
  return certifyEpisode(model, entries, limits);
}

template <class Result, class Work>
Result guarded(Work&& work) {
  Result result;
  try {
    work(result);
  } catch (const Failure& failure) {
    result.status = failure.status;
    result.detail = failure.detail;
  } catch (const Limit& limit) {
    result.status = CertificationStatus::ResourceLimit;
    result.detail = limit.what();
  } catch (const std::bad_alloc&) {
    result.status = CertificationStatus::ResourceLimit;
    result.detail = "Memory exhausted while constructing the complete certificate";
  } catch (const std::exception& error) {
    result.status = CertificationStatus::Invalid;
    result.detail = error.what();
  }
  return result;
}

void requireCertificate(const EpisodeResult& episode, const std::string& context) {
  if (!episode.certified()) {
    throw Failure{episode.status, context + ": " + episode.detail};
  }
}

}  // namespace

const char* certificationStatusName(CertificationStatus status) {
  switch (status) {
    case CertificationStatus::Certified:
      return "certified";
    case CertificationStatus::NonSettling:
      return "non-settling";
    case CertificationStatus::OrderDependent:
      return "order-dependent";
    case CertificationStatus::Invalid:
      return "invalid-reference";
    case CertificationStatus::ResourceLimit:
      return "resource-limit";
    case CertificationStatus::UnprovedBound:
      return "unproved-bound";
  }
  return "invalid-status";
}

EpisodeResult certifyEpisode(const std::vector<State>& entries,
                             const EpisodeRelation& relation,
                             const CompilerLimits& limits) {
  return guarded<EpisodeResult>([&](EpisodeResult& result) {
    exploreEpisode(entries, relation, limits, result);
  });
}

EpisodeResult certifyEpisode(const EventModel& model,
                             const std::vector<State>& entries,
                             const CompilerLimits& limits) {
  return certifyEpisode(entries,
                        {[&](const State& state) { return model.stable(state); },
                         [&](const State& state) { return model.successors(state); }},
                        limits);
}

EpisodeResult certifyBootstrap(const EventModel& model, const Bits& inputs,
                               const std::vector<Bits>& initialStorage,
                               const CompilerLimits& limits) {
  return guarded<EpisodeResult>([&](EpisodeResult& result) {
    result = bootstrapImpl(model, inputs, initialStorage, limits);
  });
}

CompileResult compileTransitionTable(const EventModel& model,
                                     const CompileOptions& options) {
  return guarded<CompileResult>([&](CompileResult& result) {
    const auto& network = model.network();
    const auto& limits = options.limits;
    checkNetworkSize(network, limits);
    checkBits(options.initialInputs);
    if (options.initialInputs.size() != network.externalInputs.size()) {
      invalid("Initial external valuation has the wrong width");
    }
    const size_t inputCount = valuationCount(network.externalInputs.size(),
                                             limits.maxExternalBits,
                                             "External input space");
    if (!options.initialStorage.empty() &&
        options.initialStorage.size() != network.primitives.size()) {
      invalid("Initial storage relation has the wrong primitive count");
    }
    std::vector<Bits> baseStorage;
    std::vector<std::pair<size_t, size_t>> unknown;
    for (size_t primitive = 0; primitive < network.primitives.size(); ++primitive) {
      const size_t width = network.primitives[primitive].storageBits;
      baseStorage.emplace_back(width, 0);
      if (!options.initialStorage.empty() &&
          options.initialStorage[primitive].size() != width) {
        invalid("Initial storage relation has the wrong primitive width");
      }
      for (size_t bit = 0; bit < width; ++bit) {
        const auto value = options.initialStorage.empty()
                               ? std::optional<uint8_t>{}
                               : options.initialStorage[primitive][bit];
        if (!value) {
          unknown.emplace_back(primitive, bit);
        } else if (*value > 1) {
          invalid("Initial storage relation contains a non-Boolean bit");
        } else {
          baseStorage.back()[bit] = *value;
        }
      }
    }
    const size_t initialCount = valuationCount(unknown.size(),
                                               limits.maxInitialStorageBits,
                                               "Initial storage space");
    const size_t seedCount = valuationCount(bootstrapSeedNets(network).size(),
                                            limits.maxBootstrapSeedBits,
                                            "Auxiliary bootstrap seed space");
    if (initialCount > limits.maxInitialConfigurations ||
        seedCount > limits.maxInitialConfigurations / initialCount) {
      resource("Complete initialization relation exceeds the configuration limit");
    }

    TransitionTable table;
    table.initialInputs = options.initialInputs;
    table.singleExternalInputChange = options.singleExternalInputChange;
    std::unordered_map<std::string, size_t> byBoundary;
    const auto insertBoundary = [&](const State& boundary) {
      if (!model.stable(boundary)) {
        invalid("Only complete stable states may become macro boundaries");
      }
      auto key = boundary.key();
      const auto found = byBoundary.find(key);
      if (found != byBoundary.end()) {
        return found->second;
      }
      if (table.boundaries.size() >= limits.maxBoundaryStates) {
        resource("Reachable complete-boundary set exceeds its state limit");
      }
      const size_t index = table.boundaries.size();
      byBoundary.emplace(std::move(key), index);
      table.boundaries.push_back(boundary);
      return index;
    };

    for (size_t code = 0; code < initialCount; ++code) {
      auto storage = baseStorage;
      for (size_t bit = 0; bit < unknown.size(); ++bit) {
        const auto [primitive, position] = unknown[bit];
        storage[primitive][position] =
            static_cast<uint8_t>((code >> bit) & size_t{1});
      }
      const auto episode = certifyBootstrap(model, options.initialInputs, storage, limits);
      requireCertificate(episode, "Bootstrap origin " + std::to_string(code));
      table.maxWaves = std::max(table.maxWaves, episode.maxWaves);
      const size_t boundary = insertBoundary(episode.stableStates.front());
      table.initials.push_back({std::move(storage), boundary});
    }

    // BFS closure is the boundary invariant. No externally allowed input may be
    // skipped because it fails to settle or yields an inconvenient next state.
    for (size_t from = 0; from < table.boundaries.size(); ++from) {
      const State boundary = table.boundaries[from];
      for (size_t code = 0; code < inputCount; ++code) {
        Bits input = valuation(code, network.externalInputs.size());
        if (options.singleExternalInputChange) {
          size_t changes = 0;
          for (size_t bit = 0; bit < input.size(); ++bit) {
            changes += input[bit] != boundary.current[network.externalInputs[bit]];
          }
          if (changes > 1) {
            continue;
          }
        }
        if (table.rows.size() >= limits.maxTransactions) {
          resource("Reachable macro transition table exceeds its transaction limit");
        }
        const State entry = model.admit(boundary, input);
        const auto episode = certifyEpisode(model, {entry}, limits);
        requireCertificate(episode, "Boundary " + std::to_string(from) +
                                        ", external valuation " + std::to_string(code));
        table.maxWaves = std::max(table.maxWaves, episode.maxWaves);
        const size_t to = insertBoundary(episode.stableStates.front());
        table.rows.push_back({from, std::move(input), to});
      }
    }
    result.status = CertificationStatus::Certified;
    result.table = std::move(table);
  });
}

}  // namespace KEPLER_FORMAL::SEC::LATCH
