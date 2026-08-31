// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "phase/PeriodicSignalAnalysis.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace KEPLER_FORMAL::SEC::PHASE {

namespace {

size_t effectiveConcurrency(size_t requested) {
  if (std::getenv("KEPLER_NO_MT") != nullptr) {
    return 1;
  }
  return requested;
}

template <typename Function>
void parallelIndexes(size_t count, size_t concurrency, Function&& function) {
  if (count == 0) {
    return;
  }
  const size_t effective = effectiveConcurrency(concurrency);
  if (effective == 1 || count == 1) {
    function(0, count);
    return;
  }
  const int arenaConcurrency =
      effective == 0
          ? tbb::task_arena::automatic
          : static_cast<int>(std::min<size_t>(
                effective,
                static_cast<size_t>(std::numeric_limits<int>::max())));
  tbb::task_arena arena(arenaConcurrency);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, count),
                      [&](const tbb::blocked_range<size_t>& range) {
                        function(range.begin(), range.end());
                      });
  });
}

size_t minimumPeriod(const TernarySimulationResult& simulation,
                     size_t stateIndex) {
  const size_t cycle = simulation.cycleLength;
  for (size_t period = 1; period <= cycle; ++period) {
    if (cycle % period != 0) {
      continue;
    }
    bool matches = true;
    // Generator words are reset-anchored. They must explain the transient
    // stem as well as the eventual cycle, not just the post-convergence cycle.
    for (size_t step = period; step < simulation.trace.size(); ++step) {
      if (simulation.trace[step][stateIndex] !=
          simulation.trace[step % period][stateIndex]) {
        matches = false;
        break;
      }
    }
    if (matches) {
      return period;
    }
  }
  return 0;
}

}  // namespace

PeriodicSignalAnalysisResult analyzePeriodicSignals(
    const NormalizedTransitionSystem& system,
    const TernarySimulationResult& simulation,
    const PeriodicSignalAnalysisOptions& options) {
  PeriodicSignalAnalysisResult result;
  result.statistics.candidateStateCount = system.states().size();
  if (options.maxPhases == 0) {
    result.complete = false;
    result.diagnostic = "maxPhases must be at least one";
    return result;
  }
  if (simulation.status != TernarySimulationStatus::Complete ||
      simulation.cycleLength == 0) {
    result.complete = false;
    result.diagnostic =
        "periodic-signal analysis requires a completed reset-to-cycle trace";
    return result;
  }
  if (simulation.stemLength + simulation.cycleLength >
      simulation.trace.size()) {
    result.complete = false;
    result.diagnostic = "ternary trace has an invalid stem/cycle range";
    return result;
  }
  for (const auto& state : simulation.trace) {
    if (state.size() != system.states().size()) {
      result.complete = false;
      result.diagnostic = "ternary trace state width does not match the model";
      return result;
    }
  }

  struct Candidate {
    bool valid = false;
    size_t period = 0;
    std::vector<TernaryValue> word;
  };
  std::vector<Candidate> candidates(system.states().size());
  parallelIndexes(
      candidates.size(), options.maxConcurrency, [&](size_t begin, size_t end) {
        for (size_t stateIndex = begin; stateIndex < end; ++stateIndex) {
          // Bjesse/Kukula Section III removes a generator candidate when it is
          // X anywhere in the reset stem or cycle.
          bool boolean = true;
          for (const auto& traceState : simulation.trace) {
            if (traceState[stateIndex] == TernaryValue::Unknown) {
              boolean = false;
              break;
            }
          }
          if (!boolean) {
            continue;
          }
          Candidate candidate;
          candidate.period = minimumPeriod(simulation, stateIndex);
          if (candidate.period == 0) {
            continue;
          }
          candidate.valid = true;
          candidate.word.reserve(candidate.period);
          for (size_t offset = 0; offset < candidate.period; ++offset) {
            candidate.word.push_back(simulation.trace[offset][stateIndex]);
          }
          candidates[stateIndex] = std::move(candidate);
        }
      });

  for (size_t stateIndex = 0; stateIndex < candidates.size(); ++stateIndex) {
    auto& candidate = candidates[stateIndex];
    if (!candidate.valid) {
      continue;
    }
    result.generators.push_back(
        {stateIndex, candidate.period, std::move(candidate.word), {}, false});
  }
  result.statistics.resetPeriodicStateCount = result.generators.size();

  size_t bestPhaseCount = 1;
  size_t bestUsable = 0;
  // Section IV's practical heuristic chooses a small global unfolding factor
  // that makes as many proven generator words usable as possible.
  // Every candidate period divides the detected cycle, so considering an N
  // above the cycle length can never improve on the cycle itself. This also
  // keeps an accidentally huge caller cap from turning into a huge loop.
  const size_t phaseSearchLimit =
      std::min(options.maxPhases, simulation.cycleLength);
  for (size_t phaseCount = 1; phaseCount <= phaseSearchLimit; ++phaseCount) {
    const size_t usable = static_cast<size_t>(
        std::count_if(result.generators.begin(), result.generators.end(),
                      [&](const PeriodicGenerator& generator) {
                        return generator.minimumPeriod != 0 &&
                               phaseCount % generator.minimumPeriod == 0;
                      }));
    if (usable > bestUsable) {
      bestUsable = usable;
      bestPhaseCount = phaseCount;
    }
  }
  result.phaseCount = bestPhaseCount;
  result.statistics.usableGeneratorCount = bestUsable;
  for (auto& generator : result.generators) {
    generator.usableAtSelectedPhaseCount =
        generator.minimumPeriod != 0 &&
        result.phaseCount % generator.minimumPeriod == 0;
    if (!generator.usableAtSelectedPhaseCount) {
      continue;
    }
    generator.expandedWord.reserve(result.phaseCount);
    for (size_t phase = 0; phase < result.phaseCount; ++phase) {
      generator.expandedWord.push_back(
          generator.generatorWord[phase % generator.minimumPeriod]);
    }
  }
  return result;
}

}  // namespace KEPLER_FORMAL::SEC::PHASE
