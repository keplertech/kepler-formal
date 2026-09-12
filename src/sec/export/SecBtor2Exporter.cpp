// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "export/SecBtor2Exporter.h"

#include <algorithm>
#include <bit>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "Btor2Writer.h"
#include "kinduction/KInductionProblem.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {
namespace {

using NodeId = Btor2Writer::NodeId;

// A saturating binary counter made exclusively from one-bit nodes. Even a long
// reset prefix needs only logarithmically many monitor state bits.
class ExportTimeline {
 public:
  ExportTimeline(Btor2Writer& writer, size_t lastFrame)
      : writer_(writer) {
    for (unsigned bit = 0; bit < std::bit_width(lastFrame); ++bit) {
      const auto node = writer_.state("kf_export_cycle_" + std::to_string(bit));
      bits_.push_back(node);
      writer_.init(node, writer_.constant(false));
    }
    if (!bits_.empty()) {
      const auto advance = writer_.logicalNot(at(lastFrame));
      auto carry = writer_.constant(true);
      for (const auto node : bits_) {
        writer_.next(node, writer_.logicalXor(
            node, writer_.logicalAnd(carry, advance)));
        carry = writer_.logicalAnd(carry, node);
      }
    }
  }

  NodeId at(size_t frame) {
    auto result = writer_.constant(true);
    for (size_t bit = 0; bit < bits_.size(); ++bit) {
      result = writer_.logicalAnd(result, (frame >> bit) & 1
          ? bits_[bit] : writer_.logicalNot(bits_[bit]));
    }
    return result;
  }

  NodeId before(size_t frame) {
    auto less = writer_.constant(false);
    auto equal = writer_.constant(true);
    for (size_t bit = bits_.size(); bit-- > 0;) {
      if ((frame >> bit) & 1) {
        less = writer_.logicalOr(less, writer_.logicalAnd(
            equal, writer_.logicalNot(bits_[bit])));
        equal = writer_.logicalAnd(equal, bits_[bit]);
      } else {
        equal = writer_.logicalAnd(equal, writer_.logicalNot(bits_[bit]));
      }
    }
    return less;
  }

  NodeId firstFrame() {
    if (!bits_.empty()) {
      return at(0);
    }
    if (first_ == 0) {
      first_ = writer_.state("kf_export_initial_frame");
      writer_.init(first_, writer_.constant(true));
      writer_.next(first_, writer_.constant(false));
    }
    return first_;
  }

 private:
  Btor2Writer& writer_;
  std::vector<NodeId> bits_;
  NodeId first_ = 0;
};

void constrainWhen(Btor2Writer& writer, NodeId when, NodeId condition) {
  writer.constraint(writer.logicalOr(writer.logicalNot(when), condition));
}

class TemporaryExport {
 public:
  explicit TemporaryExport(const std::filesystem::path& destination) {
    const auto parent = destination.has_parent_path()
        ? destination.parent_path() : std::filesystem::path(".");
    std::random_device random;
    for (unsigned attempt = 0; attempt < 32; ++attempt) {
      directory_ = parent / (".kf-btor2-" + std::to_string(random()) + "-" +
                             std::to_string(random()));
      if (std::filesystem::create_directory(directory_)) {
        return;
      }
    }
    throw std::runtime_error("Cannot create temporary BTOR2 export file");
  }
  ~TemporaryExport() {
    std::error_code error;
    std::filesystem::remove_all(directory_, error);
  }
  std::filesystem::path path() const { return directory_ / "model.btor2"; }

 private:
  std::filesystem::path directory_;
};

}  // namespace

void exportSecBtor2(const KInductionProblem& problem,
                   std::ostream& output,
                   const SecBtor2Metadata& metadata) {
  if (problem.bad == nullptr || problem.property == nullptr) {
    throw std::runtime_error("BTOR2 export requires a prepared SEC property");
  }
  Btor2Writer writer(output);
  // BTOR2 initialization values must precede their state declarations.
  writer.constant(false);
  writer.constant(true);
  writer.comment("Kepler Formal: prepared SEC equivalence obligation; no proof result");
  writer.comment("encoding=" + std::string(problem.usesDualRailStateEncoding
      ? "dual_rail_steady" : "binary"));
  writer.comment("startup=SEC concrete base-case observation semantics");
  writer.comment("covered_outputs=" + std::to_string(problem.observedOutputNames.size()) +
      " total_outputs=" + std::to_string(metadata.totalOutputCount != 0
          ? metadata.totalOutputCount : problem.observedOutputNames.size()));
  for (const auto& name : problem.observedOutputNames) {
    writer.comment("covered_output=" + name);
  }
  for (const auto& name : metadata.skippedOutputs) {
    writer.comment("skipped_output=" + name);
  }

  std::map<size_t, std::string> states;
  const auto addStates = [&](const auto& symbols, const std::string& prefix) {
    for (const auto symbol : symbols) {
      if (symbol < 2 || !states.emplace(symbol, prefix + std::to_string(symbol)).second) {
        throw std::runtime_error("Invalid or duplicate SEC state symbol " +
                                 std::to_string(symbol));
      }
    }
  };
  addStates(problem.state0Symbols, "design0_state_");
  addStates(problem.state1Symbols, "design1_state_");
  addStates(problem.auxiliaryStateSymbols, "auxiliary_state_");

  std::map<size_t, std::string> inputs;
  for (size_t index = 0; index < problem.inputSymbols.size(); ++index) {
    const auto symbol = problem.inputSymbols[index];
    if (symbol < 2 || states.contains(symbol)) {
      throw std::runtime_error("Invalid SEC input symbol " + std::to_string(symbol));
    }
    const auto name = index < problem.environmentInputNames.size()
        ? problem.environmentInputNames[index] : std::string();
    inputs.emplace(symbol, "input_" + std::to_string(symbol) +
        (name.empty() ? "" : "_" + name));
  }
  // allSymbols also includes explicitly published environment frontiers.
  for (const auto symbol : problem.allSymbols) {
    if (symbol >= 2 && !states.contains(symbol)) {
      inputs.try_emplace(symbol, "input_" + std::to_string(symbol));
    }
  }
  std::unordered_map<size_t, NodeId> nodes;
  for (const auto& [symbol, name] : inputs) {
    const auto node = writer.input(name);
    nodes.emplace(symbol, node);
    writer.bindVariable(symbol, node);
  }
  for (const auto& [symbol, name] : states) {
    const auto node = writer.state(name);
    nodes.emplace(symbol, node);
    writer.bindVariable(symbol, node);
  }
  const auto nodeForState = [&](size_t symbol) {
    if (!states.contains(symbol)) {
      throw std::runtime_error("BTOR2 export references undeclared state " +
                               std::to_string(symbol));
    }
    return nodes.at(symbol);
  };

  TransitionExprResolver transitions(problem);
  for (const auto& [symbol, name] : states) {
    if (transitions.contains(symbol)) {
      writer.next(nodes.at(symbol), writer.expression(transitions.at(symbol)));
    } else if (const auto found = transitions.primaryByComplement().find(symbol);
               found != transitions.primaryByComplement().end() &&
               transitions.contains(found->second)) {
      writer.next(nodes.at(symbol), writer.logicalNot(
          writer.expression(transitions.at(found->second))));
    } else {
      throw std::runtime_error("BTOR2 export missing next-state expression for " + name);
    }
  }

  const auto addRelations = [&](const auto& pairs, bool complement) {
    for (const auto& [lhs, rhs] : pairs) {
      auto lhsNode = nodeForState(lhs);
      if (complement) {
        lhsNode = writer.logicalNot(lhsNode);
      }
      writer.constraint(writer.equal(lhsNode, nodeForState(rhs)));
    }
  };
  addRelations(problem.complementedStatePairs0, true);
  addRelations(problem.complementedStatePairs1, true);
  addRelations(problem.sameFrameStateEqualityPairs0, false);
  addRelations(problem.sameFrameStateEqualityPairs1, false);
  for (const auto& rails : problem.dualRailStatePairs) {
    writer.constraint(writer.logicalOr(
        nodeForState(rails.mayBeOne), nodeForState(rails.mayBeZero)));
  }

  // These rules mirror findBaseCounterexampleImpl, including the distinction
  // between a reset prefix and an observation-only uninitialized frontier.
  // Export the concrete obligation, not engine-specific inductive hypotheses.
  const size_t resetFrames = ((!problem.hasCompleteInitialState() ||
      problem.usesDualRailStateEncoding) && problem.hasResetBootstrap())
      ? problem.resetBootstrapCycles : 0;
  const bool resetObservation = resetFrames != 0 &&
      problem.usesResetBootstrapObservationFrontier();
  const bool initialObservation = resetFrames == 0 && problem.hasSequentialState() &&
      !problem.hasCompleteInitialState();
  if (resetFrames == std::numeric_limits<size_t>::max()) {
    throw std::runtime_error("BTOR2 export reset cycle count is too large");
  }
  const size_t firstBadFrame = resetFrames != 0
      ? resetFrames + (resetObservation ? 1 : 0) : (initialObservation ? 1 : 0);
  writer.comment("reset_prefix_frames=" + std::to_string(resetFrames) +
                 " first_bad_frame=" + std::to_string(firstBadFrame));
  // A bootstrap assignment is a fact at one frontier, not an invariant.
  ExportTimeline timeline(writer, resetFrames != 0
      ? std::max(firstBadFrame, resetFrames + 1) : firstBadFrame);

  const bool initialize = resetFrames != 0 ? problem.usesDualRailStateEncoding
      : problem.hasSequentialState() && problem.hasExplicitInitialState();
  if (initialize && (resetFrames != 0 || problem.initialCondition != nullptr)) {
    if (!problem.initialStateAssignments.empty()) {
      std::map<size_t, bool> values;
      for (const auto& [symbol, value] : problem.initialStateAssignments) {
        const auto [entry, inserted] = values.emplace(symbol, value);
        if (!inserted && entry->second != value) {
          throw std::runtime_error("Conflicting SEC initial assignments");
        }
      }
      for (const auto& [symbol, value] : values) {
        writer.init(nodeForState(symbol), writer.constant(value));
      }
    } else if (resetFrames == 0 && problem.initialCondition != BoolExpr::createTrue()) {
      constrainWhen(writer, timeline.firstFrame(), writer.expression(problem.initialCondition));
    }
  }
  const auto property = writer.expression(problem.property);
  if (initialObservation) {
    constrainWhen(writer, timeline.firstFrame(), property);
  }
  if (resetFrames != 0) {
    const auto active = timeline.before(resetFrames);
    for (const auto& [symbol, asserted] : problem.resetBootstrapInputs) {
      if (!inputs.contains(symbol)) {
        throw std::runtime_error("BTOR2 export reset port is not an input");
      }
      writer.constraint(writer.equal(nodes.at(symbol),
          asserted ? active : writer.logicalNot(active)));
    }
    const auto frontier = timeline.at(resetFrames);
    if (resetObservation) {
      constrainWhen(writer, frontier, property);
    }
    if (!problem.usesDualRailStateEncoding) {
      for (const auto& [symbol, value] : problem.bootstrapStateAssignments) {
        const auto node = nodeForState(symbol);
        constrainWhen(writer, frontier, value ? node : writer.logicalNot(node));
      }
    }
  }
  auto bad = writer.expression(problem.bad);
  if (firstBadFrame != 0) {
    bad = writer.logicalAnd(writer.logicalNot(timeline.before(firstBadFrame)), bad);
  }
  writer.bad(bad, "equivalence_mismatch");
}

void exportSecBtor2File(const KInductionProblem& problem,
                       const std::string& path,
                       const SecBtor2Metadata& metadata) {
  if (path.empty()) {
    throw std::runtime_error("BTOR2 export path must not be empty");
  }
  const std::filesystem::path destination(path);
  TemporaryExport temporary(destination);
  std::ofstream output(temporary.path(), std::ios::binary);
  output.exceptions(std::ios::badbit | std::ios::failbit);
  exportSecBtor2(problem, output, metadata);
  output.close();
  std::filesystem::rename(temporary.path(), destination);
}

}  // namespace KEPLER_FORMAL::SEC
