// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "BoolExpr.h"
#include "BoolExprCache.h"
#include "export/SecBtor2Exporter.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {
namespace {

class SecBtor2ExporterTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
};

// An independent exhaustive evaluator for the small bit-vector transition
// systems used below. Tests inspect reachable behavior, without depending on
// exporter node IDs, expression shape, or SAT encodings.
class TinyBtor2Model {
 public:
  explicit TinyBtor2Model(const std::string& source) {
    std::istringstream input(source);
    std::string line;
    size_t previous = 0;
    while (std::getline(input, line)) {
      line = line.substr(0, line.find(';'));
      std::istringstream fields(line);
      Node node;
      if (!(fields >> node.id)) {
        continue;
      }
      if (!(fields >> node.op) || node.id <= previous) {
        throw std::runtime_error("BTOR2 IDs must increase");
      }
      previous = node.id;
      std::string token;
      while (fields >> token) {
        node.args.push_back(token);
      }
      if (node.op == "sort") {
        if (node.args.size() != 2 || node.args[0] != "bitvec") {
          throw std::runtime_error("Expected bit-vector sort");
        }
        sorts_.emplace(node.id, std::stoul(node.args[1]));
      } else if (node.op == "bad" || node.op == "constraint" ||
                 node.op == "output") {
        reference(node.args.at(0));
      } else {
        node.width = sorts_.at(std::stoul(node.args.at(0)));
        if (node.width == 0 || node.width >= 63) {
          throw std::runtime_error("Tiny evaluator requires widths 1..62");
        }
        if (node.op == "state") {
          states_.push_back(node.id);
          stateOffsets_[node.id] = stateBits_;
          stateBits_ += node.width;
        } else if (node.op == "input") {
          inputs_.push_back(node.id);
          inputOffsets_[node.id] = inputBits_;
          inputBits_ += node.width;
        } else if (node.op == "init" || node.op == "next") {
          const size_t state = reference(node.args.at(1));
          reference(node.args.at(2));
          if (nodes_.at(state).op != "state") {
            throw std::runtime_error("Expected state assignment");
          }
          auto& assignments = node.op == "init" ? inits_ : nexts_;
          if (!assignments.emplace(state, std::stoul(node.args[2])).second) {
            throw std::runtime_error("Duplicate state assignment");
          }
        } else if (node.op != "zero" && node.op != "one" &&
                   node.op != "ones" && node.op != "const" &&
                   node.op != "constd" && node.op != "consth") {
          const size_t operands =
              node.op == "ite" ? 3 :
              node.op == "not" || node.op == "inc" || node.op == "dec" ||
                      node.op == "uext" || node.op == "slice" ? 1 : 2;
          for (size_t index = 1; index <= operands; ++index) {
            reference(node.args.at(index));
          }
        }
      }
      order_.push_back(node.id);
      nodes_.emplace(node.id, std::move(node));
    }
    if (stateBits_ > 14 || inputBits_ > 8) {
      throw std::runtime_error("Test model is too large for exhaustive evaluation");
    }
    for (size_t state : states_) {
      if (nexts_.count(state) == 0) {
        throw std::runtime_error("Exported state has no next definition");
      }
    }
  }

  struct Reachability {
    std::vector<bool> bad;
    std::vector<bool> legal;
  };

  Reachability explore(size_t frames) const {
    Reachability result{std::vector<bool>(frames), std::vector<bool>(frames)};
    std::set<uint64_t> reachable;
    for (uint64_t state = 0; state < (uint64_t{1} << stateBits_); ++state) {
      reachable.insert(state);
    }
    for (size_t frame = 0; frame < frames; ++frame) {
      std::set<uint64_t> successors;
      for (uint64_t state : reachable) {
        for (uint64_t input = 0; input < (uint64_t{1} << inputBits_); ++input) {
          const auto values = evaluate(state, input);
          bool legal = true;
          if (frame == 0) {
            for (const auto& [target, expression] : inits_) {
              legal &= values.at(target) == values.at(expression);
            }
          }
          bool bad = false;
          for (size_t id : order_) {
            const auto& node = nodes_.at(id);
            if (node.op == "constraint") {
              legal &= values.at(std::stoul(node.args[0])) != 0;
            } else if (node.op == "bad") {
              bad |= values.at(std::stoul(node.args[0])) != 0;
            }
          }
          if (!legal) {
            continue;
          }
          result.legal[frame] = true;
          result.bad[frame] = result.bad[frame] || bad;
          uint64_t next = 0;
          for (const auto& [target, expression] : nexts_) {
            next |= values.at(expression) << stateOffsets_.at(target);
          }
          successors.insert(next);
        }
      }
      reachable = std::move(successors);
    }
    return result;
  }

 private:
  struct Node {
    size_t id = 0;
    std::string op;
    size_t width = 0;
    std::vector<std::string> args;
  };

  size_t reference(const std::string& token) const {
    const size_t id = std::stoul(token);
    if (nodes_.count(id) == 0 || nodes_.at(id).op == "sort") {
      throw std::runtime_error("BTOR2 reference is not an earlier expression");
    }
    return id;
  }

  static uint64_t mask(size_t width) {
    return (uint64_t{1} << width) - 1;
  }

  std::unordered_map<size_t, uint64_t> evaluate(
      uint64_t state, uint64_t input) const {
    std::unordered_map<size_t, uint64_t> values;
    for (size_t id : order_) {
      const auto& n = nodes_.at(id);
      if (n.op == "sort" || n.op == "init" || n.op == "next" ||
          n.op == "bad" || n.op == "constraint" || n.op == "output") {
        continue;
      }
      const auto arg = [&](size_t index) {
        return values.at(std::stoul(n.args.at(index)));
      };
      uint64_t value = 0;
      if (n.op == "state") value = state >> stateOffsets_.at(id);
      else if (n.op == "input") value = input >> inputOffsets_.at(id);
      else if (n.op == "zero") value = 0;
      else if (n.op == "one") value = 1;
      else if (n.op == "ones") value = mask(n.width);
      else if (n.op == "const") value = std::stoull(n.args[1], nullptr, 2);
      else if (n.op == "constd") value = std::stoull(n.args[1], nullptr, 10);
      else if (n.op == "consth") value = std::stoull(n.args[1], nullptr, 16);
      else if (n.op == "not") value = ~arg(1);
      else if (n.op == "and") value = arg(1) & arg(2);
      else if (n.op == "or") value = arg(1) | arg(2);
      else if (n.op == "xor") value = arg(1) ^ arg(2);
      else if (n.op == "xnor") value = ~(arg(1) ^ arg(2));
      else if (n.op == "eq" || n.op == "iff") value = arg(1) == arg(2);
      else if (n.op == "neq") value = arg(1) != arg(2);
      else if (n.op == "implies") value = !arg(1) || arg(2);
      else if (n.op == "ite") value = arg(1) ? arg(2) : arg(3);
      else if (n.op == "add") value = arg(1) + arg(2);
      else if (n.op == "sub") value = arg(1) - arg(2);
      else if (n.op == "inc") value = arg(1) + 1;
      else if (n.op == "dec") value = arg(1) - 1;
      else if (n.op == "ult") value = arg(1) < arg(2);
      else if (n.op == "ulte") value = arg(1) <= arg(2);
      else if (n.op == "uext") value = arg(1);
      else if (n.op == "slice") value = arg(1) >> std::stoul(n.args[3]);
      else if (n.op == "concat") {
        value = (arg(1) << nodes_.at(std::stoul(n.args[2])).width) | arg(2);
      } else {
        throw std::runtime_error("Unsupported tiny-evaluator operator: " + n.op);
      }
      values.emplace(id, value & mask(n.width));
    }
    return values;
  }

  std::map<size_t, Node> nodes_;
  std::vector<size_t> order_;
  std::unordered_map<size_t, size_t> sorts_;
  std::vector<size_t> states_;
  std::vector<size_t> inputs_;
  std::unordered_map<size_t, size_t> stateOffsets_;
  std::unordered_map<size_t, size_t> inputOffsets_;
  std::map<size_t, size_t> inits_;
  std::map<size_t, size_t> nexts_;
  size_t stateBits_ = 0;
  size_t inputBits_ = 0;
};

std::string dump(const KInductionProblem& problem) {
  std::ostringstream output;
  exportSecBtor2(problem, output);
  return output.str();
}

void setBad(KInductionProblem& problem, BoolExpr* bad) {
  problem.bad = bad;
  problem.property = BoolExpr::Not(bad);
  problem.observedOutputExprs0 = {bad};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"difference"};
}

KInductionProblem stateProblem(BoolExpr* next, bool initialized = true) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0 = {{2, next}};
  problem.totalStateCount = 1;
  if (initialized) {
    problem.initializedStateCount = 1;
    problem.initialStateAssignments = {{2, false}};
    problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  }
  setBad(problem, BoolExpr::Var(2));
  return problem;
}

TEST_F(SecBtor2ExporterTests, CombinationalIdenticalAndDifferentOutputs) {
  KInductionProblem problem;
  problem.inputSymbols = {2};
  problem.allSymbols = {2};
  problem.environmentInputNames = {"input"};
  setBad(problem, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(2)));
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(1).bad,
            std::vector<bool>({false}));
  setBad(problem, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::createFalse()));
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(1).bad,
            std::vector<bool>({true}));
}

TEST_F(SecBtor2ExporterTests, ExplicitInitializationAllowsFrameZeroViolation) {
  auto problem = stateProblem(BoolExpr::createFalse());
  problem.initialStateAssignments = {{2, true}};
  problem.initialCondition = BoolExpr::Var(2);
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(3).bad,
            std::vector<bool>({true, false, false}));
}

TEST_F(SecBtor2ExporterTests, EagerTransitionsPreserveSequentialDivergence) {
  auto problem = stateProblem(BoolExpr::Not(BoolExpr::Var(2)));
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(4).bad,
            std::vector<bool>({false, true, false, true}));
}

TEST_F(SecBtor2ExporterTests, UninitializedObservationAssumesEqualityAtFrameZero) {
  auto problem = stateProblem(BoolExpr::Var(2), false);
  const auto behavior = TinyBtor2Model(dump(problem)).explore(3);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true}));
  problem.transitions0 = {{2, BoolExpr::createTrue()}};
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(3).bad,
            std::vector<bool>({false, true, true}));
}

TEST_F(SecBtor2ExporterTests, PartialInitializationPreservesKnownBitAndObservation) {
  auto problem = stateProblem(BoolExpr::Var(2));
  problem.state1Symbols = {3};
  problem.allSymbols.push_back(3);
  problem.transitions1 = {{3, BoolExpr::Var(3)}};
  problem.totalStateCount = 2;
  setBad(problem, BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3)));
  const auto behavior = TinyBtor2Model(dump(problem)).explore(3);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true}));
}

TEST_F(SecBtor2ExporterTests, LazyTransitionsUseCombinedSymbolRemapping) {
  auto problem = stateProblem(BoolExpr::createFalse());
  problem.transitions0.clear();
  problem.inputSymbols = {3};
  problem.allSymbols.push_back(3);
  problem.environmentInputNames = {"data"};
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  problem.lazyTransitions->sourceByStateSymbol.emplace(
      2, LazyTransitionSource{0, BoolExpr::Var(20), LazyTransitionRail::Binary});
  problem.lazyTransitions->localToCombinedByDesign[0].emplace(20, 3);
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(3).bad,
            std::vector<bool>({false, true, true}));
}

TEST_F(SecBtor2ExporterTests, ComplementedStateWithoutIndependentNextIsPreserved) {
  auto problem = stateProblem(BoolExpr::Not(BoolExpr::Var(2)));
  problem.state0Symbols.push_back(3);
  problem.allSymbols.push_back(3);
  problem.complementedStatePairs0 = {{2, 3}};
  problem.initialStateAssignments.push_back({3, true});
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Var(3));
  problem.initializedStateCount = problem.totalStateCount = 2;
  setBad(problem, BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3))));
  const auto behavior = TinyBtor2Model(dump(problem)).explore(4);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, false, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true, true}));
}

TEST_F(SecBtor2ExporterTests, SameFrameEqualityConstrainsBothDesignStates) {
  auto problem = stateProblem(BoolExpr::Var(2));
  problem.state1Symbols = {3};
  problem.allSymbols.push_back(3);
  problem.transitions1 = {{3, BoolExpr::Var(3)}};
  problem.sameFrameStateEqualityPairs1 = {{2, 3}};
  problem.initialStateAssignments.clear();
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = problem.totalStateCount = 2;
  setBad(problem, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3)));
  const auto behavior = TinyBtor2Model(dump(problem)).explore(3);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true}));
}

TEST_F(SecBtor2ExporterTests, DualRailValidityExcludesEmptyRailValue) {
  auto problem = stateProblem(BoolExpr::Var(2));
  problem.state0Symbols.push_back(3);
  problem.allSymbols.push_back(3);
  problem.transitions0.push_back({3, BoolExpr::Var(3)});
  problem.usesDualRailStateEncoding = true;
  problem.dualRailStatePairs = {{2, 3}};
  problem.initialStateAssignments.clear();
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = problem.totalStateCount = 2;
  setBad(problem, BoolExpr::Not(BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3))));
  const auto behavior = TinyBtor2Model(dump(problem)).explore(3);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true}));
}

TEST_F(SecBtor2ExporterTests, ResetPrefixIsForcedThenPermanentlyDeasserted) {
  auto problem = stateProblem(BoolExpr::Var(3), false);
  problem.inputSymbols = {3};
  problem.allSymbols.push_back(3);
  problem.environmentInputNames = {"reset"};
  problem.resetBootstrapInputs = {{3, true}};
  problem.resetBootstrapCycles = 2;
  problem.bootstrapStateAssignments = {{2, true}};
  auto behavior = TinyBtor2Model(dump(problem)).explore(5);
  EXPECT_EQ(behavior.bad,
            std::vector<bool>({false, false, true, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true, true, true}));
  problem.resetBootstrapInputs = {{3, false}};
  problem.transitions0 = {{2, BoolExpr::Not(BoolExpr::Var(3))}};
  behavior = TinyBtor2Model(dump(problem)).explore(5);
  EXPECT_EQ(behavior.bad,
            std::vector<bool>({false, false, true, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true, true, true}));
}

TEST_F(SecBtor2ExporterTests, IncompleteResetUsesSafeObservationFrontier) {
  auto problem = stateProblem(BoolExpr::Not(BoolExpr::Var(3)), false);
  problem.inputSymbols = {3};
  problem.allSymbols.push_back(3);
  problem.resetBootstrapInputs = {{3, true}};
  problem.resetBootstrapCycles = 1;
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(4).bad,
            std::vector<bool>({false, false, true, true}));
  // If the observation itself is unsafe, it must exclude that trace rather
  // than merely hide its first bad output.
  problem.transitions0 = {{2, BoolExpr::createTrue()}};
  const auto behavior = TinyBtor2Model(dump(problem)).explore(4);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, false, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, false, false, false}));
}

TEST_F(SecBtor2ExporterTests, DualRailResetChecksFirstPostResetFrame) {
  auto problem = stateProblem(BoolExpr::Var(3));
  problem.state0Symbols.push_back(4);
  problem.inputSymbols = {3};
  problem.allSymbols = {2, 3, 4};
  problem.transitions0.push_back({4, BoolExpr::Not(BoolExpr::Var(3))});
  problem.usesDualRailStateEncoding = true;
  problem.dualRailStatePairs = {{2, 4}};
  problem.initialStateAssignments = {{2, true}, {4, true}};
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = problem.totalStateCount = 2;
  problem.resetBootstrapInputs = {{3, true}};
  problem.resetBootstrapCycles = 1;
  setBad(problem, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Not(BoolExpr::Var(4))));
  const auto behavior = TinyBtor2Model(dump(problem)).explore(4);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, true, false, false}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true, true}));
  // The concrete dual-rail reset prefix consumes structured initial facts even
  // when the optional initial-condition formula is absent.
  problem.initialCondition = nullptr;
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(4).bad, behavior.bad);
}

TEST_F(SecBtor2ExporterTests, CompleteBinaryInitializationBypassesResetBootstrap) {
  auto problem = stateProblem(BoolExpr::Var(3));
  problem.inputSymbols = {3};
  problem.allSymbols.push_back(3);
  problem.resetBootstrapInputs = {{3, true}};
  problem.resetBootstrapCycles = 2;
  setBad(problem, BoolExpr::Not(BoolExpr::Var(3)));
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(3).bad,
            std::vector<bool>({true, true, true}));
}

TEST_F(SecBtor2ExporterTests, RelationalInitialConditionOnlyConstrainsInitialFrame) {
  auto problem = stateProblem(BoolExpr::createTrue());
  problem.state1Symbols = {3};
  problem.allSymbols.push_back(3);
  problem.transitions1 = {{3, BoolExpr::createTrue()}};
  problem.initialStateAssignments.clear();
  problem.initialCondition = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
  problem.initializedStateCount = problem.totalStateCount = 2;
  setBad(problem, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  const auto behavior = TinyBtor2Model(dump(problem)).explore(3);
  EXPECT_EQ(behavior.bad, std::vector<bool>({false, true, true}));
  EXPECT_EQ(behavior.legal, std::vector<bool>({true, true, true}));
}

TEST_F(SecBtor2ExporterTests, AuxiliaryStateTransitionIsIncluded) {
  auto problem = stateProblem(BoolExpr::createFalse());
  problem.auxiliaryStateSymbols = {3};
  problem.auxiliaryTransitions = {{3, BoolExpr::createTrue()}};
  problem.allSymbols.push_back(3);
  problem.initialStateAssignments.push_back({3, false});
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3)));
  problem.initializedStateCount = problem.totalStateCount = 2;
  setBad(problem, BoolExpr::Var(3));
  EXPECT_EQ(TinyBtor2Model(dump(problem)).explore(3).bad,
            std::vector<bool>({false, true, true}));
}

TEST_F(SecBtor2ExporterTests, RejectsUndeclaredSymbolsAndMissingStateTransitions) {
  auto problem = stateProblem(BoolExpr::Var(999));
  EXPECT_THROW(dump(problem), std::exception);
  problem.transitions0.clear();
  EXPECT_THROW(dump(problem), std::exception);
  problem.transitions0 = {{2, BoolExpr::Var(2)}};
  setBad(problem, BoolExpr::Var(999));
  EXPECT_THROW(dump(problem), std::exception);
  setBad(problem, BoolExpr::Var(2));
  problem.state1Symbols = {2};
  EXPECT_THROW(dump(problem), std::exception);
}

TEST_F(SecBtor2ExporterTests, RepeatedLazyExportIsByteForByteDeterministic) {
  auto problem = stateProblem(BoolExpr::createFalse());
  problem.transitions0.clear();
  problem.inputSymbols = {3};
  problem.allSymbols.push_back(3);
  problem.environmentInputNames = {"data"};
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  problem.lazyTransitions->sourceByStateSymbol.emplace(
      2, LazyTransitionSource{0, BoolExpr::Var(20), LazyTransitionRail::Binary});
  problem.lazyTransitions->localToCombinedByDesign[0].emplace(20, 3);
  const auto first = dump(problem);
  EXPECT_EQ(first, dump(problem));
  EXPECT_EQ(first, dump(problem));
}

TEST_F(SecBtor2ExporterTests, FileReplacementPreservesDestinationOnFailure) {
  struct TemporaryDirectory {
    std::filesystem::path path;
    TemporaryDirectory() {
      std::random_device random;
      for (unsigned attempt = 0; attempt < 32; ++attempt) {
        path = std::filesystem::temp_directory_path() /
            ("kf-btor2-test-" + std::to_string(random()) + "-" +
             std::to_string(random()));
        if (std::filesystem::create_directory(path)) {
          return;
        }
      }
      throw std::runtime_error("Cannot create test directory");
    }
    ~TemporaryDirectory() {
      std::error_code error;
      std::filesystem::remove_all(path, error);
    }
  } directory;
  const auto destination = directory.path / "model.btor2";
  {
    std::ofstream original(destination);
    original << "previous export\n";
  }
  const auto readDestination = [&]() {
    std::ifstream file(destination);
    std::ostringstream contents;
    contents << file.rdbuf();
    return contents.str();
  };

  auto invalid = stateProblem(BoolExpr::Var(999));
  EXPECT_THROW(exportSecBtor2File(invalid, destination.string()), std::exception);
  EXPECT_EQ(readDestination(), "previous export\n");
  EXPECT_EQ(std::distance(std::filesystem::directory_iterator(directory.path),
                          std::filesystem::directory_iterator()), 1);

  auto valid = stateProblem(BoolExpr::Not(BoolExpr::Var(2)));
  EXPECT_NO_THROW(exportSecBtor2File(valid, destination.string()));
  EXPECT_EQ(readDestination(), dump(valid));
  EXPECT_EQ(TinyBtor2Model(readDestination()).explore(3).bad,
            std::vector<bool>({false, true, false}));
  EXPECT_EQ(std::distance(std::filesystem::directory_iterator(directory.path),
                          std::filesystem::directory_iterator()), 1);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
