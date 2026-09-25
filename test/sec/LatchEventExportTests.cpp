// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB.h"
#include "NLDB0.h"
#include "NLLibrary.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "latch/LatchSupportOptions.h"
#include "latch/LatchResetAdapter.h"
#include "model/SequentialDesignModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using namespace naja::NL;

struct Event {
  size_t selector;
  bool value;
  bool sampleData = false;
  bool sampleEnable = false;
  size_t order = 0;
};

bool eventInput(const std::string& name, Event event) {
  if (name.find("$event.value") != std::string::npos) return event.value;
  const std::string orderPrefix = "$event.reset.order[";
  if (const auto position = name.find(orderPrefix); position != std::string::npos) {
    const auto stage = std::stoul(name.substr(position + orderPrefix.size()));
    return (event.order >> stage) & 1;  // These circuits have exactly two sampled pins.
  }
  const std::string prefix = "$event.select[";
  if (const auto position = name.find(prefix); position != std::string::npos) {
    return (event.selector >> std::stoul(name.substr(position + prefix.size()))) & 1;
  }
  // Original PI sentinels supply desired levels during reset, unused afterward.
  if (name.ends_with("data[0]")) return event.sampleData;
  if (name.ends_with("enable[0]")) return event.sampleEnable;
  if (name.ends_with("clock[0]") || name.ends_with("reset[0]")) return false;
  throw std::runtime_error("Unexpected event input: " + name);
}

// A trace interpreter, independent of BoolExpr and the SEC transition builder.
// It interprets the one-bit BTOR2 operations emitted by Btor2Writer, including
// simultaneous state updates and constraints. No state-space-size limit or
// assumption about exported state numbering is needed for these concrete traces.
class BtorTrace {
 public:
  explicit BtorTrace(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot read exported BTOR2");
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream fields(line);
      size_t id = 0;
      if (!(fields >> id)) continue;
      Node node;
      fields >> node.operation;
      std::string field;
      while (fields >> field) node.arguments.push_back(field);
      if (!nodes_.emplace(id, std::move(node)).second)
        throw std::runtime_error("Duplicate BTOR2 node");
    }
    for (const auto& [id, node] : nodes_) {
      if (node.operation == "state") state_[id] = false;
      else if (node.operation == "init")
        initial_.emplace(number(node, 1), number(node, 2));
      else if (node.operation == "next")
        next_.emplace(number(node, 1), number(node, 2));
    }
    if (initial_.size() != state_.size() || next_.size() != state_.size())
      throw std::runtime_error("Event export must initialize and update every state");
    const auto values = evaluate({0, false});
    for (auto& [id, value] : state_) value = values.at(initial_.at(id));
  }

  bool step(Event event) {
    const auto values = evaluate(event);
    bool bad = false;
    size_t badCount = 0;
    for (const auto& [id, node] : nodes_) {
      if (node.operation == "constraint" && !values.at(number(node, 0)))
        throw std::runtime_error("Export excluded an allowed concrete event trace");
      if (node.operation == "bad") {
        bad |= values.at(number(node, 0));
        ++badCount;
      }
    }
    if (badCount != 1) throw std::runtime_error("Expected one SEC bad property");
    for (auto& [id, value] : state_) value = values.at(next_.at(id));
    return bad;
  }

 private:
  struct Node {
    std::string operation;
    std::vector<std::string> arguments;
  };
  static size_t number(const Node& node, size_t index) {
    return std::stoull(node.arguments.at(index));
  }
  std::map<size_t, bool> evaluate(Event event) const {
    std::map<size_t, bool> values;
    for (const auto& [id, node] : nodes_) {
      const auto& op = node.operation;
      const auto operand = [&](size_t index) { return values.at(number(node, index)); };
      if (op == "sort") {
        if (node.arguments != std::vector<std::string>{"bitvec", "1"})
          throw std::runtime_error("Expected one-bit BTOR2 sort");
      } else if (op == "input") values[id] = eventInput(node.arguments.at(1), event);
      else if (op == "state") values[id] = state_.at(id);
      else if (op == "const") {
        if (node.arguments.at(1) != "0" && node.arguments.at(1) != "1")
          throw std::runtime_error("Expected Boolean BTOR2 constant");
        values[id] = node.arguments.at(1) == "1";
      } else if (op == "not") values[id] = !operand(1);
      else if (op == "and") values[id] = operand(1) && operand(2);
      else if (op == "or") values[id] = operand(1) || operand(2);
      else if (op == "xor") values[id] = operand(1) != operand(2);
      else if (op == "eq") values[id] = operand(1) == operand(2);
      else if (op != "init" && op != "next" && op != "constraint" &&
               op != "bad" && op != "output")
        throw std::runtime_error("Unhandled BTOR2 operation: " + op);
    }
    return values;
  }
  std::map<size_t, Node> nodes_;
  std::map<size_t, bool> state_;
  std::map<size_t, size_t> initial_, next_;
};

struct TemporaryDirectory {
  std::filesystem::path path;
  TemporaryDirectory() {
    for (size_t attempt = 0; attempt < 32; ++attempt) {
      path = std::filesystem::temp_directory_path() /
          ("kf-latch-export-" + std::to_string(std::random_device{}()) + "-" + std::to_string(attempt));
      if (std::filesystem::create_directory(path)) return;
    }
    throw std::runtime_error("Cannot create latch export test directory");
  }
  ~TemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
};

using State = std::unordered_map<size_t, bool>;
State initialState(const SequentialDesignModel& model) {
  State result;
  for (const auto& key : model.stateBits)
    result[model.inputVarByKey.at(key)] = model.initialStateValueByKey.at(key);
  return result;
}

bool nativeStep(const SequentialDesignModel& model, State& state, Event event) {
  auto values = state;
  for (const auto& key : model.environmentInputs)
    values[model.inputVarByKey.at(key)] = eventInput(model.displayNameByKey.at(key), event);
  const bool output = model.observedOutputExprByKey.at(model.observedOutputs.at(0))->evaluate(values);
  for (const auto& key : model.stateBits)
    state[model.inputVarByKey.at(key)] = model.nextStateExprByStateKey.at(key)->evaluate(values);
  return output;
}

class LatchEventExportTests : public ::testing::Test {
 protected:
  void SetUp() override {
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    designs_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    primitives_ = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
  }
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    BoolExprCache::destroy();
  }

  SNLDesign* chain(const char* name, size_t length, bool resetFlop = false) {
    auto* top = SNLDesign::create(designs_, SNLDesign::Type::Standard, NLName(name));
    const auto input = [&](const char* portName) {
      auto* term = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName(portName));
      auto* net = SNLScalarNet::create(top, NLName(portName));
      term->setNet(net);
      return net;
    };
    auto* data = input("data");
    auto* enable = input("enable");
    if (resetFlop) {
      auto* clock = input("clock");
      auto* reset = input("reset");
      auto* primitive = SNLDesign::create(primitives_, SNLDesign::Type::Primitive, NLName(std::string(name) + "_ff"));
      auto* d = SNLScalarTerm::create(primitive, SNLTerm::Direction::Input, NLName("D"));
      auto* c = SNLScalarTerm::create(primitive, SNLTerm::Direction::Input, NLName("C"));
      auto* r = SNLScalarTerm::create(primitive, SNLTerm::Direction::Input, NLName("R"));
      auto* q = SNLScalarTerm::create(primitive, SNLTerm::Direction::Output, NLName("Q"));
      using Modeling = SNLDesignModeling;
      using Operator = Modeling::BooleanExpression::Operator;
      Modeling::SequentialModel sequential;
      sequential.kind = Modeling::SequentialModel::Kind::FlipFlop;
      sequential.clockedOn.root = sequential.clockedOn.addTerm(c);
      Modeling::SequentialState storage;
      auto& expression = storage.nextState;
      const auto dNode = expression.addTerm(d);
      const auto notReset = expression.addOperation(Operator::Not, {expression.addTerm(r)});
      expression.root = expression.addOperation(Operator::And, {dNode, notReset});
      sequential.states.push_back(storage);
      Modeling::BooleanExpression output;
      output.root = output.addState(0);
      sequential.outputs.push_back({q, output});
      Modeling::setSequentialModel(primitive, sequential);
      auto* ff = SNLInstance::create(top, primitive, NLName("reset_ff"));
      auto* captured = SNLScalarNet::create(top, NLName("captured"));
      ff->getInstTerm(d)->setNet(data);
      ff->getInstTerm(c)->setNet(clock);
      ff->getInstTerm(r)->setNet(reset);
      ff->getInstTerm(q)->setNet(captured);
      data = captured;
    }
    for (size_t index = 0; index < length; ++index) {
      auto* latch = SNLInstance::create(top, NLDB0::getDLatch(), NLName("latch" + std::to_string(index)));
      auto* output = SNLScalarNet::create(top, NLName("q" + std::to_string(index)));
      latch->getInstTerm(NLDB0::getDLatchData())->setNet(data);
      latch->getInstTerm(NLDB0::getDLatchEnable())->setNet(enable);
      latch->getInstTerm(NLDB0::getDLatchOutput())->setNet(output);
      data = output;
    }
    SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"))->setNet(data);
    return top;
  }

  void compareResetTrace(bool rightIsWire) {
    auto* leftTop = chain("reset_single", 1, true);
    auto* rightTop = chain("reset_reference", rightIsWire ? 0 : 4, true);
    SupportOptions options;
    options.enabled = true;
    options.singleInputChange = true;
    options.initialInputs = false;
    options.initialStorage = true;
    options.workers = 2;
    ScopedSupportOptions scope(options);
    const auto left = SequentialDesignModel::extract(leftTop);
    const auto right = SequentialDesignModel::extract(rightTop);
    ASSERT_FALSE(left.hasUnsupportedFeatures());
    ASSERT_FALSE(right.hasUnsupportedFeatures());
    ASSERT_EQ(left.observedOutputs.size(), 1u);
    ASSERT_EQ(right.observedOutputs.size(), 1u);
    const SecResetSpec reset{2, {{"reset", true}}};
    const auto adaptedLeft = adaptResetCycles(left, reset);
    const auto adaptedRight = adaptResetCycles(right, reset);
    ASSERT_TRUE(adaptedLeft.model) << adaptedLeft.error;
    ASSERT_TRUE(adaptedRight.model) << adaptedRight.error;
    TemporaryDirectory directory;
    for (const auto engine : {SecEngine::KInduction, SecEngine::Imc, SecEngine::Pdr}) {
      for (const auto encoding : {SecEncoding::Binary, SecEncoding::DualRailSteady}) {
        SCOPED_TRACE(::testing::Message() << "engine=" << static_cast<int>(engine)
            << " encoding=" << static_cast<int>(encoding));
        const auto path = directory.path / "reset.btor2";
        SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
            engine, encoding, reset, Btor2ExportOptions{path.string(), true});
        const auto result = strategy.runExtractedModels(left, right, 16);
        ASSERT_EQ(result.status, SequentialEquivalenceStatus::Exported) << result.reason;
        BtorTrace exported(path);
        auto leftState = initialState(*adaptedLeft.model);
        auto rightState = initialState(*adaptedRight.model);
        bool data = false, enable = false, clock = false, captured = true, held = true;
        // Sorted input positions: clock=0,data=1,enable=2,reset=3; 4 is stutter.
        const std::vector<Event> events = {
            {0, true, true, false, 0}, {3, false, false, false, 3}, // two reset clock cycles
            {4, false}, {2, true}, {1, true}, {0, true}, {2, false},
            {1, false}, {0, false}, {0, true}, {3, true}, {4, false}, {2, true}};
        for (size_t frame = 0; frame < events.size(); ++frame) {
          SCOPED_TRACE(::testing::Message() << "reset/event step=" << frame);
          const auto event = events[frame];
          if (frame < reset.cycles) {
            data = event.sampleData;
            enable = event.sampleEnable;
            captured = false;  // asserted synchronous reset sees a complete source-clock pulse
            clock = false;
            if (enable) held = captured;
          } else {
            const bool previousClock = clock;
            if (event.selector == 0) clock = event.value;
            if (event.selector == 1) data = event.value;
            if (event.selector == 2) enable = event.value;
            // A request to reassert reset is a stutter after release.
            if (!previousClock && clock) captured = data;
            if (enable) held = captured;
          }
          const auto leftOutput = nativeStep(*adaptedLeft.model, leftState, event);
          const auto rightOutput = nativeStep(*adaptedRight.model, rightState, event);
          ASSERT_EQ(leftOutput, frame < reset.cycles ? false : held);
          ASSERT_EQ(rightOutput, frame < reset.cycles ? false : rightIsWire ? captured : held);
          const auto bad = exported.step(event);
          EXPECT_EQ(bad, leftOutput != rightOutput);
          EXPECT_EQ(bad, frame >= reset.cycles && rightIsWire && held != captured);
        }
      }
    }
  }

  void compareTrace(bool rightIsWire) {
    auto* leftTop = chain("single_latch", 1);
    auto* rightTop = chain("reference", rightIsWire ? 0 : 4);
    TemporaryDirectory directory;
    for (const bool initial : {false, true}) {
      SupportOptions options;
      options.enabled = true;
      options.singleInputChange = true;
      options.initialInputs = initial;
      options.initialStorage = initial;
      options.workers = 2;
      ScopedSupportOptions scope(options);
      const auto left = SequentialDesignModel::extract(leftTop);
      const auto right = SequentialDesignModel::extract(rightTop);
      ASSERT_FALSE(left.hasUnsupportedFeatures());
      ASSERT_FALSE(right.hasUnsupportedFeatures());
      ASSERT_EQ(left.observedOutputs.size(), 1u);
      ASSERT_EQ(right.observedOutputs.size(), 1u);
      ASSERT_TRUE(left.skippedObservedOutputs.empty());
      ASSERT_TRUE(right.skippedObservedOutputs.empty());
      for (const auto engine : {SecEngine::KInduction, SecEngine::Imc, SecEngine::Pdr}) {
        for (const auto encoding : {SecEncoding::Binary, SecEncoding::DualRailSteady}) {
          SCOPED_TRACE(::testing::Message() << "initial=" << initial << " engine="
              << static_cast<int>(engine) << " encoding=" << static_cast<int>(encoding));
          const auto path = directory.path / "trace.btor2";
          SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
              engine, encoding, {}, Btor2ExportOptions{path.string(), true});
          const auto result = strategy.runExtractedModels(left, right, 16);
          ASSERT_EQ(result.status, SequentialEquivalenceStatus::Exported) << result.reason;
          BtorTrace exported(path);
          auto leftState = initialState(left);
          auto rightState = initialState(right);
          bool data = initial, enable = initial, held = initial;
          // Includes first-frame mismatch (no phantom bootstrap suppression),
          // hold, opening, transparent updates, closing, and reserved-selector stutter.
          const std::vector<Event> events = {{0, !initial}, {1, true}, {0, false},
              {0, true}, {1, false}, {0, false}, {3, true}, {1, true},
              {1, false}, {0, true}, {3, false}, {1, true}};
          size_t frame = 0;
          for (const auto event : events) {
            SCOPED_TRACE(::testing::Message() << "event=" << frame++);
            if (event.selector == 0) data = event.value;
            if (event.selector == 1) enable = event.value;
            if (enable) held = data;
            const auto leftOutput = nativeStep(left, leftState, event);
            const auto rightOutput = nativeStep(right, rightState, event);
            ASSERT_EQ(leftOutput, held);
            ASSERT_EQ(rightOutput, rightIsWire ? data : held);
            const bool exportedBad = exported.step(event);
            EXPECT_EQ(exportedBad, leftOutput != rightOutput);
            EXPECT_EQ(exportedBad, rightIsWire && held != data);
          }
        }
      }
    }
  }
 private:
  NLLibrary* designs_ = nullptr;
  NLLibrary* primitives_ = nullptr;
};

TEST_F(LatchEventExportTests, BadPropertyTracksHoldOpenAndCaptureWithoutHiddenFrames) {
  compareTrace(true);
}

TEST_F(LatchEventExportTests, DifferentInternalSettlingDepthsShareExternalObservations) {
  compareTrace(false);
}

TEST_F(LatchEventExportTests, ResetPrefixMasksOnlyResetCyclesAndRetainsPostreleaseLatchMismatch) {
  compareResetTrace(true);
}

TEST_F(LatchEventExportTests, ResetPrefixPreservesEquivalentChainsWithDifferentSettlingDepths) {
  compareResetTrace(false);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
