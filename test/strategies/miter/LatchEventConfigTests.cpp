// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <zlib.h>

#include "KeplerFormalDriver.h"
#include "LatchEventConfig.h"
#include "latch/LatchSupportOptions.h"

namespace {

using KEPLER_FORMAL::LatchEventConfig;
using Result = LatchEventConfig::ArgumentResult;

bool parseYaml(LatchEventConfig& config, const std::string& body, std::string& error) {
  return config.parseYaml(YAML::Load("latch_support: true\nsec_latch_events: " + body), error);
}

Result parseArgument(LatchEventConfig& config, std::vector<std::string> args,
                     std::string& error) {
  std::vector<char*> argv;
  for (auto& arg : args) argv.push_back(arg.data());
  int index = 0;
  return config.parseArgument(static_cast<int>(argv.size()), argv.data(), index, error);
}

constexpr auto complete = "{input_changes: single, initial_inputs: 0, initial_storage: 0}";

TEST(LatchEventConfigTests, DisabledByDefaultWithoutAnImplicitContract) {
  LatchEventConfig config;
  std::string error;
  EXPECT_TRUE(config.parseYaml(YAML::Load("{}"), error));
  EXPECT_FALSE(config.options().enabled);
  EXPECT_FALSE(config.options().initialInputs.has_value());
  EXPECT_TRUE(config.validate(false, false, false, error));
  EXPECT_EQ(parseArgument(config, {"--unrelated"}, error), Result::NotHandled);
}

TEST(LatchEventConfigTests, ExplicitMasterOffPreservesLegacyWithoutTuning) {
  LatchEventConfig config;
  std::string error;
  ASSERT_TRUE(config.parseYaml(YAML::Load("latch_support: false"), error));
  EXPECT_FALSE(config.options().enabled);
  EXPECT_TRUE(config.validate(false, true, true, error));
}

TEST(LatchEventConfigTests, EventTuningNeverEnablesTheMasterGate) {
  for (const auto* gate : {"", "latch_support: false\n"}) {
    LatchEventConfig config;
    std::string error;
    ASSERT_TRUE(config.parseYaml(YAML::Load(std::string(gate) + "sec_latch_events: " + complete), error));
    EXPECT_FALSE(config.options().enabled);
    EXPECT_FALSE(config.validate(true, false, false, error));
    EXPECT_NE(error.find("latch_support"), std::string::npos);
  }
  LatchEventConfig config;
  std::string error;
  ASSERT_EQ(parseArgument(config, {"--sec-latch-events", "single"}, error), Result::Parsed);
  EXPECT_FALSE(config.options().enabled);
  EXPECT_FALSE(config.validate(true, false, false, error));
}

TEST(LatchEventConfigTests, MasterGateRequiresExactBooleanYamlValues) {
  for (const auto* value : {"0", "1", "yes", "off", "null", "[]", "{enabled: true}"}) {
    LatchEventConfig config;
    std::string error;
    EXPECT_FALSE(config.parseYaml(YAML::Load(std::string("latch_support: ") + value), error));
    EXPECT_FALSE(config.options().enabled);
    EXPECT_NE(error.find("latch_support"), std::string::npos);
  }
}

TEST(LatchEventConfigTests, MasterEnableAloneDoesNotInventAnEventContract) {
  LatchEventConfig config;
  std::string error;
  ASSERT_EQ(parseArgument(config, {"--latch_support"}, error), Result::Parsed);
  EXPECT_TRUE(config.options().enabled);
  EXPECT_FALSE(config.validate(true, false, false, error));
  EXPECT_NE(error.find("explicit input_changes"), std::string::npos);
  ASSERT_EQ(parseArgument(config, {"--sec-latch-initial-inputs", "0"}, error), Result::Parsed);
  ASSERT_EQ(parseArgument(config, {"--sec-latch-initial-storage", "0"}, error), Result::Parsed);
  EXPECT_FALSE(config.validate(true, false, false, error));
  ASSERT_EQ(parseArgument(config, {"--sec-latch-events", "single"}, error), Result::Parsed);
  EXPECT_TRUE(config.validate(true, false, false, error));
}

TEST(LatchEventConfigTests, MasterEnableCanFollowTuningFlags) {
  LatchEventConfig config;
  std::string error;
  ASSERT_TRUE(config.parseYaml(YAML::Load(std::string("sec_latch_events: ") + complete), error));
  ASSERT_EQ(parseArgument(config, {"--latch_support"}, error), Result::Parsed);
  EXPECT_TRUE(config.validate(true, false, false, error));
}

TEST(LatchEventConfigTests, RequiresExplicitInputChangeAndBothInitializationChoices) {
  for (const auto* body : {"{}", "{input_changes: single}",
      "{input_changes: single, initial_inputs: 0}",
      "{input_changes: single, initial_storage: 0}",
      "{initial_inputs: 0, initial_storage: 0}"}) {
    SCOPED_TRACE(body);
    LatchEventConfig config;
    std::string error;
    ASSERT_TRUE(parseYaml(config, body, error));
    EXPECT_FALSE(config.validate(true, false, false, error));
    EXPECT_NE(error.find("explicit"), std::string::npos);
  }
}

TEST(LatchEventConfigTests, BooleanZeroIsPresentRatherThanMissing) {
  LatchEventConfig config;
  std::string error;
  ASSERT_TRUE(parseYaml(config, complete, error));
  EXPECT_TRUE(config.validate(true, false, false, error));
  EXPECT_TRUE(config.options().enabled);
  EXPECT_TRUE(config.options().singleInputChange);
  EXPECT_EQ(config.options().initialInputs, false);
  EXPECT_EQ(config.options().initialStorage, false);
}

TEST(LatchEventConfigTests, AnyChangesAndExplicitHighInitializationAreAccepted) {
  LatchEventConfig config;
  std::string error;
  ASSERT_TRUE(parseYaml(config,
      "{input_changes: any, initial_inputs: 1, initial_storage: 1, workers: 2, "
      "max_waves: 7, max_states: 9, max_transactions: 11, max_symbolic_nodes: 123, "
      "max_sat_conflicts: 45, max_sat_decisions: 67}", error));
  EXPECT_TRUE(config.validate(true, false, false, error));
  EXPECT_FALSE(config.options().singleInputChange);
  EXPECT_EQ(config.options().initialInputs, true);
  EXPECT_EQ(config.options().initialStorage, true);
  EXPECT_EQ(config.options().workers, 2u);
  EXPECT_EQ(config.options().limits.maxWaves, 7u);
  EXPECT_EQ(config.options().limits.maxBoundaryStates, 9u);
  EXPECT_EQ(config.options().limits.maxTransactions, 11u);
  EXPECT_EQ(config.options().maxSymbolicNodes, 123u);
  EXPECT_EQ(config.options().maxSatConflicts, 45u);
  EXPECT_EQ(config.options().maxSatDecisions, 67u);
}

TEST(LatchEventConfigTests, RejectsUnknownKeysAndMalformedYamlShapes) {
  for (const auto* body : {"null", "true", "[]", "[single]",
      "{unknown: 1}", "{workers: [1]}", "{input_changes: {mode: single}}"}) {
    SCOPED_TRACE(body);
    LatchEventConfig config;
    std::string error;
    EXPECT_FALSE(parseYaml(config, body, error));
    EXPECT_FALSE(error.empty());
  }
}

TEST(LatchEventConfigTests, RejectsBadEnumsAndNonBinaryInitialization) {
  for (const auto* argument : {"--sec-latch-events", "--sec-latch-initial-inputs",
                              "--sec-latch-initial-storage"}) {
    for (const auto* value : {"", "maybe", "true", "-1", "2"}) {
      SCOPED_TRACE(std::string(argument) + " " + value);
      LatchEventConfig config;
      std::string error;
      EXPECT_EQ(parseArgument(config, {argument, value}, error), Result::Error);
    }
  }
}

TEST(LatchEventConfigTests, RejectsNegativeOverflowAndZeroResourceLimits) {
  for (const auto* argument : {"--sec-latch-max-waves", "--sec-latch-max-states",
                              "--sec-latch-max-transactions", "--sec-latch-max-nodes",
                              "--sec-latch-sat-conflicts", "--sec-latch-sat-decisions"}) {
    for (const auto* value : {"0", "-1", "184467440737095516160", "1.2", "3junk"}) {
      LatchEventConfig config;
      std::string error;
      EXPECT_EQ(parseArgument(config, {argument, value}, error), Result::Error);
    }
  }
  LatchEventConfig config;
  std::string error;
  EXPECT_EQ(parseArgument(config, {"--sec-latch-workers", "-1"}, error), Result::Error);
  EXPECT_EQ(parseArgument(config, {"--sec-latch-workers", "2147483648"}, error), Result::Error);
  EXPECT_EQ(parseArgument(config, {"--sec-latch-workers", "0"}, error), Result::Parsed);
}

TEST(LatchEventConfigTests, SymbolicBudgetsAreCheckedAndDoNotEnableGate) {
  for (const auto* flag : {"--sec-latch-sat-conflicts", "--sec-latch-sat-decisions"}) {
    LatchEventConfig config;
    std::string error;
    EXPECT_EQ(parseArgument(config, {flag, "4294967296"}, error), Result::Error);
    EXPECT_EQ(parseArgument(config, {flag, "42"}, error), Result::Parsed);
    EXPECT_FALSE(config.options().enabled);
    EXPECT_FALSE(config.validate(true, false, false, error));
  }
  for (const auto* key : {"max_symbolic_nodes", "max_sat_conflicts", "max_sat_decisions"}) {
    LatchEventConfig config;
    std::string error;
    EXPECT_FALSE(parseYaml(config, std::string("{") + key + ": 0}", error));
  }
  LatchEventConfig config;
  std::string error;
  ASSERT_EQ(parseArgument(config, {"--sec-latch-max-nodes", "99"}, error), Result::Parsed);
  EXPECT_EQ(config.options().maxSymbolicNodes, 99u);
  EXPECT_FALSE(config.validate(true, false, false, error));
}

TEST(LatchEventConfigTests, RejectsMissingFlagValueAndIncompatibleWorkflows) {
  LatchEventConfig config;
  std::string error;
  EXPECT_EQ(parseArgument(config, {"--sec-latch-events"}, error), Result::Error);
  ASSERT_TRUE(parseYaml(config, complete, error));
  EXPECT_FALSE(config.validate(false, false, false, error));
  EXPECT_TRUE(config.validate(true, true, false, error));
  EXPECT_FALSE(config.validate(true, false, true, error));
  EXPECT_NE(error.find("leaf boundaries"), std::string::npos);
}

class LatchEventCliTests : public ::testing::Test {
 protected:
  void SetUp() override {
    oldDirectory_ = std::filesystem::current_path();
    directory_ = std::filesystem::temp_directory_path() /
        ("kepler_latch_event_cli_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory_);
    std::filesystem::current_path(directory_);
    library_ = "library(test) { cell(LATCH) { latch(IQ, IQN) { enable : E; data_in : D; }"
        "pin(D) { direction : input; } pin(E) { direction : input; }"
        "pin(Q) { direction : output; function : IQ; } } }";
    write("cells.lib", library_);
    write("chain.v", "module top(input d, input e, output q); wire a,b,c;"
        "LATCH l0(.D(d),.E(e),.Q(a)); LATCH l1(.D(a),.E(e),.Q(b));"
        "LATCH l2(.D(b),.E(e),.Q(c)); LATCH l3(.D(c),.E(e),.Q(q)); endmodule\n");
    write("self.v", "module top(input e, output q); LATCH l0(.D(q),.E(e),.Q(q)); endmodule\n");
    write("race.v", "module top(input d, input e, output q, output good);"
        "LATCH l0(.D(d),.E(e),.Q(q)); assign good = d; endmodule\n");
  }
  void TearDown() override {
    std::filesystem::current_path(oldDirectory_);
    std::filesystem::remove_all(directory_);
  }
  void write(const std::string& name, const std::string& contents) {
    std::ofstream(directory_ / name) << contents;
  }
  KEPLER_FORMAL::RunResult run(std::vector<std::string> arguments) {
    arguments.insert(arguments.begin(), "kepler-formal");
    std::vector<char*> argv;
    for (auto& argument : arguments) argv.push_back(argument.data());
    KEPLER_FORMAL::RunResult result;
    const int code = KEPLER_FORMAL::runKeplerFormal(
        static_cast<int>(argv.size()), argv.data(), result);
    EXPECT_EQ(code, result.exitCode);
    return result;
  }
  KEPLER_FORMAL::RunResult config(const std::string& design, const std::string& mode,
                                 const std::string& extra = "",
                                 const std::string& library = "cells.lib") {
    write("run.yaml", "format: verilog\nverification: sec\nsec_encoding: binary\n"
        "input_paths: [" + design + ", " + design + "]\nliberty_files: [" + library + "]\n"
        "latch_support: true\nsec_latch_events: {input_changes: " + mode +
        ", initial_inputs: 0, initial_storage: 0, workers: 2}\n" + extra);
    return run({"--config", "run.yaml"});
  }
  std::vector<std::string> options() {
    return {"--latch_support", "--sec-latch-events", "single", "--sec-latch-initial-inputs", "0",
            "--sec-latch-initial-storage", "0"};
  }
  std::filesystem::path oldDirectory_, directory_;
  std::string library_;
};

TEST_F(LatchEventCliTests, YamlModelsFourTransparentLatchesInAChain) {
  const auto result = config("chain.v", "single");
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.exitCode, 0);
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchEventCliTests, RejectsTuningWithoutMasterGateInBothFrontends) {
  auto arguments = options();
  arguments.erase(arguments.begin());  // Keep tuning, deliberately omit master.
  arguments.insert(arguments.end(), {"-verilog", "-v", "sec", "self.v", "self.v", "cells.lib"});
  EXPECT_NE(run(arguments).exitCode, 0);
  write("disabled.yaml", "format: verilog\nverification: sec\n"
        "input_paths: [self.v, self.v]\nliberty_files: [cells.lib]\n"
        "latch_support: false\nsec_latch_events: " + std::string(complete) + "\n");
  EXPECT_NE(run({"--config", "disabled.yaml"}).exitCode, 0);
}

TEST_F(LatchEventCliTests, MasterDisabledRetainsOpaqueLatchesAndIndependentStrictPolicy) {
  const std::string base = "format: verilog\nverification: sec\n"
      "input_paths: [race.v, race.v]\nliberty_files: [cells.lib]\nlatch_support: false\n";
  write("disabled.yaml", base);
  const auto legacy = run({"--config", "disabled.yaml"});
  EXPECT_EQ(legacy.coveredOutputs, 1u);
  EXPECT_EQ(legacy.totalOutputs, 2u);
  write("disabled.yaml", base + "error_on_opaque: true\n");
  const auto strict = run({"--config", "disabled.yaml"});
  EXPECT_NE(strict.exitCode, 0);
  EXPECT_EQ(strict.status, KEPLER_FORMAL::RunStatus::Unsupported);
  EXPECT_NE(strict.reason.find("error-on-opaque"), std::string::npos);
}

TEST_F(LatchEventCliTests, SelfFeedbackRetainsExplicitInitialStorage) {
  const auto result = config("self.v", "any");
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchEventCliTests, AnyInputChangesKeepRaceOpaqueWhileSingleChangesModelIt) {
  const auto any = config("race.v", "any");
  EXPECT_EQ(any.coveredOutputs, 1u) << any.reason;
  EXPECT_EQ(any.totalOutputs, 2u);
  EXPECT_EQ(any.skippedObservedOutputs.size(), 1u);
  const auto single = config("race.v", "single");
  EXPECT_EQ(single.status, KEPLER_FORMAL::RunStatus::Equivalent) << single.reason;
  EXPECT_EQ(single.coveredOutputs, 2u);
}

TEST_F(LatchEventCliTests, StrictOpacityPolicyStopsAnUncertifiedLatchComponent) {
  const auto result = config("race.v", "any", "error_on_opaque: true\n");
  EXPECT_NE(result.exitCode, 0);
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported);
  EXPECT_NE(result.reason.find("error-on-opaque"), std::string::npos);
}

TEST_F(LatchEventCliTests, FlagsAreAcceptedBeforeAndAfterFormat) {
  const std::vector<std::string> base{
      "-verilog", "-v", "sec", "--sec-encoding", "binary", "self.v", "self.v", "cells.lib"};
  for (const bool before : {true, false}) {
    auto arguments = before ? options() : base;
    const auto suffix = before ? base : options();
    arguments.insert(arguments.end(), suffix.begin(), suffix.end());
    const auto result = run(arguments);
    EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  }
}

TEST_F(LatchEventCliTests, RejectsLecClocklessResetCyclesAndSelectedLeafBoundaries) {
  const std::vector<std::vector<std::string>> incompatible{
      {"-v", "lec"},
      {"-v", "sec", "--sec-reset-cycles", "1", "--sec-reset-port", "e=1"},
      {"-v", "sec", "--set-as-boundary", "l0", "l0"}};
  for (const auto& suffix : incompatible) {
    auto arguments = options();
    arguments.insert(arguments.end(), {"-verilog", "self.v", "self.v", "cells.lib"});
    arguments.insert(arguments.end(), suffix.begin(), suffix.end());
    EXPECT_NE(run(arguments).exitCode, 0);
  }
}

TEST_F(LatchEventCliTests, GzipLibraryIsModeledInCompactMode) {
  auto* file = gzopen((directory_ / "cells.lib.gz").string().c_str(), "wb");
  ASSERT_NE(file, nullptr);
  ASSERT_EQ(gzwrite(file, library_.data(), library_.size()), library_.size());
  ASSERT_EQ(gzclose(file), Z_OK);
  const auto result = config("chain.v", "single", "compact_mode: true\n", "cells.lib.gz");
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchEventCliTests, WorkflowRestoresOuterEventContractAndDisablesItByDefault) {
  namespace Latch = KEPLER_FORMAL::SEC::LATCH;
  Latch::SupportOptions original;
  original.enabled = true;
  original.initialInputs = true;
  original.initialStorage = true;
  original.workers = 3;
  Latch::ScopedSupportOptions outer(original);
  EXPECT_EQ(config("self.v", "single").status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_TRUE(Latch::supportOptions().enabled);
  EXPECT_EQ(Latch::supportOptions().initialInputs, true);
  EXPECT_EQ(Latch::supportOptions().workers, 3u);
  const auto legacy = run({"-verilog", "-v", "sec", "self.v", "self.v", "cells.lib"});
  EXPECT_EQ(legacy.coveredOutputs, 0u);
  EXPECT_TRUE(Latch::supportOptions().enabled);
  EXPECT_EQ(Latch::supportOptions().initialStorage, true);
}

}  // namespace
