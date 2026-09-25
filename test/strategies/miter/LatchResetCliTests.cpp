// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "KeplerFormalDriver.h"
#include "latch/LatchSupportOptions.h"

namespace {

class LatchResetCliTests : public ::testing::Test {
 protected:
  void SetUp() override {
    previousDirectory_ = std::filesystem::current_path();
    directory_ = std::filesystem::temp_directory_path() /
        ("kepler_latch_reset_cli_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    ASSERT_TRUE(std::filesystem::create_directory(directory_));
    std::filesystem::current_path(directory_);
    write("cells.lib", R"(
library(test) {
  cell(LATCH) {
    latch(IQ, IQN) { enable : E; data_in : D; }
    pin(D) { direction : input; } pin(E) { direction : input; }
    pin(Q) { direction : output; function : IQ; }
  }
  cell(DFF) {
    ff(IQ, IQN) { clocked_on : C; next_state : "D & !R"; }
    pin(D) { direction : input; } pin(C) { direction : input; } pin(R) { direction : input; }
    pin(Q) { direction : output; function : IQ; }
  }
  cell(DFFN) {
    ff(IQ, IQN) { clocked_on : C; next_state : D; clear : "!R"; }
    pin(D) { direction : input; } pin(C) { direction : input; } pin(R) { direction : input; }
    pin(Q) { direction : output; function : IQ; }
  }
  cell(BUF) { pin(A) { direction : input; } pin(Y) { direction : output; function : A; } }
  cell(INV) { pin(A) { direction : input; } pin(Y) { direction : output; function : "!A"; } }
  cell(XOR2) {
    pin(A) { direction : input; } pin(B) { direction : input; }
    pin(Y) { direction : output; function : "A ^ B"; }
  }
})");
    write("chain.v", chain("DFF", "data", "clock"));
    write("different.v", chain("DFF", "inverted_data", "clock",
        "wire inverted_data; INV invert_data(.A(data),.Y(inverted_data));"));
    write("active_low.v", chain("DFFN", "data", "clock"));
    write("buffered.v", chain("DFF", "data", "local_clock",
        "wire local_clock; BUF route(.A(clock),.Y(local_clock));"));
    write("inverted.v", chain("DFF", "data", "local_clock",
        "wire local_clock; INV route(.A(clock),.Y(local_clock));"));
    write("latch_only.v", "module top(input clock, data, enable, reset, output out);"
        " LATCH l(.D(data),.E(enable),.Q(out)); endmodule\n");
    write("two_clocks.v", "module top(input clock, other_clock, data, enable, reset, output out, other_out);"
        " DFF f(.D(data),.C(clock),.R(reset),.Q(out));"
        " DFF g(.D(data),.C(other_clock),.R(reset),.Q(other_out)); endmodule\n");
    write("masked_left.v", "module top(input clock, data, enable, reset, output out);"
        " DFF f(.D(data),.C(clock),.R(reset),.Q(out)); endmodule\n");
    write("masked_right.v", "module top(input clock, data, enable, reset, output out); wire q;"
        " DFF f(.D(data),.C(clock),.R(reset),.Q(q)); XOR2 mask(.A(q),.B(reset),.Y(out)); endmodule\n");
    write("wire.v", "module top(input clock, data, enable, reset, output out);"
        " assign out = data; endmodule\n");
    write("wire_different.v", "module top(input clock, data, enable, reset, output out);"
        " INV invert_data(.A(data),.Y(out)); endmodule\n");
  }
  void TearDown() override {
    std::filesystem::current_path(previousDirectory_);
    std::filesystem::remove_all(directory_);
  }
  std::string chain(const std::string& flop, const std::string& data, const std::string& clock,
                    const std::string& route = "") {
    return "module top(input clock, data, enable, reset, output out); wire q, middle; " + route +
        flop + " f(.D(" + data + "),.C(" + clock + "),.R(reset),.Q(q));"
        " LATCH a(.D(q),.E(enable),.Q(middle)); LATCH b(.D(middle),.E(enable),.Q(out)); endmodule\n";
  }
  void write(const std::string& file, const std::string& text) { std::ofstream(directory_ / file) << text; }
  KEPLER_FORMAL::RunResult run(std::vector<std::string> arguments) {
    arguments.insert(arguments.begin(), "kepler-formal");
    std::vector<char*> argv;
    for (auto& argument : arguments) argv.push_back(argument.data());
    KEPLER_FORMAL::RunResult result;
    const int code = KEPLER_FORMAL::runKeplerFormal(static_cast<int>(argv.size()), argv.data(), result);
    EXPECT_EQ(code, result.exitCode);
    return result;
  }
  KEPLER_FORMAL::RunResult yaml(const std::string& first, const std::string& second,
      bool compact = false, bool activeHigh = true, const std::string& extraPorts = "") {
    write("run.yaml", "format: verilog\nverification: sec\nsec_engine: pdr\nsec_encoding: binary\n"
        "max_k: 48\ninput_paths: [" + first + ", " + second + "]\nliberty_files: [cells.lib]\n"
        "compact_mode: " + (compact ? "true" : "false") + "\n"
        "latch_support: true\n"
        "sec_latch_events: {input_changes: single, initial_inputs: 0, initial_storage: 0, workers: 2}\n"
        "sec_reset:\n  cycles: 2\n  ports:\n    - name: reset\n      active_value: " +
        (activeHigh ? "1" : "0") + "\n" + extraPorts);
    return run({"--config", "run.yaml"});
  }
  std::filesystem::path previousDirectory_, directory_;
};

TEST_F(LatchResetCliTests, YamlResetCyclesWorkWithIndependentLatchEnableAndCompactExtraction) {
  for (const bool compact : {false, true}) {
    const auto result = yaml("chain.v", "chain.v", compact);
    EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
    EXPECT_EQ(result.exitCode, 0);
    EXPECT_EQ(result.coveredOutputs, 1u);
  }
}

TEST_F(LatchResetCliTests, ResetFlagsComposeWithLatchSupportFlags) {
  const auto result = run({"--latch_support", "--sec-latch-events", "single",
      "--sec-latch-initial-inputs", "0", "--sec-latch-initial-storage", "0",
      "--sec-reset-cycles", "2", "--sec-reset-port", "reset=1",
      "-verilog", "-v", "sec", "--sec-encoding", "binary", "chain.v", "chain.v", "cells.lib"});
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchResetCliTests, BufferedAndInvertedClockRoutesAreDiscovered) {
  for (const auto* design : {"buffered.v", "inverted.v"}) {
    const auto result = yaml(design, design);
    EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
    EXPECT_EQ(result.coveredOutputs, 1u);
  }
}

TEST_F(LatchResetCliTests, DifferentAfterResetRemainsDifferentInBothExtractionModes) {
  for (const bool compact : {false, true}) {
    const auto result = yaml("chain.v", "different.v", compact);
    EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Different) << result.reason;
    EXPECT_NE(result.exitCode, 0);
    EXPECT_EQ(result.coveredOutputs, 1u);
  }
}

TEST_F(LatchResetCliTests, ResetOnlyMismatchIsNotComparedAndResetCannotBeReasserted) {
  const auto result = yaml("masked_left.v", "masked_right.v", true);
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchResetCliTests, ActiveLowAsynchronousResetIsAccepted) {
  const auto result = yaml("active_low.v", "active_low.v", true, false);
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchResetCliTests, MissingEdgeClockDoesNotTreatLatchEnableAsAClock) {
  const auto result = yaml("latch_only.v", "latch_only.v");
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("flip-flop clock"), std::string::npos) << result.reason;
  EXPECT_EQ(result.reason.find("cannot use clock-cycle"), std::string::npos) << result.reason;
}

TEST_F(LatchResetCliTests, MultipleClockRootsHaveASpecificDiagnostic) {
  const auto result = yaml("two_clocks.v", "two_clocks.v");
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("multiple independent"), std::string::npos) << result.reason;
}

TEST_F(LatchResetCliTests, MultipleResetPortsAreRejectedExplicitly) {
  const auto result = yaml("chain.v", "chain.v", false, true,
      "    - name: enable\n      active_value: 0\n");
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("reset"), std::string::npos) << result.reason;
  EXPECT_EQ(result.reason.find("cannot use clock-cycle"), std::string::npos) << result.reason;
}

TEST_F(LatchResetCliTests, OffOnOffPreservesLegacyResetSemanticsForEveryEngine) {
  const auto legacy = [&](const std::string& engine, const std::string& second) {
    // Legacy reset supports multiple reset ports and needs no discovered FF
    // carrier. Either leaked event dispatch or leaked reset sampling rejects
    // this deliberately clockless design instead of preserving that behavior.
    write("legacy.yaml", "format: verilog\nverification: sec\nsec_engine: " + engine +
        "\nsec_encoding: binary\nmax_k: 8\nlatch_support: false\n"
        "input_paths: [wire.v, " + second + "]\nliberty_files: [cells.lib]\n"
        "sec_reset:\n  cycles: 2\n  ports:\n"
        "    - name: reset\n      active_value: 1\n"
        "    - name: enable\n      active_value: 0\n");
    return run({"--config", "legacy.yaml"});
  };
  for (const auto* engine : {"k_induction", "imc", "pdr"}) {
    SCOPED_TRACE(engine);
    for (const auto* second : {"wire.v", "wire_different.v"}) {
      SCOPED_TRACE(second);
      const auto before = legacy(engine, second);
      const auto enabled = yaml("chain.v", "chain.v");
      EXPECT_EQ(enabled.status, KEPLER_FORMAL::RunStatus::Equivalent) << enabled.reason;
      const auto after = legacy(engine, second);
      EXPECT_EQ(before.status, std::string(second) == "wire.v"
          ? KEPLER_FORMAL::RunStatus::Equivalent : KEPLER_FORMAL::RunStatus::Different) << before.reason;
      EXPECT_EQ(before.coveredOutputs, 1u) << before.reason;
      EXPECT_EQ(after.status, before.status) << after.reason;
      EXPECT_EQ(after.exitCode, before.exitCode);
      EXPECT_EQ(after.coveredOutputs, before.coveredOutputs);
      EXPECT_EQ(after.totalOutputs, before.totalOutputs);
      EXPECT_EQ(after.skippedObservedOutputs, before.skippedObservedOutputs);
      EXPECT_EQ(after.reason, before.reason);
      EXPECT_EQ(after.reason.find("Event contract"), std::string::npos);
      EXPECT_EQ(after.reason.find("boolean-epochs"), std::string::npos);
      EXPECT_EQ(after.reason.find("reset/event step"), std::string::npos);
      EXPECT_FALSE(KEPLER_FORMAL::SEC::LATCH::supportOptions().enabled);
      EXPECT_FALSE(KEPLER_FORMAL::SEC::LATCH::supportOptions().initialInputs.has_value());
    }
  }
}

}  // namespace
