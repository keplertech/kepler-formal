// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "Config.h"
#include "KeplerFormalDriver.h"

namespace {

class OpaquePolicyCliTests : public ::testing::Test {
 protected:
  void SetUp() override {
    oldDirectory_ = std::filesystem::current_path();
    oldPolicy_ = KEPLER_FORMAL::Config::getErrorOnOpaque();
    directory_ = std::filesystem::temp_directory_path() /
        ("kepler_opaque_policy_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory_);
    std::filesystem::current_path(directory_);
    write("supported.v", "module top(input a, output y); assign y = a; endmodule\n");
    write("opaque.v", "module top(input a, output y); wire unused; "
          "UNMODELED hidden(.D(a), .Q(unused)); assign y = a; endmodule\n");
    write("cells.lib", "library(test) { cell(UNMODELED) { "
          "pin(D) { direction : input; } pin(Q) { direction : output; } } }\n");
  }
  void TearDown() override {
    KEPLER_FORMAL::Config::setErrorOnOpaque(oldPolicy_);
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
  KEPLER_FORMAL::RunResult config(const std::string& extra,
                                 const std::string& mode = "sec") {
    write("run.yaml", "format: verilog\nverification: " + mode +
          "\ninput_paths: [supported.v, opaque.v]\n"
          "liberty_files: [cells.lib]\n" + extra);
    return run({"--config", "run.yaml"});
  }
  std::filesystem::path oldDirectory_;
  std::filesystem::path directory_;
  bool oldPolicy_ = false;
};

TEST_F(OpaquePolicyCliTests, DefaultStillProvesSupportedOutputWithUnusedOpaqueCell) {
  const auto result = config("");
  EXPECT_EQ(result.exitCode, 0);
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(OpaquePolicyCliTests, YamlEnabledRejectsSecondDesignAndDisabledPreservesBehavior) {
  const auto result = config("error_on_opaque: true\n");
  EXPECT_NE(result.exitCode, 0);
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported);
  EXPECT_NE(result.reason.find("design 2"), std::string::npos);
  EXPECT_NE(result.reason.find("hidden"), std::string::npos);
  EXPECT_EQ(config("error_on_opaque: false\n").status,
            KEPLER_FORMAL::RunStatus::Equivalent);
}

TEST_F(OpaquePolicyCliTests, CompactModeEnforcesPolicyInEitherDesign) {
  for (const bool first : {true, false}) {
    const auto result = run({"-verilog", "-v", "sec", "--compact",
        "--error-on-opaque", first ? "opaque.v" : "supported.v",
        first ? "supported.v" : "opaque.v", "cells.lib"});
    EXPECT_NE(result.exitCode, 0);
    EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported);
    EXPECT_NE(result.reason.find(first ? "design 1" : "design 2"), std::string::npos);
  }
}

TEST_F(OpaquePolicyCliTests, FlagBeforeFormatAndExplicitInputListsAreAccepted) {
  const auto result = run({"--error-on-opaque", "-verilog", "-v", "sec",
      "--design1", "opaque.v", "--design2", "supported.v", "--liberty", "cells.lib"});
  EXPECT_EQ(result.status, KEPLER_FORMAL::RunStatus::Unsupported);
  EXPECT_NE(result.reason.find("error-on-opaque"), std::string::npos);
}

TEST_F(OpaquePolicyCliTests, ScopedRunsResetAndRestorePolicy) {
  KEPLER_FORMAL::Config::setErrorOnOpaque(true);
  EXPECT_EQ(config("").status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getErrorOnOpaque());
  KEPLER_FORMAL::Config::setErrorOnOpaque(false);
  EXPECT_EQ(config("error_on_opaque: true\n").status,
            KEPLER_FORMAL::RunStatus::Unsupported);
  EXPECT_FALSE(KEPLER_FORMAL::Config::getErrorOnOpaque());
}

TEST_F(OpaquePolicyCliTests, RejectsInvalidConfigurationAndNonSecEnablement) {
  for (const auto* value : {"null", "[true]", "{enabled: true}", "perhaps"}) {
    EXPECT_NE(config(std::string("error_on_opaque: ") + value + "\n").exitCode, 0);
  }
  EXPECT_NE(config("error_on_opaque: true\n", "lec").exitCode, 0);
  EXPECT_NE(run({"-verilog", "supported.v", "supported.v",
                 "--error-on-opaque"}).exitCode, 0);
}

}  // namespace
