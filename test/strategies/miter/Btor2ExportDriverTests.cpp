// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "KeplerFormalDriver.h"
#include "NLUniverse.h"

namespace {

class Btor2ExportDriverTests : public ::testing::Test {
 protected:
  void SetUp() override {
    directory_ = std::filesystem::temp_directory_path() /
                 ("kepler_btor2_driver_" + std::to_string(
                     std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory_);
    std::ofstream(directory_ / "reference.v")
        << "module top(input a, input b, output y); or (y, a, b); endmodule\n";
    std::ofstream(directory_ / "different.v")
        << "module top(input a, input b, output y); and (y, a, b); endmodule\n";
  }

  void TearDown() override {
    // Also clean up on an assertion or an unexpected exception before the
    // in-process driver has returned normally.
    if (naja::NL::NLUniverse::get() != nullptr) {
      naja::NL::NLUniverse::get()->destroy();
    }
    KEPLER_FORMAL::cleanupKeplerFormalState();
    std::filesystem::remove_all(directory_);
  }

  int run(bool dumpOnly, KEPLER_FORMAL::RunResult& result) {
    std::vector<std::string> args = {
        "kepler-formal", "-verilog", "-v", "sec", "--sec-encoding", "binary",
        (directory_ / "reference.v").string(),
        (directory_ / "different.v").string(),
        "--dump-btor2", (directory_ / "problem.btor2").string()};
    if (dumpOnly) {
      args.push_back("--dump-only");
    }
    std::vector<char*> argv;
    for (auto& arg : args) {
      argv.push_back(arg.data());
    }
    return KEPLER_FORMAL::runKeplerFormal(
        static_cast<int>(argv.size()), argv.data(), result);
  }

  std::filesystem::path directory_;
};

TEST_F(Btor2ExportDriverTests, DumpOnlyReportsExportedWithoutClaimingProof) {
  KEPLER_FORMAL::RunResult result;
  ASSERT_EQ(0, run(true, result));
  EXPECT_EQ(0, result.exitCode);
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Exported, result.status);
  EXPECT_STREQ("exported", KEPLER_FORMAL::runStatusName(result.status));
  EXPECT_EQ("sec", result.verification);
  EXPECT_EQ(1, result.coveredOutputs);
  EXPECT_EQ(1, result.totalOutputs);
  EXPECT_EQ(0, result.provenOutputs);
  EXPECT_NE(std::string::npos, result.reason.find("proof not run"));
  EXPECT_EQ(nullptr, naja::NL::NLUniverse::get());
  ASSERT_TRUE(std::filesystem::is_regular_file(directory_ / "problem.btor2"));
  std::ifstream file(directory_ / "problem.btor2");
  const std::string contents{std::istreambuf_iterator<char>(file),
                             std::istreambuf_iterator<char>()};
  EXPECT_NE(std::string::npos, contents.find(" bad "));
}

TEST_F(Btor2ExportDriverTests, ExportThenSolveReplacesPreviousExportedResult) {
  KEPLER_FORMAL::RunResult result;
  ASSERT_EQ(0, run(true, result));
  ASSERT_EQ(KEPLER_FORMAL::RunStatus::Exported, result.status);
  ASSERT_EQ(3, run(false, result));
  EXPECT_EQ(3, result.exitCode);
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Different, result.status);
  EXPECT_EQ(0, result.provenOutputs);
  EXPECT_EQ(nullptr, naja::NL::NLUniverse::get());
  EXPECT_TRUE(std::filesystem::is_regular_file(directory_ / "problem.btor2"));
}

}  // namespace
