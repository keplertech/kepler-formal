// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>

#include "BoolExprCache.h"
#include "Tree2BoolExpr.h"

extern int KeplerFormalMain(int argc, char** argv);

namespace {

class Btor2ExportCliTests : public ::testing::Test {
 protected:
  void SetUp() override {
    originalDirectory_ = std::filesystem::current_path();
    directory_ = std::filesystem::temp_directory_path() /
                 ("kepler_btor2_cli_" + std::to_string(
                     std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory_);
    std::filesystem::current_path(directory_);
    // This pair differs at frame zero, so export-only success cannot be
    // confused with the normal proof result (counterexample, exit code 3).
    write("design0.v", "module top(input a, input b, output y);\n"
                       "or (y, a, b);\nendmodule\n");
    write("design1.v", "module top(input a, input b, output y);\n"
                       "and (y, a, b);\nendmodule\n");
  }

  void TearDown() override {
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
    std::filesystem::current_path(originalDirectory_);
    std::filesystem::remove_all(directory_);
  }

  void write(const std::string& name, const std::string& contents) {
    std::ofstream(directory_ / name) << contents;
  }

  std::string read(const std::string& name) const {
    std::ifstream file(directory_ / name);
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
  }

  int run(std::vector<std::string> args) {
    args.insert(args.begin(), "kepler-formal");
    std::vector<char*> argv;
    for (auto& argument : args) {
      argv.push_back(argument.data());
    }
    const int result = KeplerFormalMain(static_cast<int>(argv.size()), argv.data());
    if (auto logger = spdlog::get("kepler_formal_main_logger")) {
      logger->flush();
    }
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
    return result;
  }

  int runConfig(const std::string& options, const std::string& mode = "sec") {
    write("config.yaml",
          "format: verilog\nverification: " + mode + "\n"
          "input_paths: [design0.v, design1.v]\n"
          "log_file: run.log\n" + options);
    return run({"--config", "config.yaml"});
  }

  void expectExport(const std::string& path) const {
    ASSERT_TRUE(std::filesystem::exists(directory_ / path));
    const auto contents = read(path);
    EXPECT_NE(contents.find("sort bitvec 1"), std::string::npos);
    EXPECT_NE(contents.find(" bad "), std::string::npos);
  }

  std::filesystem::path directory_;
  std::filesystem::path originalDirectory_;
};

TEST_F(Btor2ExportCliTests, YamlExportsToDefaultPathWithoutRunningProof) {
  ASSERT_EQ(runConfig("btor2_export: true\ndump_only: true\n"), EXIT_SUCCESS);
  expectExport("miter.btor2");
  const auto log = read("run.log");
  EXPECT_NE(log.find("SEC BTOR2 exported to miter.btor2; proof not run"),
            std::string::npos);
  EXPECT_NE(log.find("Export covers 1/1 observed outputs"), std::string::npos);
  EXPECT_EQ(log.find("SEC proved equivalence"), std::string::npos);
  EXPECT_EQ(log.find("SEC found a counterexample"), std::string::npos);
}

TEST_F(Btor2ExportCliTests, YamlCustomPathExportsThenSolves) {
  EXPECT_EQ(runConfig("btor2_export: true\n"
                      "btor2_export_path: compared.btor2\n"
                      "dump_only: false\n"), 3);
  expectExport("compared.btor2");
  EXPECT_NE(read("run.log").find("SEC found a counterexample"), std::string::npos);
}

TEST_F(Btor2ExportCliTests, YamlCompactModeExports) {
  ASSERT_EQ(runConfig("btor2_export: true\n"
                      "btor2_export_path: compact.btor2\n"
                      "dump_only: true\ncompact_mode: true\n"), EXIT_SUCCESS);
  expectExport("compact.btor2");
}

TEST_F(Btor2ExportCliTests, CompactIdenticalInputReuseExports) {
  ASSERT_EQ(run({"-verilog", "-v", "sec", "--compact",
                 "design0.v", "design0.v", "--dump-btor2", "reused.btor2",
                 "--dump-only"}), EXIT_SUCCESS);
  expectExport("reused.btor2");
}

TEST_F(Btor2ExportCliTests, YamlDisabledDoesNotExport) {
  EXPECT_EQ(runConfig("btor2_export: false\ndump_only: false\n"), 3);
  EXPECT_FALSE(std::filesystem::exists(directory_ / "miter.btor2"));
}

TEST_F(Btor2ExportCliTests, YamlRejectsInvalidTypesAndCombinations) {
  const std::vector<std::string> invalidOptions = {
      "btor2_export: [true]\n",
      "btor2_export: perhaps\n",
      "btor2_export: null\n",
      "btor2_export: true\nbtor2_export_path: [output.btor2]\n",
      "btor2_export: true\nbtor2_export_path: ''\n",
      "btor2_export: true\nbtor2_export_path: null\n",
      "btor2_export: true\ndump_only: {enabled: true}\n",
      "btor2_export: true\ndump_only: perhaps\n",
      "dump_only: true\n",
      "btor2_export: false\ndump_only: true\n",
      "btor2_export_path: output.btor2\n",
  };
  for (const auto& options : invalidOptions) {
    SCOPED_TRACE(options);
    EXPECT_EQ(runConfig(options), EXIT_FAILURE);
  }
  EXPECT_FALSE(std::filesystem::exists(directory_ / "miter.btor2"));
}

TEST_F(Btor2ExportCliTests, RejectsExportOptionsInLec) {
  EXPECT_EQ(runConfig("btor2_export: true\n", "lec"), EXIT_FAILURE);
  EXPECT_EQ(run({"-verilog", "design0.v", "design1.v",
                 "--dump-btor2", "output.btor2"}), EXIT_FAILURE);
  EXPECT_FALSE(std::filesystem::exists(directory_ / "output.btor2"));
}

TEST_F(Btor2ExportCliTests, CliExportOptionsBeforeInputFormat) {
  ASSERT_EQ(run({"--dump-only", "--dump-btor2", "before.btor2",
                 "-v", "sec", "-verilog", "design0.v", "design1.v"}),
            EXIT_SUCCESS);
  expectExport("before.btor2");
}

TEST_F(Btor2ExportCliTests, CliExportOptionsAfterInputFiles) {
  ASSERT_EQ(run({"-verilog", "design0.v", "design1.v", "-v", "sec",
                 "--dump-btor2", "after.btor2", "--dump-only"}), EXIT_SUCCESS);
  expectExport("after.btor2");
}

TEST_F(Btor2ExportCliTests, CliExplicitDesignListsAndExportOptions) {
  ASSERT_EQ(run({"-verilog", "-v", "sec", "--design1", "design0.v",
                 "--dump-btor2", "lists.btor2", "--design2", "design1.v",
                 "--dump-only"}), EXIT_SUCCESS);
  expectExport("lists.btor2");
}

TEST_F(Btor2ExportCliTests, CliRejectsMissingPathAndDumpOnlyWithoutExport) {
  EXPECT_EQ(run({"--dump-btor2"}), EXIT_FAILURE);
  EXPECT_EQ(run({"--dump-btor2", "--dump-only", "-verilog", "-v", "sec",
                 "design0.v", "design1.v"}), EXIT_FAILURE);
  EXPECT_EQ(run({"-verilog", "-v", "sec", "design0.v", "design1.v",
                 "--dump-btor2"}), EXIT_FAILURE);
  EXPECT_EQ(run({"-verilog", "-v", "sec", "design0.v", "design1.v",
                 "--dump-btor2", ""}), EXIT_FAILURE);
  EXPECT_EQ(run({"-verilog", "-v", "sec", "design0.v", "design1.v",
                 "--dump-only"}), EXIT_FAILURE);
}

TEST_F(Btor2ExportCliTests, ExportWriteFailureDoesNotReportExported) {
  std::filesystem::create_directory(directory_ / "blocked.btor2");
  EXPECT_NE(runConfig("btor2_export: true\n"
                      "btor2_export_path: blocked.btor2\n"
                      "dump_only: true\n"), EXIT_SUCCESS);
  EXPECT_EQ(read("run.log").find("SEC BTOR2 exported to"), std::string::npos);
}

}  // namespace
