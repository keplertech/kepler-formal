// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "KeplerFormalDriver.h"
#include "NLUniverse.h"

namespace {

class DriverSeparationTests : public ::testing::Test {
 protected:
  void SetUp() override {
    directory_ = std::filesystem::temp_directory_path() /
        ("kepler_driver_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory_);
    std::ofstream(directory_ / "design.v")
        << "module top(input a, output y); assign y = a; endmodule\n";
    std::ofstream(directory_ / "run.yaml")
        << "format: verilog\ninput_paths:\n  - "
        << (directory_ / "design.v").string() << "\n  - "
        << (directory_ / "design.v").string()
        << "\npy_tech_files: [host_primitives.py]\n";
  }

  void TearDown() override {
    KEPLER_FORMAL::cleanupKeplerFormalState();
    if (auto* universe = naja::NL::NLUniverse::get()) universe->destroy();
    std::filesystem::remove_all(directory_);
  }

  template<class Run>
  KEPLER_FORMAL::RunResult run(Run invoke) {
    std::vector<std::string> arguments = {
        "kepler-formal", "--config", (directory_ / "run.yaml").string()};
    std::vector<char*> argv;
    for (auto& argument : arguments) argv.push_back(argument.data());
    KEPLER_FORMAL::RunResult result;
    const int code = invoke(static_cast<int>(argv.size()), argv.data(), result);
    EXPECT_EQ(code, result.exitCode);
    EXPECT_EQ(nullptr, naja::NL::NLUniverse::get());
    return result;
  }

  std::filesystem::path directory_;
};

// A host can supply primitives without embedding Python in the common core.
// These designs use built-in assign primitives, so this host has none to add.
class HostPrimitiveLoader final : public KEPLER_FORMAL::PrimitiveLibraryLoader {
 public:
  void prepare(const char*) const override {
    EXPECT_EQ(nullptr, naja::NL::NLUniverse::get());
    if (reject) throw std::runtime_error("Host rejected primitive loading");
  }
  void load(naja::NL::NLLibrary*, const std::filesystem::path&) const override {
    loaded = true;
  }
  bool reject = false;
  mutable bool loaded = false;
};

TEST_F(DriverSeparationTests, SharedRunUsesEachDriversPrimitivePolicy) {
  HostPrimitiveLoader host;
  auto runHost = [&](int argc, char** argv, KEPLER_FORMAL::RunResult& result) {
    return KEPLER_FORMAL::runKeplerFormal(argc, argv, result, host);
  };
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Equivalent, run(runHost).status);
  EXPECT_TRUE(host.loaded);

  const auto defaultHost = run([](int argc, char** argv, KEPLER_FORMAL::RunResult& result) {
    return KEPLER_FORMAL::runKeplerFormal(argc, argv, result);
  });
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Error, defaultHost.status);
  EXPECT_NE(std::string::npos, defaultHost.reason.find("in-process file API"));

  // One host's rejection must not affect a following call from another host.
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Equivalent, run(runHost).status);
}

TEST_F(DriverSeparationTests, PreparationFailureLeavesTheSharedRunReusable) {
  HostPrimitiveLoader host;
  host.reject = true;
  auto runHost = [&](int argc, char** argv, KEPLER_FORMAL::RunResult& result) {
    return KEPLER_FORMAL::runKeplerFormal(argc, argv, result, host);
  };
  const auto failure = run(runHost);
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Error, failure.status);
  EXPECT_EQ("Host rejected primitive loading", failure.reason);
  EXPECT_FALSE(host.loaded);

  host.reject = false;
  EXPECT_EQ(KEPLER_FORMAL::RunStatus::Equivalent, run(runHost).status);
}

TEST_F(DriverSeparationTests, VersionRequestReturnsNoResultWithoutLoading) {
  HostPrimitiveLoader host;
  host.reject = true;
  for (const auto* flag : {"--version", "-V"}) {
    std::vector<std::string> arguments = {"kepler-formal", flag};
    std::vector<char*> argv;
    for (auto& argument : arguments) argv.push_back(argument.data());
    KEPLER_FORMAL::RunResult result;
    result.status = KEPLER_FORMAL::RunStatus::Equivalent;
    result.reason = "previous verification";
    result.bound = 10;

    testing::internal::CaptureStdout();
    const int code = KEPLER_FORMAL::runKeplerFormal(
        static_cast<int>(argv.size()), argv.data(), result, host);
    const auto output = testing::internal::GetCapturedStdout();

    EXPECT_EQ(EXIT_SUCCESS, code);
    EXPECT_EQ(code, result.exitCode);
    EXPECT_EQ(KEPLER_FORMAL::RunStatus::NoResult, result.status);
    EXPECT_TRUE(result.reason.empty());
    EXPECT_EQ(0, result.bound);
    EXPECT_EQ(0, output.find("kepler-formal version: "));
    EXPECT_FALSE(host.loaded);
    EXPECT_EQ(nullptr, naja::NL::NLUniverse::get());
  }
}

}  // namespace
