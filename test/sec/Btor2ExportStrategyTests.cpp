// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <stdexcept>

#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC {
namespace {

SequentialEquivalenceStrategy makeStrategy(Btor2ExportOptions options) {
  return SequentialEquivalenceStrategy(
      nullptr, nullptr, Config::SolverType::KISSAT, SecEngine::Pdr,
      SecEncoding::Binary, {}, options);
}

TEST(Btor2ExportStrategyTests, RejectsDumpOnlyWithoutExportPath) {
  try {
    makeStrategy({"", true});
    FAIL() << "Dump-only without an export path must be rejected";
  } catch (const std::invalid_argument& error) {
    EXPECT_STREQ(error.what(), "BTOR2 dump-only requires an export path");
  }
}

TEST(Btor2ExportStrategyTests, AcceptsDisabledExport) {
  EXPECT_NO_THROW(makeStrategy({}));
}

TEST(Btor2ExportStrategyTests, AcceptsExportBeforeSolving) {
  EXPECT_NO_THROW(makeStrategy({"equivalence.btor2", false}));
}

TEST(Btor2ExportStrategyTests, AcceptsDumpOnlyWithExportPath) {
  EXPECT_NO_THROW(makeStrategy({"equivalence.btor2", true}));
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
