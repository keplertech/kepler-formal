// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <string>
#include "latch/LatchSupportOptions.h"
namespace YAML { class Node; }
namespace KEPLER_FORMAL {
class LatchEventConfig {
 public:
  enum class ArgumentResult { NotHandled, Parsed, Error };
  bool parseYaml(const YAML::Node& config, std::string& error);
  ArgumentResult parseArgument(int argc, char** argv, int& index, std::string& error);
  bool validate(bool isSec, bool hasResetCycles, bool hasLeafBoundaries, std::string& error) const;
  const SEC::LATCH::SupportOptions& options() const { return options_; }
 private:
  SEC::LATCH::SupportOptions options_;
  bool explicitTuning_ = false;
  bool inputChangesExplicit_ = false;
  bool set(const std::string& key, const std::string& value, std::string& error);
};
}  // namespace KEPLER_FORMAL
