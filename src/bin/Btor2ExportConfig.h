// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <string>

#include "export/Btor2ExportOptions.h"

namespace YAML {
class Node;
}

namespace KEPLER_FORMAL {

// Keeps export-specific CLI/YAML parsing outside the main input loader.
class Btor2ExportConfig {
 public:
  enum class ArgumentResult { NotHandled, Parsed, Error };

  bool parseYaml(const YAML::Node& config, std::string& error);

  // On success, index points to the last consumed argument. The caller
  // advances it once before parsing the next option.
  ArgumentResult parseArgument(
      int argc, char** argv, int& index, std::string& error);

  bool validate(bool isSec, std::string& error) const;
  const SEC::Btor2ExportOptions& options() const { return options_; }

 private:
  SEC::Btor2ExportOptions options_;
  bool explicit_ = false;
};

}  // namespace KEPLER_FORMAL
