// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "Btor2ExportConfig.h"

#include <yaml-cpp/yaml.h>

namespace KEPLER_FORMAL {
namespace {

bool readBoolean(
    const YAML::Node& config, const char* key, bool& value, std::string& error) {
  const auto node = config[key];
  if (!node) {
    return true;
  }
  if (node.IsScalar()) {
    try {
      value = node.as<bool>();
      return true;
    } catch (const YAML::Exception&) {
    }
  }
  error = std::string(key) + " must be a boolean";
  return false;
}

}  // namespace

bool Btor2ExportConfig::parseYaml(
    const YAML::Node& config, std::string& error) {
  const auto enabledNode = config["btor2_export"];
  const auto pathNode = config["btor2_export_path"];
  const auto dumpOnlyNode = config["dump_only"];
  explicit_ = enabledNode.IsDefined() || pathNode.IsDefined() ||
              dumpOnlyNode.IsDefined();
  bool enabled = false;
  if (!readBoolean(config, "btor2_export", enabled, error) ||
      !readBoolean(config, "dump_only", options_.dumpOnly, error)) {
    return false;
  }
  std::string path = "miter.btor2";
  if (pathNode) {
    if (!pathNode.IsScalar() || pathNode.as<std::string>().empty()) {
      error = "btor2_export_path must be a non-empty string";
      return false;
    }
    path = pathNode.as<std::string>();
    if (!enabled) {
      error = "btor2_export_path requires btor2_export: true";
      return false;
    }
  }
  options_.path = enabled ? path : std::string{};
  return true;
}

Btor2ExportConfig::ArgumentResult Btor2ExportConfig::parseArgument(
    int argc, char** argv, int& index, std::string& error) {
  const std::string argument = argv[index];
  if (argument == "--dump-only") {
    explicit_ = true;
    options_.dumpOnly = true;
    return ArgumentResult::Parsed;
  }
  if (argument != "--dump-btor2") {
    return ArgumentResult::NotHandled;
  }
  explicit_ = true;
  if (index + 1 >= argc || std::string(argv[index + 1]).empty() ||
      std::string(argv[index + 1]).rfind("--", 0) == 0) {
    error = "--dump-btor2 requires a non-empty output file path";
    return ArgumentResult::Error;
  }
  options_.path = argv[++index];
  return ArgumentResult::Parsed;
}

bool Btor2ExportConfig::validate(bool isSec, std::string& error) const {
  if (explicit_ && !isSec) {
    error = "BTOR2 export options are only supported with SEC verification";
    return false;
  }
  if (options_.dumpOnly && !options_.enabled()) {
    error = "dump_only/--dump-only requires btor2_export: true or --dump-btor2 <file>";
    return false;
  }
  return true;
}

}  // namespace KEPLER_FORMAL
