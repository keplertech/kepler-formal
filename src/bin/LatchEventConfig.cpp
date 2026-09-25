// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "LatchEventConfig.h"
#include <charconv>
#include <limits>
#include <map>
#include <string_view>
#include <yaml-cpp/yaml.h>

namespace KEPLER_FORMAL {
namespace {
bool number(const std::string& text, size_t& value) {
  if (text.empty()) return false;
  const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size();
}
}
bool LatchEventConfig::set(const std::string& key, const std::string& value, std::string& error) {
  options_.explicitConfiguration = true;
  if (key == "input_changes") {
    if (value == "any" || value == "single") {
      options_.singleInputChange = value == "single";
      return true;
    }
    error = "sec_latch_events.input_changes must be any or single";
    return false;
  }
  if (key == "initial_inputs" || key == "initial_storage") {
    if (value != "0" && value != "1") {
      error = "sec_latch_events." + key + " must explicitly be Boolean 0 or 1";
      return false;
    }
    (key == "initial_inputs" ? options_.initialInputs : options_.initialStorage) = value == "1";
    return true;
  }
  size_t count = 0;
  if (!number(value, count)) {
    error = "sec_latch_events." + key + " must be a nonnegative integer";
    return false;
  }
  if (key == "workers" && count <= size_t(std::numeric_limits<int>::max())) options_.workers = count;
  else if (key == "max_waves" && count) options_.limits.maxWaves = count;
  else if (key == "max_states" && count) options_.limits.maxBoundaryStates = count;
  else if (key == "max_transactions" && count) options_.limits.maxTransactions = count;
  else if (key == "max_symbolic_nodes" && count) options_.maxSymbolicNodes = count;
  else if (key == "max_sat_conflicts" && count && count <= std::numeric_limits<unsigned>::max()) options_.maxSatConflicts = count;
  else if (key == "max_sat_decisions" && count && count <= std::numeric_limits<unsigned>::max()) options_.maxSatDecisions = count;
  else {
    error = "unknown or invalid sec_latch_events option: " + key;
    return false;
  }
  return true;
}

bool LatchEventConfig::parseYaml(const YAML::Node& config, std::string& error) {
  if (const auto gate = config["latch_support"]) {
    if (!gate.IsScalar() ||
        (gate.as<std::string>() != "true" && gate.as<std::string>() != "false")) {
      error = "latch_support must be true or false";
      return false;
    }
    options_.enabled = gate.as<std::string>() == "true";
  }
  const auto node = config["sec_latch_events"];
  if (!node) return true;
  options_.explicitConfiguration = true;
  if (!node.IsMap()) { error = "sec_latch_events must be a map"; return false; }
  for (auto item : node) {
    if (!item.first.IsScalar() || !item.second.IsScalar()) {
      error = "sec_latch_events entries must be scalar key/value pairs";
      return false;
    }
    if (!set(item.first.as<std::string>(), item.second.as<std::string>(), error)) return false;
  }
  return true;
}

LatchEventConfig::ArgumentResult LatchEventConfig::parseArgument(
    int argc, char** argv, int& index, std::string& error) {
  const std::string_view argument(argv[index]);
  if (argument == "--latch_support" || argument == "--no-latch_support") {
    options_.enabled = argument == "--latch_support";
    return ArgumentResult::Parsed;
  }
  static const std::map<std::string, std::string> names{
      {"--sec-latch-events", "input_changes"},
      {"--sec-latch-initial-inputs", "initial_inputs"},
      {"--sec-latch-initial-storage", "initial_storage"},
      {"--sec-latch-workers", "workers"},
      {"--sec-latch-max-waves", "max_waves"},
      {"--sec-latch-max-states", "max_states"},
      {"--sec-latch-max-transactions", "max_transactions"},
      {"--sec-latch-max-nodes", "max_symbolic_nodes"},
      {"--sec-latch-sat-conflicts", "max_sat_conflicts"},
      {"--sec-latch-sat-decisions", "max_sat_decisions"}};
  const auto option = names.find(argv[index]);
  if (option == names.end()) return ArgumentResult::NotHandled;
  if (index + 1 == argc) { error = option->first + " requires a value"; return ArgumentResult::Error; }
  return set(option->second, argv[++index], error) ? ArgumentResult::Parsed : ArgumentResult::Error;
}

bool LatchEventConfig::validate(bool isSec, bool hasResetCycles, bool hasLeafBoundaries,
                                std::string& error) const {
  (void)hasResetCycles;  // Clock discovery and cycle expansion need the extracted models.
  if (!options_.enabled) {
    if (!options_.explicitConfiguration) return true;
    error = "latch event tuning requires latch_support: true or --latch_support";
    return false;
  }
  // No tuning is required: unrestricted input changes and symbolic initial
  // values are the default. LEC and selected-leaf workflows stay unchanged
  // unless the user explicitly supplies SEC event options.
  if (!options_.explicitConfiguration) return true;
  if (!isSec) error = "latch event options require SEC verification";
  else if (hasLeafBoundaries)
    error = "latch events currently require the complete top interface, not selected leaf boundaries";
  else return true;
  return false;
}
}  // namespace KEPLER_FORMAL
