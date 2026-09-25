// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "LibertyLatchModels.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <spdlog/spdlog.h>
#include <zlib.h>

#include "NLLibrary.h"
#include "SNLBooleanTree.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLLibertyConstructor.h"
#include "SNLScalarTerm.h"
#include "YosysLibertyParser.h"

namespace KEPLER_FORMAL {
namespace {

using namespace naja::NL;
using Ast = Yosys::LibertyAst;
using Modeling = SNLDesignModeling;
using Expression = Modeling::BooleanExpression;
using Conflict = Modeling::SequentialState::ClearPresetValue;

const Ast* child(const Ast* parent, const std::string& name) {
  const Ast* found = nullptr;
  for (const auto* candidate : parent->children) {
    if (candidate->id != name) continue;
    if (found) throw std::runtime_error("duplicate `" + name + "` attribute");
    found = candidate;
  }
  return found;
}

bool contains(const Ast* parent, const std::string& name) {
  if (parent->id == name) return true;
  return std::any_of(parent->children.begin(), parent->children.end(),
      [&](const Ast* candidate) { return contains(candidate, name); });
}

const Ast* requireChild(const Ast* parent, const std::string& name) {
  const auto* result = child(parent, name);
  if (!result || result->value.empty()) {
    throw std::runtime_error("missing `" + name + "` expression");
  }
  return result;
}

Conflict conflictValue(const Ast* latch, const std::string& name) {
  const auto* value = child(latch, name);
  if (!value) return Conflict::Unknown;
  if (value->value == "L") return Conflict::Zero;
  if (value->value == "H") return Conflict::One;
  if (value->value == "N") return Conflict::Hold;
  if (value->value == "T") return Conflict::Toggle;
  if (value->value == "X") return Conflict::Unknown;
  throw std::runtime_error("unsupported `" + name + "` value");
}

Conflict complement(Conflict value) {
  if (value == Conflict::Zero) return Conflict::One;
  if (value == Conflict::One) return Conflict::Zero;
  return value;
}

Expression expression(SNLDesign* primitive, const std::string& text,
                       const SNLBooleanTree::StateIdentifiers& states) {
  SNLBooleanTree tree;
  tree.parse(primitive, text, states);
  auto result = tree.getBooleanExpression();
  if (!result.isValid()) throw std::runtime_error("invalid Boolean expression");
  for (const auto& node : result.nodes) {
    if (node.operation == Expression::Operator::Term &&
        (!node.term || node.term->getDirection() != SNLTerm::Direction::Input)) {
      throw std::runtime_error("expression refers to a non-input pin");
    }
  }
  return result;
}

std::set<SNLBitTerm*> terms(const Expression& expression) {
  std::set<SNLBitTerm*> result;
  for (const auto& node : expression.nodes) {
    if (node.operation == Expression::Operator::Term) result.insert(node.term);
  }
  return result;
}

void addArcs(const Modeling::SequentialModel& model) {
  std::set<SNLBitTerm*> updateInputs;
  for (const auto& state : model.states) {
    for (const auto* expression : {&state.nextState,
         state.clear ? &*state.clear : nullptr,
         state.preset ? &*state.preset : nullptr}) {
      if (expression) {
        const auto inputs = terms(*expression);
        updateInputs.insert(inputs.begin(), inputs.end());
      }
    }
  }
  for (auto* enable : terms(model.clockedOn)) {
    Modeling::setTermRole(enable, Modeling::SNLTermRole::Clock);
    for (auto* input : updateInputs) {
      const auto existing = Modeling::getInputRelatedClocks(input);
      if (std::find(existing.begin(), existing.end(), enable) == existing.end()) {
        Modeling::addInputsToClockArcs({input}, enable);
      }
    }
    for (const auto& output : model.outputs) {
      Modeling::setTermRole(output.term, Modeling::SNLTermRole::DataOutput);
      const auto existing = Modeling::getOutputRelatedClocks(output.term);
      if (std::find(existing.begin(), existing.end(), enable) == existing.end()) {
        Modeling::addClockToOutputsArcs(enable, {output.term});
      }
    }
  }
}

void populate(SNLDesign* primitive, const Ast* cell) {
  for (const auto* forbidden : {"ff", "ff_bank", "latch_bank", "memory",
                               "statetable", "state_function", "bus", "bundle",
                               "power_down_function"}) {
    if (contains(cell, forbidden)) {
      throw std::runtime_error(std::string("unsupported latch cell containing `") +
                               forbidden + "`");
    }
  }
  std::vector<const Ast*> latches;
  SNLBooleanTree::StateIdentifiers identifiers;
  std::string sharedEnable;
  for (const auto* group : cell->children) {
    if (group->id != "latch") continue;
    for (const auto* attribute : group->children) {
      if (attribute->id != "enable" && attribute->id != "data_in" &&
          attribute->id != "clear" && attribute->id != "preset" &&
          attribute->id != "clear_preset_var1" && attribute->id != "clear_preset_var2") {
        throw std::runtime_error("unsupported latch attribute `" + attribute->id + "`");
      }
    }
    if (group->args.empty() || group->args.size() > 2) {
      throw std::runtime_error("latch requires one or two state identifiers");
    }
    if (group->args.size() == 1 && child(group, "clear_preset_var2")) {
      throw std::runtime_error("clear_preset_var2 requires a second state identifier");
    }
    const auto& enable = requireChild(group, "enable")->value;
    requireChild(group, "data_in");
    if (!latches.empty() && enable != sharedEnable) {
      throw std::runtime_error("multiple latch groups require identical enable expressions");
    }
    sharedEnable = enable;
    for (size_t i = 0; i < group->args.size(); ++i) {
      if (group->args[i].empty() || primitive->getScalarTerm(NLName(group->args[i])) ||
          !identifiers.emplace(group->args[i],
              SNLBooleanTree::StateIdentifier{latches.size(), i != 0}).second) {
        throw std::runtime_error("ambiguous or duplicate latch state identifier");
      }
    }
    latches.push_back(group);
  }
  if (latches.empty()) return;

  Modeling::SequentialModel model;
  model.kind = Modeling::SequentialModel::Kind::Latch;
  model.clockedOn = expression(primitive, sharedEnable, identifiers);
  for (const auto* latch : latches) {
    Modeling::SequentialState state;
    state.nextState = expression(primitive, requireChild(latch, "data_in")->value, identifiers);
    if (const auto* clear = child(latch, "clear")) {
      state.clear = expression(primitive, clear->value, identifiers);
    }
    if (const auto* preset = child(latch, "preset")) {
      state.preset = expression(primitive, preset->value, identifiers);
    }
    state.clearPresetValue = conflictValue(latch, "clear_preset_var1");
    if (latch->args.size() == 2 &&
        conflictValue(latch, "clear_preset_var2") != complement(state.clearPresetValue)) {
      // SequentialModel represents the second variable as the complement of
      // the first; Liberty can instead specify independently forced values.
      throw std::runtime_error("clear/preset conflict does not preserve complementary state outputs");
    }
    model.states.push_back(std::move(state));
  }
  for (auto* term : primitive->getBitTerms()) {
    if (!dynamic_cast<SNLScalarTerm*>(term) || term->getDirection() == SNLTerm::Direction::InOut) {
      throw std::runtime_error("only scalar input/output latch pins are supported");
    }
    if (term->getDirection() != SNLTerm::Direction::Output) continue;
    const Ast* pin = nullptr;
    for (const auto* candidate : cell->children) {
      if (candidate->id == "pin" && candidate->args.size() == 1 &&
          candidate->args[0] == term->getName().getString()) {
        if (pin) throw std::runtime_error("duplicate output pin definition");
        pin = candidate;
      }
    }
    if (!pin || child(pin, "three_state")) {
      throw std::runtime_error("missing output pin or unsupported tri-state output");
    }
    model.outputs.push_back({term,
        expression(primitive, requireChild(pin, "function")->value, identifiers)});
  }
  if (!model.isValid()) throw std::runtime_error("incomplete latch model");
  Modeling::setSequentialModel(primitive, model);
  addArcs(model);
}

std::string readText(const std::filesystem::path& path) {
  // gzopen accepts uncompressed streams as well, so gzip and ordinary files
  // follow the same parser without choosing behavior from the filename.
  const auto close = [](gzFile file) { if (file) gzclose(file); };
  std::unique_ptr<gzFile_s, decltype(close)> input(gzopen(path.string().c_str(), "rb"), close);
  if (!input) throw std::runtime_error("cannot open Liberty input");
  std::string contents;
  std::array<char, 65536> buffer;
  int count;
  while ((count = gzread(input.get(), buffer.data(), buffer.size())) > 0) {
    contents.append(buffer.data(), static_cast<size_t>(count));
  }
  if (count < 0) throw std::runtime_error("cannot decompress Liberty input");
  if (contents.starts_with("PK\003\004")) {
    throw std::runtime_error("supplemental latch models do not support ZIP archives");
  }
  return contents;
}

}  // namespace

void constructLibertyWithLatchModels(NLLibrary* library,
                                     const std::filesystem::path& path) {
  if (!library) throw std::invalid_argument("Liberty latch loading requires a library");
  std::unordered_set<std::string> seen;
  for (const auto* design : library->getSNLDesigns()) {
    seen.insert(design->getName().getString());
  }
  SNLLibertyConstructor(library).construct(path);
  try {
    std::istringstream input(readText(path));
    // Full parsing retains clear_preset_var2 and unsupported behavior markers
    // which the frontend's structural parse may omit.
    Yosys::LibertyParser parser(input);
    if (!parser.ast || parser.ast->id != "library") {
      throw std::runtime_error("expected a Liberty library group");
    }
    for (const auto* cell : parser.ast->children) {
      if (cell->id != "cell" || cell->args.empty() ||
          !seen.insert(cell->args.front()).second) continue;
      auto* primitive = library->getSNLDesign(NLName(cell->args.front()));
      if (!primitive || !contains(cell, "latch") || Modeling::hasSequentialModel(primitive)) continue;
      try {
        populate(primitive, cell);
      } catch (const std::exception& error) {
        SPDLOG_WARN("Liberty latch cell `{}` in {} remains opaque: {}",
                    cell->args.front(), path.string(), error.what());
      }
    }
  } catch (const std::exception& error) {
    SPDLOG_WARN("Supplemental Liberty latch modeling unavailable for {}: {}",
                path.string(), error.what());
  }
}

}  // namespace KEPLER_FORMAL
