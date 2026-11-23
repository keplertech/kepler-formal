// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "ClockGateModels.h"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

#include "NLDB.h"
#include "NLLibrary.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLBitNet.h"
#include "SNLBooleanTree.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstTerm.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"

namespace KEPLER_FORMAL {

namespace {

using naja::NL::NLLibrary;
using naja::NL::NLName;
using naja::NL::SNLBitTerm;
using naja::NL::SNLBooleanTree;
using naja::NL::SNLDesign;
using naja::NL::SNLDesignModeling;
using naja::NL::SNLInstance;
using naja::NL::SNLScalarNet;
using naja::NL::SNLScalarTerm;
using naja::NL::SNLTerm;
using naja::NL::SNLTruthTable;

constexpr const char* kClockGateAttribute = "clock_gating_integrated_cell";
constexpr const char* kLatchSuffix = "__kf_icg_latch";
constexpr const char* kGateSuffix = "__kf_icg_gate";
constexpr const char* kStateNetSuffix = "__kf_icg_state";
// A statetable with more inputs than this is not an ICG worth enumerating.
constexpr size_t kMaxStatetableInputs = 8;
// Instances left opaque are listed one by one up to this many.
constexpr size_t kMaxReportedKeptInstances = 10;

// ---------------------------------------------------------------------------
// Liberty reader. It keeps only the `cell` groups that carry a given
// attribute and steps over everything else without building anything, so the
// timing tables that make up most of a Liberty file cost one pass over the
// text.
// ---------------------------------------------------------------------------

struct LibertyStatement {
  std::string name;
  std::vector<std::string> args;      // group or complex attribute arguments
  std::string value;                  // simple attribute value
  std::vector<LibertyStatement> body; // statements of a group
  bool isGroup = false;
};

class LibertyReader {
 public:
  explicit LibertyReader(std::string_view text) : text_(text) {}

  std::vector<LibertyStatement> readCellsWith(std::string_view attribute) {
    std::vector<LibertyStatement> cells;
    skipBlanks();
    while (!atEnd()) {
      const std::string name = readWord();
      if (name.empty()) {
        ++pos_;
      } else {
        skipBlanks();
        if (peek() == '(') {
          skipParentheses();
          skipBlanks();
        }
        if (peek() == '{') {
          ++pos_;
          if (name == "library") {
            readLibraryBody(attribute, cells);
          } else {
            skipBody();
          }
        } else {
          skipToStatementEnd();
        }
      }
      skipBlanks();
    }
    return cells;
  }

 private:
  bool atEnd() const { return pos_ >= text_.size(); }
  char peek(size_t ahead = 0) const {
    return pos_ + ahead < text_.size() ? text_[pos_ + ahead] : '\0';
  }
  static bool isSpace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' ||
           c == '\v';
  }
  static bool isPunctuation(char c) {
    return c == '(' || c == ')' || c == '{' || c == '}' || c == ':' ||
           c == ';' || c == ',' || c == '"';
  }
  bool atLineContinuation() const {
    return peek() == '\\' &&
           (peek(1) == '\n' || (peek(1) == '\r' && peek(2) == '\n'));
  }
  void skipLineContinuation() { pos_ += peek(1) == '\r' ? 3 : 2; }
  bool atComment() const {
    return peek() == '/' && (peek(1) == '*' || peek(1) == '/');
  }
  void skipComment() {
    if (peek(1) == '*') {
      const size_t end = text_.find("*/", pos_ + 2);
      pos_ = end == std::string_view::npos ? text_.size() : end + 2;
    } else {
      const size_t end = text_.find('\n', pos_);
      pos_ = end == std::string_view::npos ? text_.size() : end + 1;
    }
  }
  void skipBlanks() {
    while (!atEnd()) {
      if (isSpace(peek())) {
        ++pos_;
      } else if (atLineContinuation()) {
        skipLineContinuation();
      } else if (atComment()) {
        skipComment();
      } else {
        return;
      }
    }
  }

  std::string readWord() {
    const size_t begin = pos_;
    while (!atEnd() && !isSpace(peek()) && !isPunctuation(peek()) &&
           !atComment() && !atLineContinuation()) {
      ++pos_;
    }
    return std::string(text_.substr(begin, pos_ - begin));
  }

  // A quoted string; a backslash before a line break continues it.
  std::string readString() {
    std::string value;
    ++pos_;  // opening quote
    while (!atEnd()) {
      if (atLineContinuation()) {
        skipLineContinuation();
      } else if (peek() == '\\' && peek(1) == '"') {
        value.push_back('"');
        pos_ += 2;
      } else if (peek() == '"') {
        ++pos_;
        break;
      } else {
        value.push_back(peek());
        ++pos_;
      }
    }
    return value;
  }
  void skipString() {
    ++pos_;  // opening quote
    while (!atEnd()) {
      if (peek() == '\\' && peek(1) == '"') {
        pos_ += 2;
      } else if (peek() == '"') {
        ++pos_;
        return;
      } else {
        ++pos_;
      }
    }
  }

  std::vector<std::string> readArguments() {
    std::vector<std::string> args;
    ++pos_;  // opening parenthesis
    while (!atEnd()) {
      skipBlanks();
      if (peek() == ')') {
        ++pos_;
        break;
      }
      if (peek() == ',') {
        ++pos_;
      } else if (peek() == '"') {
        args.push_back(readString());
      } else {
        std::string word = readWord();
        if (word.empty()) {
          ++pos_;  // a character that belongs to no argument
        } else {
          args.push_back(std::move(word));
        }
      }
    }
    return args;
  }
  void skipParentheses() {
    ++pos_;  // opening parenthesis
    while (!atEnd()) {
      if (peek() == '"') {
        skipString();
      } else if (peek() == ')') {
        ++pos_;
        return;
      } else {
        ++pos_;
      }
    }
  }

  // The value of a simple attribute: a string, or the text up to the end of
  // the statement.
  std::string readValue() {
    while (!atEnd() && (peek() == ' ' || peek() == '\t' ||
                        atLineContinuation())) {
      if (atLineContinuation()) {
        skipLineContinuation();
      } else {
        ++pos_;
      }
    }
    if (peek() == '"') {
      std::string value = readString();
      while (!atEnd() && (peek() == ' ' || peek() == '\t')) {
        ++pos_;
      }
      return value;
    }
    const size_t begin = pos_;
    while (!atEnd() && peek() != ';' && peek() != '\n' && peek() != '}' &&
           !atComment()) {
      ++pos_;
    }
    size_t end = pos_;
    while (end > begin && isSpace(text_[end - 1])) {
      --end;
    }
    return std::string(text_.substr(begin, end - begin));
  }

  void skipToStatementEnd() {
    while (!atEnd() && peek() != '}') {
      if (peek() == '"') {
        skipString();
      } else if (peek() == ';' || peek() == '\n') {
        ++pos_;
        return;
      } else {
        ++pos_;
      }
    }
  }

  // Steps over a group body whose opening brace has been read.
  void skipBody() {
    size_t depth = 1;
    while (!atEnd()) {
      const char c = peek();
      if (c == '"') {
        skipString();
      } else if (atComment()) {
        skipComment();
      } else {
        ++pos_;
        if (c == '{') {
          ++depth;
        } else if (c == '}' && --depth == 0) {
          return;
        }
      }
    }
  }

  void parseStatement(LibertyStatement& statement) {
    skipBlanks();
    statement.name = readWord();
    skipBlanks();
    if (peek() == '(') {
      statement.args = readArguments();
      skipBlanks();
      if (peek() == '{') {
        ++pos_;
        statement.isGroup = true;
        while (true) {
          skipBlanks();
          if (atEnd()) {
            break;
          }
          if (peek() == '}') {
            ++pos_;
            break;
          }
          const size_t before = pos_;
          LibertyStatement child;
          parseStatement(child);
          if (pos_ == before) {
            ++pos_;  // a character that starts no statement
          } else if (!child.name.empty()) {
            statement.body.push_back(std::move(child));
          }
        }
      } else if (peek() == ';') {
        ++pos_;
      }
    } else if (peek() == ':') {
      ++pos_;
      statement.value = readValue();
      if (peek() == ';') {
        ++pos_;
      }
    } else {
      skipToStatementEnd();
    }
  }

  void readLibraryBody(std::string_view attribute,
                       std::vector<LibertyStatement>& cells) {
    while (true) {
      skipBlanks();
      if (atEnd()) {
        return;
      }
      if (peek() == '}') {
        ++pos_;
        return;
      }
      const size_t begin = pos_;
      const std::string name = readWord();
      if (name.empty()) {
        ++pos_;
        continue;
      }
      skipBlanks();
      if (peek() == '(') {
        skipParentheses();
        skipBlanks();
        if (peek() == '{') {
          ++pos_;
          skipBody();
          const std::string_view group = text_.substr(begin, pos_ - begin);
          if (name == "cell" &&
              group.find(attribute) != std::string_view::npos) {
            LibertyReader groupReader(group);
            LibertyStatement cell;
            groupReader.parseStatement(cell);
            cells.push_back(std::move(cell));
          }
        } else if (peek() == ';') {
          ++pos_;
        }
      } else {
        skipToStatementEnd();
      }
    }
  }

  std::string_view text_;
  size_t pos_ = 0;
};

const LibertyStatement* findChild(const LibertyStatement& group,
                                  std::string_view name) {
  for (const auto& child : group.body) {
    if (child.name == name) {
      return &child;
    }
  }
  return nullptr;
}

std::string attributeValue(const LibertyStatement& group,
                           std::string_view name) {
  const auto* child = findChild(group, name);
  return child != nullptr && !child->isGroup ? child->value : std::string();
}

std::vector<std::string> splitWords(const std::string& text) {
  std::istringstream stream(text);
  return {std::istream_iterator<std::string>(stream),
          std::istream_iterator<std::string>()};
}

std::vector<std::string> splitFields(const std::string& text, char separator) {
  std::vector<std::string> fields;
  std::string field;
  std::istringstream stream(text);
  while (std::getline(stream, field, separator)) {
    fields.push_back(field);
  }
  return fields;
}

// ---------------------------------------------------------------------------
// From a Liberty cell to the latch and gate it stands for.
// ---------------------------------------------------------------------------

struct ClockGateOutput {
  std::string pin;
  std::string function;
};

struct ClockGateCell {
  std::string name;
  // Pin that carries the latch state in the output functions.
  std::string statePin;
  // Pin whose level opens the latch, and that level.
  std::string enablePin;
  bool transparentHigh = false;
  // Pins the latch data depends on, and its truth table over them: entry `i`
  // is the data when pin `k` has the value of bit `k` of `i`.
  std::vector<std::string> dataPins;
  std::vector<bool> data;
  std::vector<std::string> inputPins;
  std::vector<ClockGateOutput> outputs;
};

enum class Level { Low, High, Any };

struct StatetableRow {
  std::vector<Level> inputs;
  Level current = Level::Any;
  // 0 or 1: the latch takes that value. -1: it holds.
  int next = -1;
};

bool parseLevel(const std::string& token, Level& level) {
  if (token == "L") {
    level = Level::Low;
  } else if (token == "H") {
    level = Level::High;
  } else if (token == "-") {
    level = Level::Any;
  } else {
    return false;
  }
  return true;
}

bool levelMatches(Level level, bool value) {
  return level == Level::Any || (level == Level::High) == value;
}

// The rows of a statetable with one internal node, or the reason it is not a
// plain level-sensitive table: edge-triggered rows (R, F), unknown next
// states (X) and tables with several nodes are not latches.
std::optional<std::vector<StatetableRow>> parseStatetableRows(
    const std::string& table, size_t inputCount, std::string& whyNot) {
  std::vector<StatetableRow> rows;
  for (const auto& rowText : splitFields(table, ',')) {
    if (splitWords(rowText).empty()) {
      continue;
    }
    const auto fields = splitFields(rowText, ':');
    if (fields.size() != 3) {
      whyNot = "statetable row `" + rowText + "` does not have three fields";
      return std::nullopt;
    }
    const auto inputs = splitWords(fields[0]);
    const auto current = splitWords(fields[1]);
    const auto next = splitWords(fields[2]);
    if (inputs.size() != inputCount || current.size() != 1 ||
        next.size() != 1) {
      whyNot = "statetable row `" + rowText +
               "` does not match the declared inputs and internal node";
      return std::nullopt;
    }
    StatetableRow row;
    for (const auto& token : inputs) {
      Level level = Level::Any;
      if (!parseLevel(token, level)) {
        whyNot = "statetable input value `" + token +
                 "` is not a level (L, H or -)";
        return std::nullopt;
      }
      row.inputs.push_back(level);
    }
    if (!parseLevel(current[0], row.current)) {
      whyNot = "statetable current state `" + current[0] + "` is not a level";
      return std::nullopt;
    }
    if (next[0] == "L") {
      row.next = 0;
    } else if (next[0] == "H") {
      row.next = 1;
    } else if (next[0] == "N") {
      row.next = -1;
    } else {
      whyNot = "statetable next state `" + next[0] + "` is not L, H or N";
      return std::nullopt;
    }
    rows.push_back(std::move(row));
  }
  if (rows.empty()) {
    whyNot = "statetable has no rows";
    return std::nullopt;
  }
  return rows;
}

// Reads the latch out of the statetable rows: the one input whose level opens
// it, and the data it takes while open.
bool describeLatch(const std::vector<StatetableRow>& rows,
                   const std::vector<std::string>& inputs,
                   ClockGateCell& cell,
                   std::string& whyNot) {
  const size_t inputCount = inputs.size();
  const size_t assignments = size_t{1} << inputCount;
  // Per input assignment: -1 hold, 0 or 1 the value taken.
  std::vector<int> next(assignments, -1);
  for (size_t assignment = 0; assignment < assignments; ++assignment) {
    for (const bool current : {false, true}) {
      std::optional<int> result;
      for (const auto& row : rows) {
        bool matches = levelMatches(row.current, current);
        for (size_t i = 0; matches && i < inputCount; ++i) {
          matches = levelMatches(row.inputs[i], ((assignment >> i) & 1) != 0);
        }
        if (!matches) {
          continue;
        }
        if (result.has_value() && *result != row.next) {
          whyNot = "statetable rows disagree for one input combination";
          return false;
        }
        result = row.next;
      }
      if (!result.has_value()) {
        whyNot = "statetable does not cover every input combination";
        return false;
      }
      if (current && *result != next[assignment]) {
        whyNot = "statetable next state depends on the current state";
        return false;
      }
      next[assignment] = *result;
    }
  }

  // The latch is open exactly while one input is at one level.
  std::optional<size_t> enableIndex;
  for (size_t i = 0; i < inputCount && !enableIndex.has_value(); ++i) {
    for (const bool level : {false, true}) {
      bool isEnable = true;
      for (size_t assignment = 0; isEnable && assignment < assignments;
           ++assignment) {
        const bool open = next[assignment] != -1;
        isEnable = open == ((((assignment >> i) & 1) != 0) == level);
      }
      if (isEnable) {
        enableIndex = i;
        cell.transparentHigh = level;
        break;
      }
    }
  }
  if (!enableIndex.has_value()) {
    whyNot = "statetable is not a latch opened by the level of one input";
    return false;
  }
  if (cell.transparentHigh) {
    // The latch of a negative-edge clock gate. SEC keeps a latch that is open
    // while the clock is high as state captured at the falling edge, which a
    // flop clocked on that same edge must not read.
    whyNot = "its latch is open while `" + inputs[*enableIndex] +
             "` is high; only clock gates for rising-edge flops are modeled";
    return false;
  }
  cell.enablePin = inputs[*enableIndex];

  // Data over the other inputs, taken where the latch is open.
  std::vector<size_t> others;
  for (size_t i = 0; i < inputCount; ++i) {
    if (i != *enableIndex) {
      others.push_back(i);
    }
  }
  const auto dataAt = [&](size_t otherValues) {
    size_t assignment = cell.transparentHigh ? size_t{1} << *enableIndex : 0;
    for (size_t k = 0; k < others.size(); ++k) {
      if ((otherValues >> k) & 1) {
        assignment |= size_t{1} << others[k];
      }
    }
    return next[assignment] == 1;
  };
  // Keep only the inputs the data depends on.
  std::vector<size_t> used;
  for (size_t k = 0; k < others.size(); ++k) {
    for (size_t values = 0; values < (size_t{1} << others.size()); ++values) {
      if (dataAt(values) != dataAt(values ^ (size_t{1} << k))) {
        used.push_back(k);
        break;
      }
    }
  }
  cell.dataPins.clear();
  for (const size_t k : used) {
    cell.dataPins.push_back(inputs[others[k]]);
  }
  cell.data.assign(size_t{1} << used.size(), false);
  for (size_t values = 0; values < cell.data.size(); ++values) {
    size_t otherValues = 0;
    for (size_t j = 0; j < used.size(); ++j) {
      if ((values >> j) & 1) {
        otherValues |= size_t{1} << used[j];
      }
    }
    cell.data[values] = dataAt(otherValues);
  }
  return true;
}

std::optional<ClockGateCell> describeClockGate(const LibertyStatement& cell,
                                               std::string& whyNot) {
  ClockGateCell gate;
  if (cell.args.empty()) {
    whyNot = "cell without a name";
    return std::nullopt;
  }
  gate.name = cell.args[0];

  const LibertyStatement* statetable = nullptr;
  for (const auto& child : cell.body) {
    if (child.name == "bus" || child.name == "bundle") {
      whyNot = "bus and bundle pins are not supported";
      return std::nullopt;
    }
    if (child.name == "statetable") {
      if (statetable != nullptr) {
        whyNot = "more than one statetable";
        return std::nullopt;
      }
      statetable = &child;
    }
  }
  if (statetable == nullptr) {
    whyNot = "no statetable";
    return std::nullopt;
  }
  if (statetable->args.size() != 2) {
    whyNot = "statetable does not name its inputs and its internal node";
    return std::nullopt;
  }
  const auto tableInputs = splitWords(statetable->args[0]);
  const auto tableNodes = splitWords(statetable->args[1]);
  if (tableNodes.size() != 1) {
    whyNot = "statetable has " + std::to_string(tableNodes.size()) +
             " internal nodes; one latch is supported";
    return std::nullopt;
  }
  if (tableInputs.empty() || tableInputs.size() > kMaxStatetableInputs) {
    whyNot = "statetable has " + std::to_string(tableInputs.size()) +
             " inputs";
    return std::nullopt;
  }

  gate.statePin = tableNodes[0];
  for (const auto& child : cell.body) {
    if (child.name != "pin" || child.args.empty()) {
      continue;
    }
    const std::string direction = attributeValue(child, "direction");
    if (direction == "input") {
      gate.inputPins.push_back(child.args[0]);
    } else if (direction == "output") {
      std::string function = attributeValue(child, "state_function");
      if (function.empty()) {
        function = attributeValue(child, "function");
      }
      if (function.empty()) {
        whyNot = "output pin `" + child.args[0] + "` has no function";
        return std::nullopt;
      }
      gate.outputs.push_back({child.args[0], std::move(function)});
    } else if (direction == "internal") {
      if (attributeValue(child, "internal_node") == tableNodes[0]) {
        gate.statePin = child.args[0];
      }
    } else {
      whyNot = "pin `" + child.args[0] + "` has direction `" + direction + "`";
      return std::nullopt;
    }
  }
  if (gate.outputs.empty()) {
    whyNot = "no output pin";
    return std::nullopt;
  }
  for (const auto& input : tableInputs) {
    if (std::find(gate.inputPins.begin(), gate.inputPins.end(), input) ==
        gate.inputPins.end()) {
      whyNot = "statetable input `" + input + "` is not an input pin";
      return std::nullopt;
    }
  }
  if (std::find(gate.inputPins.begin(), gate.inputPins.end(), gate.statePin) !=
      gate.inputPins.end()) {
    whyNot = "internal node `" + gate.statePin + "` is also an input pin";
    return std::nullopt;
  }

  const auto rows = parseStatetableRows(attributeValue(*statetable, "table"),
                                        tableInputs.size(), whyNot);
  if (!rows.has_value() || !describeLatch(*rows, tableInputs, gate, whyNot)) {
    return std::nullopt;
  }
  return gate;
}

// ---------------------------------------------------------------------------
// Naja primitives for the latch and the gate.
// ---------------------------------------------------------------------------

using Expression = SNLDesignModeling::BooleanExpression;

// Shannon expansion of `table` over `terms`, last term first. Entry `i` of
// the table is the value when term `k` has the value of bit `k` of `i`.
Expression::NodeID addFunction(Expression& expression,
                               const std::vector<SNLBitTerm*>& terms,
                               const std::vector<bool>& table,
                               size_t begin,
                               size_t termCount) {
  const size_t size = size_t{1} << termCount;
  const auto all = [&](size_t from, size_t count, bool value) {
    return std::all_of(table.begin() + from, table.begin() + from + count,
                       [value](bool entry) { return entry == value; });
  };
  if (all(begin, size, false)) {
    return expression.addConstant(false);
  }
  if (all(begin, size, true)) {
    return expression.addConstant(true);
  }
  const size_t half = size / 2;
  const size_t low = begin;
  const size_t high = begin + half;
  if (std::equal(table.begin() + low, table.begin() + low + half,
                 table.begin() + high)) {
    return addFunction(expression, terms, table, low, termCount - 1);
  }
  using Operator = Expression::Operator;
  const auto term = expression.addTerm(terms[termCount - 1]);
  const auto notTerm = [&]() {
    return expression.addOperation(Operator::Not, {term});
  };
  const auto cofactor = [&](size_t from) {
    return addFunction(expression, terms, table, from, termCount - 1);
  };
  if (all(low, half, false)) {
    return all(high, half, true)
               ? term
               : expression.addOperation(Operator::And, {term, cofactor(high)});
  }
  if (all(high, half, false)) {
    const auto inverted = notTerm();
    return all(low, half, true)
               ? inverted
               : expression.addOperation(Operator::And,
                                         {inverted, cofactor(low)});
  }
  if (all(high, half, true)) {
    return expression.addOperation(Operator::Or, {term, cofactor(low)});
  }
  if (all(low, half, true)) {
    return expression.addOperation(Operator::Or, {notTerm(), cofactor(high)});
  }
  const auto whenHigh =
      expression.addOperation(Operator::And, {term, cofactor(high)});
  const auto whenLow =
      expression.addOperation(Operator::And, {notTerm(), cofactor(low)});
  return expression.addOperation(Operator::Or, {whenHigh, whenLow});
}

// The latch: open while `enablePin` is at its transparent level, output on a
// pin named after the internal node.
SNLDesign* createLatchPrimitive(NLLibrary* library, const ClockGateCell& cell) {
  auto* latch = SNLDesign::create(library, SNLDesign::Type::Primitive,
                                  NLName(cell.name + kLatchSuffix));
  auto* enable = SNLScalarTerm::create(latch, SNLTerm::Direction::Input,
                                       NLName(cell.enablePin));
  std::vector<SNLBitTerm*> dataTerms;
  for (const auto& pin : cell.dataPins) {
    dataTerms.push_back(SNLScalarTerm::create(
        latch, SNLTerm::Direction::Input, NLName(pin)));
  }
  auto* state = SNLScalarTerm::create(latch, SNLTerm::Direction::Output,
                                      NLName(cell.statePin));

  SNLDesignModeling::SequentialModel model;
  model.kind = SNLDesignModeling::SequentialModel::Kind::Latch;
  auto enableNode = model.clockedOn.addTerm(enable);
  if (!cell.transparentHigh) {
    enableNode = model.clockedOn.addOperation(Expression::Operator::Not,
                                              {enableNode});
  }
  model.clockedOn.root = enableNode;
  SNLDesignModeling::SequentialState latchState;
  latchState.nextState.root = addFunction(
      latchState.nextState, dataTerms, cell.data, 0, dataTerms.size());
  model.states.push_back(std::move(latchState));
  Expression stateOutput;
  stateOutput.root = stateOutput.addState(0);
  model.outputs.push_back({state, std::move(stateOutput)});
  SNLDesignModeling::setSequentialModel(latch, model);

  // The roles and timing arcs Naja derives from a Liberty `latch` group.
  SNLDesignModeling::setTermRole(enable,
                                 SNLDesignModeling::SNLTermRole::Clock);
  SNLDesignModeling::BitTerms updateTerms;
  for (auto* term : dataTerms) {
    updateTerms.push_back(term);
  }
  if (!updateTerms.empty()) {
    SNLDesignModeling::addInputsToClockArcs(updateTerms, enable);
  }
  SNLDesignModeling::setTermRole(state,
                                 SNLDesignModeling::SNLTermRole::DataOutput);
  SNLDesignModeling::addClockToOutputsArcs(enable, {state});
  return latch;
}

// Reading the latch as a wire is right for whatever samples the output when
// the clock rises, the instant the latch closes. That needs two properties of
// an output function, checked here over all its inputs:
// - while the clock is low and the latch open, the output does not follow the
//   latch, so the enable cannot move it between clock edges;
// - the output never falls when the clock rises, so its rising edges are
//   rising edges of the clock.
// `CLK & IQ` has both. An inverted gated clock, or a latch output passed
// through, does not; the cell then keeps no model.
void requireGatedClockOutput(const SNLBooleanTree& tree,
                             const SNLDesign* pins,
                             const ClockGateCell& cell,
                             const ClockGateOutput& output) {
  std::vector<naja::NL::SNLBooleanTreeInputNode*> nodes;
  std::optional<size_t> clockIndex;
  std::optional<size_t> stateIndex;
  for (const auto& [term, node] : tree.getInputs()) {
    if (term == pins->getScalarTerm(NLName(cell.enablePin))) {
      clockIndex = nodes.size();
    } else if (term == pins->getScalarTerm(NLName(cell.statePin))) {
      stateIndex = nodes.size();
    }
    nodes.push_back(node);
  }
  if (!clockIndex.has_value() || tree.getRoot() == nullptr) {
    return;  // the output does not read the clock: no edge comes from it
  }
  const auto valueAt = [&](size_t assignment) {
    for (size_t i = 0; i < nodes.size(); ++i) {
      nodes[i]->setValue(((assignment >> i) & 1) != 0);
    }
    return tree.getRoot()->getValue();
  };
  const size_t clockBit = size_t{1} << *clockIndex;
  for (size_t assignment = 0; assignment < (size_t{1} << nodes.size());
       ++assignment) {
    if ((assignment & clockBit) != 0) {
      continue;
    }
    const bool clockLow = valueAt(assignment);
    if (clockLow && !valueAt(assignment | clockBit)) {
      throw std::runtime_error("output `" + output.pin +
                               "` falls when `" + cell.enablePin +
                               "` rises: not a gated copy of the clock");
    }
    if (stateIndex.has_value() &&
        clockLow != valueAt(assignment ^ (size_t{1} << *stateIndex))) {
      throw std::runtime_error("output `" + output.pin +
                               "` follows the latch while `" + cell.enablePin +
                               "` is low: not a gated copy of the clock");
    }
  }
}

// The gate: the output functions of the cell, with the internal node as one
// more input. Truth tables are built the way Naja builds them for a Liberty
// `function`.
SNLDesign* createGatePrimitive(NLLibrary* library, const ClockGateCell& cell) {
  std::vector<std::string> candidates = cell.inputPins;
  candidates.push_back(cell.statePin);

  // Which pins the output functions read.
  std::set<std::string> usedPins;
  {
    auto* scratch = SNLDesign::create(library, SNLDesign::Type::Primitive,
                                      NLName(cell.name + "__kf_icg_scratch"));
    try {
      for (const auto& pin : candidates) {
        SNLScalarTerm::create(scratch, SNLTerm::Direction::Input, NLName(pin));
      }
      for (const auto& output : cell.outputs) {
        SNLBooleanTree tree;
        tree.parse(scratch, output.function);
        for (const auto& [term, _] : tree.getInputs()) {
          usedPins.insert(term->getName().getString());
        }
        requireGatedClockOutput(tree, scratch, cell, output);
      }
    } catch (...) {
      scratch->destroy();
      throw;
    }
    scratch->destroy();
  }

  auto* gate = SNLDesign::create(library, SNLDesign::Type::Primitive,
                                 NLName(cell.name + kGateSuffix));
  try {
    for (const auto& pin : candidates) {
      if (usedPins.count(pin) != 0) {
        SNLScalarTerm::create(gate, SNLTerm::Direction::Input, NLName(pin));
      }
    }
    for (const auto& output : cell.outputs) {
      SNLScalarTerm::create(gate, SNLTerm::Direction::Output,
                            NLName(output.pin));
    }
    SNLBooleanTree::Terms inputs;
    for (auto* term : gate->getBitTerms()) {
      if (term->getDirection() == SNLTerm::Direction::Input) {
        inputs.push_back(term);
      }
    }
    std::reverse(inputs.begin(), inputs.end());
    std::vector<SNLTruthTable> truthTables;
    for (const auto& output : cell.outputs) {
      SNLBooleanTree tree;
      tree.parse(gate, output.function);
      truthTables.push_back(tree.getTruthTable(inputs));
    }
    if (truthTables.size() == 1) {
      SNLDesignModeling::setTruthTable(gate, truthTables.front());
    } else {
      SNLDesignModeling::setTruthTables(gate, truthTables);
    }
  } catch (...) {
    gate->destroy();
    throw;
  }
  return gate;
}

// True when `primitive` has the scalar pins the model reads and drives.
bool hasModelPins(const SNLDesign* primitive, const ClockGateCell& cell) {
  const auto hasPin = [&](const std::string& pin,
                          SNLTerm::Direction direction) {
    auto* term = primitive->getScalarTerm(NLName(pin));
    return term != nullptr && term->getDirection() == direction;
  };
  if (!hasPin(cell.enablePin, SNLTerm::Direction::Input)) {
    return false;
  }
  for (const auto& pin : cell.inputPins) {
    if (!hasPin(pin, SNLTerm::Direction::Input)) {
      return false;
    }
  }
  for (const auto& output : cell.outputs) {
    if (!hasPin(output.pin, SNLTerm::Direction::Output)) {
      return false;
    }
  }
  return true;
}

std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return {};
  }
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

bool isCompressed(const std::string& text) {
  return text.size() >= 2 &&
         ((static_cast<unsigned char>(text[0]) == 0x1f &&
           static_cast<unsigned char>(text[1]) == 0x8b) ||
          (text[0] == 'P' && text[1] == 'K'));
}

NLName unusedInstanceName(const SNLDesign* design, const std::string& stem) {
  std::string name = stem;
  for (size_t suffix = 1; design->getInstance(NLName(name)) != nullptr;
       ++suffix) {
    name = stem + "_" + std::to_string(suffix);
  }
  return NLName(name);
}

NLName unusedNetName(const SNLDesign* design, const std::string& stem) {
  std::string name = stem;
  for (size_t suffix = 1; design->getNet(NLName(name)) != nullptr; ++suffix) {
    name = stem + "_" + std::to_string(suffix);
  }
  return NLName(name);
}

// The latch of the model is open while the clock is low, which SEC reduces to
// a wire. That is what a flop capturing on the rising edge of the gated clock
// sees; a flop on the falling edge, a latch or any other logic would read the
// value the latch held through the high phase instead. So an instance is
// replaced only when its outputs reach, through buffers and inverters of the
// same module, nothing but clock pins of flops that capture on the rising
// edge of the gated clock. Returns the reason when that is not the case.
std::optional<std::string> whyNotOnlyRisingEdgeFlops(
    const SNLInstance* instance,
    const SNLDesign* gate) {
  using Operator = Expression::Operator;
  const auto* cell = instance->getModel();
  // A net with the polarity it carries relative to the gated clock.
  std::vector<std::pair<const naja::NL::SNLBitNet*, bool>> pending;
  for (auto* term : gate->getScalarTerms()) {
    if (term->getDirection() != SNLTerm::Direction::Output) {
      continue;
    }
    auto* cellTerm = cell->getScalarTerm(term->getName());
    auto* net = cellTerm == nullptr
                    ? nullptr
                    : instance->getInstTerm(cellTerm)->getNet();
    if (net != nullptr) {
      pending.emplace_back(net, false);
    }
  }
  std::set<std::pair<const naja::NL::SNLBitNet*, bool>> seen;
  // The clock gate and the buffers already walked through: their own output
  // pins are the expected drivers of the nets visited.
  std::set<const SNLInstance*> drivers{instance};
  while (!pending.empty()) {
    const auto [net, inverted] = pending.back();
    pending.pop_back();
    if (!seen.emplace(net, inverted).second) {
      continue;
    }
    if (!net->getBitTerms().empty()) {
      return "its gated clock reaches a port of module " +
             instance->getDesign()->getName().getString();
    }
    for (auto* instTerm : net->getInstTerms()) {
      auto* load = instTerm->getInstance();
      auto* pin = instTerm->getBitTerm();
      const bool isInput = pin->getDirection() == SNLTerm::Direction::Input;
      if (!isInput && drivers.count(load) != 0) {
        continue;
      }
      const std::string loadName =
          load->getName().getString() + " (" +
          load->getModel()->getName().getString() + ")";
      if (!isInput) {
        return "its gated clock has another driver, " + loadName;
      }
      auto* model = load->getModel();
      if (!model->isPrimitive()) {
        return "its gated clock enters module instance " + loadName;
      }
      if (SNLDesignModeling::hasSequentialModel(model)) {
        const auto& sequential = SNLDesignModeling::getSequentialModel(load);
        const auto& clock = sequential.clockedOn;
        std::optional<bool> risingEdge;
        if (sequential.kind ==
                SNLDesignModeling::SequentialModel::Kind::FlipFlop &&
            clock.isValid()) {
          const auto& root = clock.nodes[clock.root];
          if (root.operation == Operator::Term && root.term == pin) {
            risingEdge = true;
          } else if (root.operation == Operator::Not &&
                     root.operands.size() == 1 &&
                     clock.nodes[root.operands[0]].operation ==
                         Operator::Term &&
                     clock.nodes[root.operands[0]].term == pin) {
            risingEdge = false;
          }
        }
        if (!risingEdge.has_value()) {
          return "its gated clock reaches pin " + pin->getName().getString() +
                 " of " + loadName + ", which is not a flop clock pin";
        }
        if (*risingEdge == inverted) {
          return "flop " + loadName +
                 " captures on the falling edge of its gated clock";
        }
        continue;
      }
      const bool isBuffer = SNLDesignModeling::isBuf(model);
      if (!isBuffer && !SNLDesignModeling::isInv(model)) {
        return "its gated clock reaches " + loadName +
               ", which is neither a flop clock pin nor a buffer";
      }
      drivers.insert(load);
      for (auto* output : model->getScalarTerms()) {
        if (output->getDirection() == SNLTerm::Direction::Input) {
          continue;
        }
        if (auto* outputNet = load->getInstTerm(output)->getNet()) {
          pending.emplace_back(outputNet, isBuffer ? inverted : !inverted);
        }
      }
    }
  }
  return std::nullopt;
}

void replaceInstance(SNLInstance* instance,
                     SNLDesign* latch,
                     SNLDesign* gate) {
  auto* parent = instance->getDesign();
  auto* cell = instance->getModel();
  const std::string stem =
      instance->isUnnamed() ? "icg_" + std::to_string(instance->getID())
                            : instance->getName().getString();
  const auto netOf = [&](const NLName& pin) -> naja::NL::SNLNet* {
    auto* term = cell->getScalarTerm(pin);
    return term == nullptr ? nullptr : instance->getInstTerm(term)->getNet();
  };

  // The latch output is the one pin the cell itself does not have.
  NLName statePin;
  for (auto* term : latch->getScalarTerms()) {
    if (term->getDirection() == SNLTerm::Direction::Output) {
      statePin = term->getName();
    }
  }
  auto* stateNet = SNLScalarNet::create(
      parent, unusedNetName(parent, stem + kStateNetSuffix));
  auto* latchInstance = SNLInstance::create(
      parent, latch, unusedInstanceName(parent, stem + kLatchSuffix));
  for (auto* term : latch->getScalarTerms()) {
    auto* net = term->getName() == statePin
                    ? static_cast<naja::NL::SNLNet*>(stateNet)
                    : netOf(term->getName());
    if (net != nullptr) {
      latchInstance->getInstTerm(term)->setNet(net);
    }
  }
  auto* gateInstance = SNLInstance::create(
      parent, gate, unusedInstanceName(parent, stem + kGateSuffix));
  for (auto* term : gate->getScalarTerms()) {
    auto* net = term->getName() == statePin
                    ? static_cast<naja::NL::SNLNet*>(stateNet)
                    : netOf(term->getName());
    if (net != nullptr) {
      gateInstance->getInstTerm(term)->setNet(net);
    }
  }
  instance->destroy();
}

struct ClockGateParts {
  SNLDesign* latch = nullptr;
  SNLDesign* gate = nullptr;
};
using ClockGatePartsByCell = std::map<const SNLDesign*, ClockGateParts>;

// The modeled cells of `library` and its child libraries, and the designs
// that may instantiate them. Naja keeps no list of the instances of a
// primitive, so the instances are found from the designs.
void collectClockGates(NLLibrary* library,
                       ClockGatePartsByCell& partsByCell,
                       std::vector<SNLDesign*>& designs) {
  for (auto* design : library->getSNLDesigns()) {
    if (!library->isPrimitives()) {
      designs.push_back(design);
      continue;
    }
    if (design->isUnnamed()) {
      continue;
    }
    const std::string name = design->getName().getString();
    auto* latch = library->getSNLDesign(NLName(name + kLatchSuffix));
    auto* gate = library->getSNLDesign(NLName(name + kGateSuffix));
    if (latch != nullptr && gate != nullptr) {
      partsByCell.emplace(design, ClockGateParts{latch, gate});
    }
  }
  for (auto* child : library->getLibraries()) {
    collectClockGates(child, partsByCell, designs);
  }
}

}  // namespace

size_t addClockGateModels(NLLibrary* primitives,
                          const std::filesystem::path& libertyPath) {
  const std::string text = readFile(libertyPath);
  if (isCompressed(text)) {
    SPDLOG_INFO(
        "Clock-gating cells are read from uncompressed Liberty only; those "
        "of {} keep no model",
        libertyPath.string());
    return 0;
  }
  if (text.find(kClockGateAttribute) == std::string::npos) {
    return 0;
  }

  size_t modeled = 0;
  for (const auto& libertyCell :
       LibertyReader(text).readCellsWith(kClockGateAttribute)) {
    if (findChild(libertyCell, kClockGateAttribute) == nullptr) {
      continue;  // the attribute name only appears in a comment or a string
    }
    std::string whyNot;
    const auto cell = describeClockGate(libertyCell, whyNot);
    const std::string cellName =
        libertyCell.args.empty() ? std::string("?") : libertyCell.args[0];
    if (!cell.has_value()) {
      SPDLOG_INFO("Clock-gating cell {} is not modeled: {}", cellName, whyNot);
      continue;
    }
    auto* primitive = primitives->getSNLDesign(NLName(cell->name));
    if (primitive == nullptr || !primitive->isPrimitive() ||
        !hasModelPins(primitive, *cell)) {
      SPDLOG_INFO(
          "Clock-gating cell {} is not modeled: its loaded pins do not match "
          "the Liberty description",
          cellName);
      continue;
    }
    if (primitives->getSNLDesign(NLName(cell->name + kLatchSuffix)) !=
            nullptr ||
        primitives->getSNLDesign(NLName(cell->name + kGateSuffix)) !=
            nullptr) {
      continue;  // already modeled from an earlier library file
    }
    SNLDesign* latch = nullptr;
    try {
      latch = createLatchPrimitive(primitives, *cell);
      createGatePrimitive(primitives, *cell);
    } catch (const std::exception& e) {
      if (latch != nullptr) {
        latch->destroy();
      }
      SPDLOG_INFO("Clock-gating cell {} is not modeled: {}", cellName,
                  e.what());
      continue;
    }
    spdlog::debug(
        "Clock-gating cell {} is modeled as a latch open while {} is low and "
        "a gate on its output",
        cellName, cell->enablePin);
    ++modeled;
  }
  if (modeled > 0) {
    SPDLOG_INFO(
        "Modeled {} clock-gating cell(s) of {} as a latch and a gate",
        modeled, libertyPath.string());
  }
  return modeled;
}

size_t expandClockGates(naja::NL::NLUniverse* universe,
                        const std::set<std::string>& keptInstanceNames) {
  ClockGatePartsByCell partsByCell;
  std::vector<SNLDesign*> designs;
  for (auto* db : universe->getUserDBs()) {
    for (auto* library : db->getLibraries()) {
      collectClockGates(library, partsByCell, designs);
    }
  }
  if (partsByCell.empty()) {
    return 0;
  }
  size_t replaced = 0;
  size_t kept = 0;
  for (auto* design : designs) {
    std::vector<SNLInstance*> instances;
    for (auto* instance : design->getInstances()) {
      if (partsByCell.count(instance->getModel()) != 0) {
        instances.push_back(instance);
      }
    }
    for (auto* instance : instances) {
      const auto& parts = partsByCell.at(instance->getModel());
      const auto whyNot =
          keptInstanceNames.count(instance->getName().getString()) != 0
              ? std::optional<std::string>(
                    "it is named as a verification boundary")
              : whyNotOnlyRisingEdgeFlops(instance, parts.gate);
      if (whyNot.has_value()) {
        // Only the first few are spelled out; the count follows.
        if (++kept <= kMaxReportedKeptInstances) {
          SPDLOG_INFO(
              "Clock gate {} ({}) in {} keeps no model: {}",
              instance->getName().getString(),
              instance->getModel()->getName().getString(),
              design->getName().getString(), *whyNot);
        }
        continue;
      }
      replaceInstance(instance, parts.latch, parts.gate);
      ++replaced;
    }
  }
  if (kept > kMaxReportedKeptInstances) {
    SPDLOG_INFO("{} more clock gate instance(s) keep no model",
                kept - kMaxReportedKeptInstances);
  }
  return replaced;
}

}  // namespace KEPLER_FORMAL
