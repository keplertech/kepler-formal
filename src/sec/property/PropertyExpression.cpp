// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "property/PropertyExpression.h"

#include "property/BoolExprOps.h"

#include <cctype>
#include <optional>
#include <stdexcept>

namespace KEPLER_FORMAL::SEC {

namespace {

using namespace PropertyLogic;

const std::string kEmptyOp;

struct Range {
  size_t low = 0;
  size_t high = 0;
  bool unbounded = false;

  bool exact() const { return !unbounded && low == high; }
};

std::optional<size_t> parseCount(const PropertyTerm& term) {
  if (!term.isAtom || term.atom.empty()) {
    return std::nullopt;
  }
  size_t value = 0;
  for (const char character : term.atom) {
    if (!std::isdigit(static_cast<unsigned char>(character))) {
      return std::nullopt;
    }
    value = value * 10 + static_cast<size_t>(character - '0');
    if (value > 1'000'000) {
      return std::nullopt;
    }
  }
  return value;
}

// Reads the `low high` pair that precedes an operand; high may be `$`.
std::optional<Range> parseRange(const PropertyTerm& low, const PropertyTerm& high) {
  Range range;
  const auto lowValue = parseCount(low);
  if (!lowValue) {
    return std::nullopt;
  }
  range.low = *lowValue;
  if (high.isAtom && high.atom == "$") {
    range.unbounded = true;
    return range;
  }
  const auto highValue = parseCount(high);
  if (!highValue || *highValue < range.low) {
    return std::nullopt;
  }
  range.high = *highValue;
  return range;
}

// A monitor chain has one bit per delay step, so keep one property from
// silently allocating an unreasonable amount of state.
constexpr size_t kMaxDelay = 4096;

[[noreturn]] void malformed(const std::string& what) {
  throw std::invalid_argument("malformed property expression: " + what);
}

class Parser {
 public:
  explicit Parser(const std::string& text) : text_(text) {}

  PropertyTerm parseAll() {
    PropertyTerm term = parseTerm();
    skipSpace();
    if (position_ != text_.size()) {
      malformed("trailing text after the expression");
    }
    return term;
  }

 private:
  void skipSpace() {
    while (position_ < text_.size() &&
           std::isspace(static_cast<unsigned char>(text_[position_]))) {
      ++position_;
    }
  }

  PropertyTerm parseTerm() {
    skipSpace();
    if (position_ >= text_.size()) {
      malformed("unexpected end of text");
    }
    PropertyTerm term;
    if (text_[position_] == '(') {
      ++position_;
      for (;;) {
        skipSpace();
        if (position_ >= text_.size()) {
          malformed("missing closing parenthesis");
        }
        if (text_[position_] == ')') {
          ++position_;
          break;
        }
        term.items.push_back(parseTerm());
      }
      if (term.items.empty()) {
        malformed("empty list");
      }
      return term;
    }
    if (text_[position_] == ')') {
      malformed("unexpected closing parenthesis");
    }
    term.isAtom = true;
    while (position_ < text_.size() && text_[position_] != '(' &&
           text_[position_] != ')' &&
           !std::isspace(static_cast<unsigned char>(text_[position_]))) {
      term.atom += text_[position_++];
    }
    return term;
  }

  const std::string& text_;
  size_t position_ = 0;
};

bool isOp(const PropertyTerm& term, const char* name, size_t itemCount) {
  return !term.isAtom && term.items.size() == itemCount && term.op() == name;
}

}  // namespace

const std::string& PropertyTerm::op() const {
  if (isAtom || items.empty() || !items.front().isAtom) {
    return kEmptyOp;
  }
  return items.front().atom;
}

PropertyTerm parsePropertyTerm(const std::string& text) {
  return Parser(text).parseAll();
}

PropertyCompiler::PropertyCompiler(std::vector<BoolExpr*> leaves,
                                   VarAllocator allocateVar)
    : leaves_(std::move(leaves)), allocateVar_(std::move(allocateVar)) {}

BoolExpr* PropertyCompiler::leaf(const PropertyTerm& term) const {
  const auto index = term.items.size() == 2 ? parseCount(term.items[1])
                                            : std::nullopt;
  if (!index || *index >= leaves_.size()) {
    malformed("leaf index out of range");
  }
  return leaves_[*index];
}

// A Boolean term has a value in the current step and needs no monitor state.
bool PropertyCompiler::isBoolean(const PropertyTerm& term) const {
  if (term.isAtom) {
    return false;
  }
  const auto& op = term.op();
  if (op == "leaf") {
    return true;
  }
  if (op == "not" || op == "weak" || op == "strong") {
    return term.items.size() == 2 && isBoolean(term.items[1]);
  }
  if (op == "and" || op == "or" || op == "implies" || op == "property_implies" ||
      op == "iff") {
    return term.items.size() == 3 && isBoolean(term.items[1]) &&
           isBoolean(term.items[2]);
  }
  if (op == "if") {
    if (term.items.size() != 3 && term.items.size() != 4) {
      return false;
    }
    for (size_t i = 1; i < term.items.size(); ++i) {
      if (!isBoolean(term.items[i])) {
        return false;
      }
    }
    return true;
  }
  if (op == "concat") {
    // A single element without delay is the element itself.
    if (term.items.size() != 2 || term.items[1].isAtom ||
        term.items[1].items.size() != 3) {
      return false;
    }
    const auto range = parseRange(term.items[1].items[0], term.items[1].items[1]);
    return range && range->exact() && range->low == 0 &&
           isBoolean(term.items[1].items[2]);
  }
  return false;
}

BoolExpr* PropertyCompiler::boolean(const PropertyTerm& term) const {
  const auto& op = term.op();
  if (op == "leaf") {
    return leaf(term);
  }
  if (op == "not") {
    return notExpr(boolean(term.items[1]));
  }
  if (op == "weak" || op == "strong") {
    return boolean(term.items[1]);
  }
  if (op == "and") {
    return andExpr(boolean(term.items[1]), boolean(term.items[2]));
  }
  if (op == "or") {
    return orExpr(boolean(term.items[1]), boolean(term.items[2]));
  }
  if (op == "implies" || op == "property_implies") {
    return orExpr(notExpr(boolean(term.items[1])), boolean(term.items[2]));
  }
  if (op == "iff") {
    return notExpr(BoolExpr::Xor(boolean(term.items[1]), boolean(term.items[2])));
  }
  if (op == "if") {
    BoolExpr* condition = boolean(term.items[1]);
    BoolExpr* result = orExpr(notExpr(condition), boolean(term.items[2]));
    if (term.items.size() == 4) {
      result = andExpr(result, orExpr(condition, boolean(term.items[3])));
    }
    return result;
  }
  // concat with one zero-delay element, as accepted by isBoolean.
  return boolean(term.items[1].items[2]);
}

// The value `signal` had one step earlier. Monitor bits start at zero, so
// nothing is seen before the first step.
BoolExpr* PropertyCompiler::delayOne(BoolExpr* signal) {
  if (signal == constFalse()) {
    return constFalse();
  }
  const size_t var = allocateVar_();
  monitor_.push_back({var, signal});
  return BoolExpr::Var(var);
}

BoolExpr* PropertyCompiler::delayExact(BoolExpr* signal, size_t steps) {
  if (steps > kMaxDelay) {
    throw Unsupported{"a delay longer than " + std::to_string(kMaxDelay) +
                      " cycles is not supported"};
  }
  for (size_t step = 0; step < steps; ++step) {
    signal = delayOne(signal);
  }
  return signal;
}

// 1 when `signal` was 1 between `low` and `high` steps ago (inclusive), or at
// least `low` steps ago when the window is unbounded.
BoolExpr* PropertyCompiler::delayWindow(BoolExpr* signal,
                                        size_t low,
                                        size_t high,
                                        bool unbounded) {
  BoolExpr* current = delayExact(signal, low);
  if (unbounded) {
    // seen: `current` was 1 in an earlier step.
    const size_t var = allocateVar_();
    BoolExpr* seen = BoolExpr::Var(var);
    BoolExpr* result = orExpr(current, seen);
    monitor_.push_back({var, result});
    return result;
  }
  if (high - low > kMaxDelay) {
    throw Unsupported{"a delay window wider than " + std::to_string(kMaxDelay) +
                      " cycles is not supported"};
  }
  BoolExpr* result = current;
  for (size_t step = low; step < high; ++step) {
    current = delayOne(current);
    result = orExpr(result, current);
  }
  return result;
}

// 1 in a step where an attempt of `sequence`, begun in a step where `start`
// was 1, matches. Attempts begun in different steps share the monitor bits:
// only the existence of a match matters to the caller.
BoolExpr* PropertyCompiler::matchSequence(const PropertyTerm& sequence,
                                          BoolExpr* start) {
  if (isBoolean(sequence)) {
    return andExpr(start, boolean(sequence));
  }
  const auto& op = sequence.op();
  if (op == "weak" || op == "strong") {
    if (sequence.items.size() != 2) {
      malformed(op);
    }
    return matchSequence(sequence.items[1], start);
  }
  if (op == "concat") {
    BoolExpr* current = start;
    for (size_t i = 1; i < sequence.items.size(); ++i) {
      const auto& element = sequence.items[i];
      if (element.isAtom || element.items.size() != 3) {
        malformed("concat element");
      }
      const auto range = parseRange(element.items[0], element.items[1]);
      if (!range) {
        malformed("concat delay");
      }
      current = delayWindow(current, range->low, range->high, range->unbounded);
      current = matchSequence(element.items[2], current);
    }
    return current;
  }
  if (op == "or") {
    if (sequence.items.size() != 3) {
      malformed("or");
    }
    return orExpr(matchSequence(sequence.items[1], start),
                  matchSequence(sequence.items[2], start));
  }
  if (op == "repeat") {
    if (sequence.items.size() != 5 || !sequence.items[1].isAtom) {
      malformed("repeat");
    }
    const auto range = parseRange(sequence.items[2], sequence.items[3]);
    if (!range) {
      malformed("repeat range");
    }
    const auto& kind = sequence.items[1].atom;
    const auto& operand = sequence.items[4];
    if (!isBoolean(operand)) {
      throw Unsupported{"repetition of a sequence that is not a Boolean "
                        "expression is not supported yet"};
    }
    if (range->low == 0) {
      throw Unsupported{"a repetition that admits an empty match is not "
                        "supported yet"};
    }
    if (range->low > kMaxDelay || (!range->unbounded && range->high > kMaxDelay)) {
      throw Unsupported{"a repetition count above " + std::to_string(kMaxDelay) +
                        " is not supported"};
    }
    BoolExpr* value = boolean(operand);
    if (kind == "consecutive") {
      // stage k: the operand held in k consecutive steps since the start.
      BoolExpr* stage = andExpr(start, value);
      BoolExpr* result = range->low == 1 ? stage : constFalse();
      const size_t last = range->unbounded ? range->low : range->high;
      for (size_t count = 2; count <= last; ++count) {
        stage = andExpr(delayOne(stage), value);
        if (count >= range->low) {
          result = orExpr(result, stage);
        }
      }
      if (range->unbounded) {
        // Once `low` repetitions matched, each further step extends the match.
        const size_t var = allocateVar_();
        BoolExpr* extended = orExpr(stage, andExpr(BoolExpr::Var(var), value));
        monitor_.push_back({var, extended});
        return extended;
      }
      return result;
    }
    if (kind == "goto") {
      if (range->unbounded) {
        throw Unsupported{"an unbounded goto repetition is not supported yet"};
      }
      // Each stage waits for the next step where the operand holds.
      BoolExpr* stageStart = start;
      BoolExpr* result = constFalse();
      for (size_t count = 1; count <= range->high; ++count) {
        const size_t var = allocateVar_();
        BoolExpr* waiting = orExpr(stageStart, BoolExpr::Var(var));
        monitor_.push_back({var, andExpr(waiting, notExpr(value))});
        BoolExpr* matched = andExpr(waiting, value);
        if (count >= range->low) {
          result = orExpr(result, matched);
        }
        stageStart = delayOne(matched);
      }
      return result;
    }
    throw Unsupported{"nonconsecutive repetition is not supported yet"};
  }
  if (op == "throughout" || op == "intersect" || op == "within" ||
      op == "first_match" || op == "and") {
    throw Unsupported{"the sequence operator `" + op + "` is not supported yet"};
  }
  throw Unsupported{"`" + (op.empty() ? std::string("?") : op) +
                    "` is not supported as a sequence"};
}

// 1 in a step where an attempt of the sequence, begun in a step where `start`
// was 1, can no longer match. Supported when no alternative has to be
// tracked: exact delays, and one final bounded delay window.
BoolExpr* PropertyCompiler::failSequence(const PropertyTerm& sequence,
                                         BoolExpr* start) {
  if (isBoolean(sequence)) {
    return andExpr(start, notExpr(boolean(sequence)));
  }
  const auto& op = sequence.op();
  if (op == "weak") {
    if (sequence.items.size() != 2) {
      malformed(op);
    }
    return failSequence(sequence.items[1], start);
  }

  BoolExpr* fail = constFalse();
  BoolExpr* current = start;
  // Consumes `count` consecutive steps in which `value` must hold.
  const auto requireSteps = [&](BoolExpr* value, size_t count) {
    for (size_t step = 0; step < count; ++step) {
      if (step != 0) {
        current = delayOne(current);
      }
      fail = orExpr(fail, andExpr(current, notExpr(value)));
      current = andExpr(current, value);
    }
  };
  const auto requireElement = [&](const PropertyTerm& element) {
    if (isBoolean(element)) {
      requireSteps(boolean(element), 1);
      return;
    }
    if (element.op() == "repeat" && element.items.size() == 5 &&
        element.items[1].isAtom && element.items[1].atom == "consecutive" &&
        isBoolean(element.items[4])) {
      const auto range = parseRange(element.items[2], element.items[3]);
      if (!range) {
        malformed("repeat range");
      }
      if (range->low == 0 || range->low > kMaxDelay) {
        throw Unsupported{"this repetition count is not supported in a "
                          "consequent"};
      }
      if (!range->exact()) {
        throw Unsupported{"a repetition range in a consequent needs alternative "
                          "tracking, which is not supported yet"};
      }
      requireSteps(boolean(element.items[4]), range->low);
      return;
    }
    throw Unsupported{"this sequence form is not supported in a consequent yet"};
  };

  if (op == "repeat") {
    requireElement(sequence);
    return fail;
  }
  if (op != "concat") {
    throw Unsupported{"`" + (op.empty() ? std::string("?") : op) +
                      "` is not supported in a consequent yet"};
  }
  for (size_t i = 1; i < sequence.items.size(); ++i) {
    const auto& element = sequence.items[i];
    if (element.isAtom || element.items.size() != 3) {
      malformed("concat element");
    }
    const auto range = parseRange(element.items[0], element.items[1]);
    if (!range) {
      malformed("concat delay");
    }
    const bool last = i + 1 == sequence.items.size();
    if (range->exact()) {
      current = delayExact(current, range->low);
      requireElement(element.items[2]);
      continue;
    }
    if (!last || !isBoolean(element.items[2])) {
      throw Unsupported{"a delay range inside a consequent needs alternative "
                        "tracking, which is not supported yet"};
    }
    if (range->unbounded) {
      // A weak sequence that may wait forever cannot be violated here.
      return fail;
    }
    // The attempt fails when its last allowed step passes and the operand was
    // false in every step of the window.
    BoolExpr* value = boolean(element.items[2]);
    BoolExpr* deadline = delayExact(current, range->high);
    BoolExpr* falseSince = notExpr(value);
    for (size_t step = range->low; step < range->high; ++step) {
      falseSince = andExpr(notExpr(value), delayOne(falseSince));
    }
    return orExpr(fail, andExpr(deadline, falseSince));
  }
  return fail;
}

// 1 in a step where an attempt of `property`, begun in a step where `start`
// was 1, is found violated.
BoolExpr* PropertyCompiler::failProperty(const PropertyTerm& property,
                                         BoolExpr* start) {
  if (isBoolean(property)) {
    return andExpr(start, notExpr(boolean(property)));
  }
  const auto& op = property.op();
  const auto operandCount = property.items.size() - 1;
  const auto requireOperands = [&](size_t expected) {
    if (property.isAtom || operandCount != expected) {
      malformed(op);
    }
  };
  const auto requireBooleans = [&](const char* what) {
    for (size_t i = 1; i < property.items.size(); ++i) {
      if (!isBoolean(property.items[i])) {
        throw Unsupported{std::string(what) +
                          " is supported on Boolean operands only"};
      }
    }
  };

  if (op == "implies") {
    requireOperands(2);
    return failProperty(property.items[2],
                        matchSequence(property.items[1], start));
  }
  if (op == "implies_next") {
    requireOperands(2);
    return failProperty(property.items[2],
                        delayOne(matchSequence(property.items[1], start)));
  }
  if (op == "and") {
    requireOperands(2);
    return orExpr(failProperty(property.items[1], start),
                  failProperty(property.items[2], start));
  }
  if (op == "not") {
    requireOperands(1);
    // `not` of a sequence: the sequence must never match.
    const auto& operandOp = property.items[1].op();
    if (operandOp == "concat" || operandOp == "repeat") {
      return matchSequence(property.items[1], start);
    }
    throw Unsupported{"negation of a property that is not a sequence is not "
                      "supported"};
  }
  if (op == "if") {
    if (property.items.size() != 3 && property.items.size() != 4) {
      malformed(op);
    }
    if (!isBoolean(property.items[1])) {
      malformed("if condition");
    }
    BoolExpr* condition = boolean(property.items[1]);
    BoolExpr* fail =
        failProperty(property.items[2], andExpr(start, condition));
    if (property.items.size() == 4) {
      fail = orExpr(fail, failProperty(property.items[3],
                                       andExpr(start, notExpr(condition))));
    }
    return fail;
  }
  if (op == "nexttime") {
    if (operandCount == 1) {
      return failProperty(property.items[1], delayOne(start));
    }
    requireOperands(3);
    const auto range = parseRange(property.items[1], property.items[2]);
    if (!range || !range->exact()) {
      throw Unsupported{"nexttime with a range is not supported"};
    }
    return failProperty(property.items[3], delayExact(start, range->low));
  }
  if (op == "always") {
    if (operandCount != 1) {
      throw Unsupported{"always with a range is not supported yet"};
    }
    // The operand must hold from the start on: keep the attempt alive.
    const size_t var = allocateVar_();
    BoolExpr* active = orExpr(start, BoolExpr::Var(var));
    monitor_.push_back({var, active});
    return failProperty(property.items[1], active);
  }
  if (op == "until" || op == "until_with") {
    requireOperands(2);
    requireBooleans("until");
    BoolExpr* hold = boolean(property.items[1]);
    BoolExpr* release = boolean(property.items[2]);
    const size_t var = allocateVar_();
    BoolExpr* active = orExpr(start, BoolExpr::Var(var));
    // Keep waiting while the left side holds and the right side has not.
    monitor_.push_back({var, andExpr(active, andExpr(hold, notExpr(release)))});
    if (op == "until_with") {
      return andExpr(active, notExpr(hold));
    }
    return andExpr(active, andExpr(notExpr(hold), notExpr(release)));
  }
  if (op == "or" || op == "property_implies" || op == "iff") {
    requireBooleans(op == "or" ? "a disjunction of properties"
                               : "this property operator");
  }
  if (op == "s_until" || op == "s_until_with" || op == "s_eventually" ||
      op == "eventually" || op == "s_always" || op == "s_nexttime" ||
      op == "strong") {
    throw Unsupported{"`" + op + "` is a liveness obligation, which is not "
                      "supported yet"};
  }
  if (op == "concat" || op == "repeat" || op == "weak") {
    return failSequence(property, start);
  }
  throw Unsupported{"the property operator `" +
                    (op.empty() ? std::string("?") : op) +
                    "` is not supported yet"};
}

CompiledProperty PropertyCompiler::finish(BoolExpr* signal) {
  CompiledProperty result;
  result.signal = signal;
  result.monitor = std::move(monitor_);
  monitor_.clear();
  return result;
}

CompiledProperty PropertyCompiler::compileViolation(const PropertyTerm& property) {
  monitor_.clear();
  try {
    return finish(failProperty(property, constTrue()));
  } catch (const Unsupported& unsupported) {
    CompiledProperty result;
    result.unsupportedReason = unsupported.reason;
    return result;
  } catch (const std::invalid_argument& error) {
    CompiledProperty result;
    result.unsupportedReason = error.what();
    return result;
  }
}

CompiledProperty PropertyCompiler::compileCover(const PropertyTerm& property) {
  monitor_.clear();
  try {
    // A cover counts matches of a sequence. An implication also holds when
    // its left side never matches, which a cover must not report as a hit.
    const std::function<bool(const PropertyTerm&)> isSequence =
        [&](const PropertyTerm& term) {
          const auto& op = term.op();
          if (op == "leaf" || op == "concat" || op == "repeat") {
            return true;
          }
          if (op == "weak" || op == "strong") {
            return term.items.size() == 2 && isSequence(term.items[1]);
          }
          if (op == "and" || op == "or") {
            return term.items.size() == 3 && isSequence(term.items[1]) &&
                   isSequence(term.items[2]);
          }
          return false;
        };
    if (!isSequence(property)) {
      throw Unsupported{"cover is supported on sequences only, not on `" +
                        property.op() + "`"};
    }
    return finish(matchSequence(property, constTrue()));
  } catch (const Unsupported& unsupported) {
    CompiledProperty result;
    result.unsupportedReason = unsupported.reason;
    return result;
  } catch (const std::invalid_argument& error) {
    CompiledProperty result;
    result.unsupportedReason = error.what();
    return result;
  }
}

}  // namespace KEPLER_FORMAL::SEC
