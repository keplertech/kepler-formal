// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
// Exercises the actual native Python parser without requiring a Naja runtime.
#include "PyLatchOptions.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using namespace KEPLER_FORMAL;
size_t checks = 0;
void check(bool value, const std::string& detail) {
  ++checks;
  if (!value) throw std::runtime_error(detail);
}

BorrowedLatchOptions contract() {
  BorrowedLatchOptions options;
  options.inputChanges = LatchInputChanges::Single;
  options.initialInputs = false;
  options.initialStorage = false;
  return options;
}

void invalid(const BorrowedLatchOptions& options, const char* message,
             bool sec = true, bool boundaries = false) {
  try { (void)options.validated(sec, boundaries); }
  catch (const std::invalid_argument& error) {
    check(std::string(error.what()).find(message) != std::string::npos, error.what());
    return;
  }
  check(false, "invalid C++ options accepted");
}

void cppOptions() {
  for (bool sec : {false, true}) {
    for (bool boundaries : {false, true}) {
      auto defaults = BorrowedLatchOptions{}.validated(sec, boundaries);
      check(defaults.enabled && !defaults.initialInputs && !defaults.initialStorage,
            "default-on support invented an event contract");
      BorrowedLatchOptions disabled;
      disabled.enabled = false;
      check(!disabled.validated(sec, boundaries).enabled, "explicit false ignored");
    }
  }
  auto options = contract();
  auto result = options.validated(true, false);
  check(result.enabled && result.singleInputChange && result.initialInputs == false &&
        result.initialStorage == false, "explicit zeros did not survive");
  options.inputChanges = LatchInputChanges::Any;
  options.initialInputs = true;
  options.initialStorage = true;
  options.workers = 0;
  options.maxWaves = 17;
  options.maxStates = 23;
  options.maxTransactions = 31;
  options.maxSymbolicNodes = 1000;
  options.maxSatConflicts = 123;
  options.maxSatDecisions = 456;
  result = options.validated(true, false);
  check(!result.singleInputChange && result.initialInputs == true && result.initialStorage == true &&
        result.workers == 0 && result.limits.maxWaves == 17 && result.limits.maxBoundaryStates == 23 &&
        result.limits.maxTransactions == 31 && result.maxSymbolicNodes == 1000 &&
        result.maxSatConflicts == 123 && result.maxSatDecisions == 456, "explicit tuning lost");
  invalid(options, "only supported for SEC", false);
  invalid(options, "complete top interface", true, true);
  options.enabled = false;
  invalid(options, "requires latch_support");
  for (int field = 0; field < 3; ++field) {
    options = contract();
    if (field == 0) options.inputChanges.reset();
    if (field == 1) options.initialInputs.reset();
    if (field == 2) options.initialStorage.reset();
    invalid(options, "requires explicit");
  }
  options = contract();
  options.inputChanges = static_cast<LatchInputChanges>(999);
  invalid(options, "any or single");
  options = contract();
  options.workers = size_t(std::numeric_limits<int>::max()) + 1;
  invalid(options, "nonnegative int");
  for (int field = 0; field < 3; ++field) {
    options = contract();
    if (field == 0) options.maxWaves = 0;
    if (field == 1) options.maxStates = 0;
    if (field == 2) options.maxTransactions = 0;
    invalid(options, "positive integers");
  }
  options = {};
  options.workers = 0;
  invalid(options, "requires explicit");
  options.enabled = false;
  invalid(options, "requires latch_support");
  for (int field = 0; field < 3; ++field) {
    options = contract();
    auto& limit = field == 0 ? options.maxSymbolicNodes : field == 1 ? options.maxSatConflicts : options.maxSatDecisions;
    limit = 0;
    invalid(options, "positive integers");
    if (field) {
      limit = size_t(std::numeric_limits<unsigned>::max()) + 1;
      invalid(options, "unsigned int");
    }
    limit = 1;
    options.enabled = false;
    invalid(options, "requires latch_support");
  }
  check(isBorrowedLatchOption("latch_support") &&
        isBorrowedLatchOption("latch_max_transactions") &&
        !isBorrowedLatchOption(std::string_view("latch_support\0suffix", 20)) &&
        !isBorrowedLatchOption("latch_unknown"), "option allowlist incorrect");
}

void parsed(const std::string& expression, bool success, PyObject* exception = nullptr,
            bool sec = true, bool boundaries = false, bool expectedEnabled = true) {
  PyObject* scope = PyDict_New();
  PyDict_SetItemString(scope, "__builtins__", PyEval_GetBuiltins());
  PyObject* dictionary = PyRun_String(expression.c_str(), Py_eval_input, scope, scope);
  Py_DECREF(scope);
  check(dictionary && PyDict_Check(dictionary), "bad Python fixture: " + expression);
  BorrowedLatchOptions options;
  const bool accepted = parseBorrowedLatchOptions(dictionary, options, sec, boundaries);
  Py_DECREF(dictionary);
  check(accepted == success, "unexpected parser verdict: " + expression);
  if (!accepted) {
    check(PyErr_Occurred() && (!exception || PyErr_ExceptionMatches(exception)),
          "wrong parser exception: " + expression);
    PyErr_Clear();
  } else {
    check(!PyErr_Occurred(), "success retained Python exception");
    check(options.enabled == expectedEnabled, "parser changed the master gate");
    check(options.inputChanges.has_value() == options.initialInputs.has_value() &&
          options.initialInputs.has_value() == options.initialStorage.has_value(),
          "parser invented incomplete contract");
  }
}

void pythonOptions() {
  const std::string good = "{'latch_input_changes': 'single', "
      "'latch_initial_inputs': 0, 'latch_initial_storage': 0}";
  parsed("{}", true, nullptr, false, true);
  parsed("{}", true);
  parsed("{'latch_support': False}", true, nullptr, true, false, false);
  parsed("{'latch_support': True}", true);
  parsed(good, true);
  parsed(good + " | {'latch_support': True}", true);
  parsed(good + " | {'latch_support': False}", false, PyExc_ValueError);
  parsed(good + " | {'latch_input_changes': 'any', 'latch_initial_inputs': 1, 'latch_initial_storage': 1}", true);
  parsed(good + " | {'latch_workers': 0, 'latch_max_waves': 8, 'latch_max_states': 16, 'latch_max_transactions': 32}", true);
  parsed(good, false, PyExc_ValueError, false);
  parsed(good, false, PyExc_ValueError, true, true);
  for (const char* value : {"None", "0", "1", "'true'", "[]"})
    parsed(good + " | {'latch_support': " + value + "}", false, PyExc_TypeError);
  for (const char* key : {"latch_input_changes", "latch_initial_inputs", "latch_initial_storage"})
    parsed(good + " | {'" + key + "': None}", false, PyExc_ValueError);
  for (const char* value : {"''", "'ANY'", "'single\\x00tail'"})
    parsed(good + " | {'latch_input_changes': " + value + "}", false, PyExc_ValueError);
  for (const char* value : {"True", "1", "[]"})
    parsed(good + " | {'latch_input_changes': " + value + "}", false, PyExc_TypeError);
  for (const char* key : {"latch_initial_inputs", "latch_initial_storage"}) {
    for (const char* value : {"True", "False", "'0'", "0.0"})
      parsed(good + " | {'" + key + "': " + value + "}", false, PyExc_TypeError);
    parsed(good + " | {'" + key + "': 2}", false, PyExc_ValueError);
    parsed(good + " | {'" + key + "': -1}", false, PyExc_OverflowError);
  }
  for (const char* key : {"latch_workers", "latch_max_waves", "latch_max_states", "latch_max_transactions",
                        "latch_max_symbolic_nodes", "latch_max_sat_conflicts", "latch_max_sat_decisions"}) {
    for (const char* value : {"True", "'1'", "1.0"})
      parsed(good + " | {'" + key + "': " + value + "}", false, PyExc_TypeError);
    parsed(good + " | {'" + key + "': 2**128}", false, PyExc_OverflowError);
    parsed(good + " | {'" + key + "': -1}", false, PyExc_OverflowError);
    parsed("{'" + std::string(key) + "': 1}", false, PyExc_ValueError);
  }
  parsed(good + " | {'latch_workers': 2**31}", false, PyExc_ValueError);
  for (const char* key : {"latch_max_waves", "latch_max_states", "latch_max_transactions",
                        "latch_max_symbolic_nodes", "latch_max_sat_conflicts", "latch_max_sat_decisions"})
    parsed(good + " | {'" + key + "': 0}", false, PyExc_ValueError);
  for (const char* key : {"latch_max_sat_conflicts", "latch_max_sat_decisions"})
    parsed(good + " | {'" + key + "': 2**32}", false, PyExc_ValueError);
  parsed("{'latch_input_changes': 'any'}", false, PyExc_ValueError);
  parsed("{'latch_initial_inputs': 0}", false, PyExc_ValueError);
  parsed("{'latch_initial_storage': 0}", false, PyExc_ValueError);
}
}  // namespace

int main() {
  Py_Initialize();
  try {
    cppOptions();
    pythonOptions();
    Py_Finalize();
    std::cout << "Borrowed latch options: " << checks << " checks passed\n";
    return 0;
  } catch (const std::exception& error) {
    if (PyErr_Occurred()) PyErr_Print();
    Py_Finalize();
    std::cerr << error.what() << '\n';
    return 1;
  }
}
