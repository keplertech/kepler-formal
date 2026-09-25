// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "PyLatchOptions.h"

#include <array>
#include <exception>
#include <string>

namespace KEPLER_FORMAL {
namespace {
constexpr std::array<std::string_view, 11> keys{
    "latch_support", "latch_input_changes", "latch_initial_inputs", "latch_initial_storage",
    "latch_workers", "latch_max_waves", "latch_max_states", "latch_max_transactions",
    "latch_max_symbolic_nodes", "latch_max_sat_conflicts", "latch_max_sat_decisions"};

bool number(PyObject* dictionary, const char* key, std::optional<size_t>& value) {
  auto* item = PyDict_GetItemString(dictionary, key);
  if (!item || item == Py_None) return true;
  if (!PyLong_Check(item) || PyBool_Check(item)) {
    PyErr_Format(PyExc_TypeError, "%s must be an integer or None", key);
    return false;
  }
  const auto parsed = PyLong_AsSize_t(item);
  if (PyErr_Occurred()) return false;
  value = parsed;
  return true;
}

bool initialBit(PyObject* dictionary, const char* key, std::optional<bool>& value) {
  std::optional<size_t> parsed;
  if (!number(dictionary, key, parsed)) return false;
  if (!parsed) return true;
  if (*parsed > 1) {
    PyErr_Format(PyExc_ValueError, "%s must explicitly be Boolean 0 or 1", key);
    return false;
  }
  value = *parsed != 0;
  return true;
}
}  // namespace

bool isBorrowedLatchOption(std::string_view key) {
  for (auto allowed : keys) if (key == allowed) return true;
  return false;
}

bool parseBorrowedLatchOptions(PyObject* dictionary, BorrowedLatchOptions& options,
                              bool isSec, bool hasLeafBoundaries) {
  if (auto* enabled = PyDict_GetItemString(dictionary, "latch_support")) {
    if (!PyBool_Check(enabled)) {
      PyErr_SetString(PyExc_TypeError, "latch_support must be a bool");
      return false;
    }
    options.enabled = enabled == Py_True;
  }
  auto* changes = PyDict_GetItemString(dictionary, "latch_input_changes");
  if (changes && changes != Py_None) {
    if (!PyUnicode_Check(changes)) {
      PyErr_SetString(PyExc_TypeError, "latch_input_changes must be a string or None");
      return false;
    }
    Py_ssize_t length = 0;
    const char* text = PyUnicode_AsUTF8AndSize(changes, &length);
    if (!text) return false;
    const std::string_view value(text, size_t(length));
    if (value == "any") options.inputChanges = LatchInputChanges::Any;
    else if (value == "single") options.inputChanges = LatchInputChanges::Single;
    else {
      PyErr_SetString(PyExc_ValueError, "latch_input_changes must be any or single");
      return false;
    }
  }
  if (!initialBit(dictionary, "latch_initial_inputs", options.initialInputs) ||
      !initialBit(dictionary, "latch_initial_storage", options.initialStorage) ||
      !number(dictionary, "latch_workers", options.workers) ||
      !number(dictionary, "latch_max_waves", options.maxWaves) ||
      !number(dictionary, "latch_max_states", options.maxStates) ||
      !number(dictionary, "latch_max_transactions", options.maxTransactions) ||
      !number(dictionary, "latch_max_symbolic_nodes", options.maxSymbolicNodes) ||
      !number(dictionary, "latch_max_sat_conflicts", options.maxSatConflicts) ||
      !number(dictionary, "latch_max_sat_decisions", options.maxSatDecisions)) return false;
  try {
    (void)options.validated(isSec, hasLeafBoundaries);
    return true;
  } catch (const std::exception& error) {
    PyErr_SetString(PyExc_ValueError, error.what());
    return false;
  }
}

}  // namespace KEPLER_FORMAL
