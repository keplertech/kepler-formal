// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "KeplerNajaRuntime.h"

#include <memory>

#include "KeplerNajaProviderBuild.h"
#include "NLUniverse.h"
#include "PyNLUniverse.h"
#include "PySNLDesign.h"

namespace KEPLER_FORMAL {
namespace {

struct PyObjectDeleter {
  void operator()(PyObject* object) const { Py_XDECREF(object); }
};
using OwnedPyObject = std::unique_ptr<PyObject, PyObjectDeleter>;

bool validateProviderBuild() {
  // The GIL serializes this check. Each process checks the provider's bytes
  // once; subsequent calls still validate the live types and universe below.
  static bool validated = false;
  if (validated) {
    return true;
  }
  OwnedPyObject helper(PyImport_ImportModule("kepler_formal._provider_check"));
  if (helper == nullptr) {
    return false;
  }
  OwnedPyObject result(PyObject_CallMethod(
      helper.get(), "validate_provider", "sss", KEPLER_NAJA_PROVIDER_MANIFEST,
      KEPLER_NAJA_PROVIDER_VERSION, KEPLER_NAJA_PROVIDER_GIT_HASH));
  if (result == nullptr) {
    return false;
  }
  validated = true;
  return true;
}

bool requireIdentityString(PyObject* module, const char* method,
                           const char* expected, PyObject* errorType) {
  OwnedPyObject value(PyObject_CallMethod(module, method, nullptr));
  if (value == nullptr) {
    return false;
  }
  if (!PyUnicode_Check(value.get()) ||
      PyUnicode_CompareWithASCIIString(value.get(), expected) != 0) {
    PyErr_Format(errorType,
                 "Kepler Formal was built against NajaEDA with %s() = %s",
                 method, expected);
    return false;
  }
  return true;
}

}  // namespace

bool validateNajaRuntime(PyObject* errorType) {
  OwnedPyObject module(PyImport_ImportModule("najaeda.naja"));
  if (module == nullptr ||
      !requireIdentityString(module.get(), "getVersion",
                             KEPLER_NAJA_PROVIDER_VERSION, errorType) ||
      !requireIdentityString(module.get(), "getGitHash",
                             KEPLER_NAJA_PROVIDER_GIT_HASH, errorType) ||
      !validateProviderBuild()) {
    return false;
  }

  OwnedPyObject designType(PyObject_GetAttrString(module.get(), "SNLDesign"));
  OwnedPyObject universeType(PyObject_GetAttrString(module.get(), "NLUniverse"));
  if (designType == nullptr || universeType == nullptr) {
    return false;
  }
  // Compare live exported objects: the type objects Python sees must be the
  // ones this extension linked, or two Naja copies are loaded.
  if (designType.get() != reinterpret_cast<PyObject*>(&PYNAJA::PyTypeSNLDesign) ||
      universeType.get() != reinterpret_cast<PyObject*>(&PYNAJA::PyTypeNLUniverse)) {
    PyErr_SetString(
        errorType,
        "NajaEDA and Kepler Formal did not load the same native runtime");
    return false;
  }
  OwnedPyObject pythonUniverse(
      PyObject_CallMethod(universeType.get(), "get", nullptr));
  if (pythonUniverse == nullptr) {
    return false;
  }
  OwnedPyObject linkedUniverse(
      PYNAJA::PyNLUniverse_Link(naja::NL::NLUniverse::get()));
  if (linkedUniverse == nullptr) {
    return false;
  }
  if (pythonUniverse.get() != linkedUniverse.get()) {
    PyErr_SetString(errorType,
                    "NajaEDA and Kepler Formal disagree on the active universe");
    return false;
  }
  return true;
}

naja::NL::SNLDesign* unwrapNajaDesign(PyObject* object) {
  if (!validateNajaRuntime(PyExc_RuntimeError)) {
    return nullptr;
  }
  if (object == nullptr || !PyObject_TypeCheck(object, &PYNAJA::PyTypeSNLDesign)) {
    PyErr_SetString(PyExc_TypeError,
                    "Expected an SNLDesign from this NajaEDA runtime");
    return nullptr;
  }
  // Naja clears the wrapper's native pointer when the design, its library,
  // database or universe is destroyed.
  auto* design = reinterpret_cast<PYNAJA::PySNLDesign*>(object)->object_;
  if (design == nullptr) {
    PyErr_SetString(PyExc_ReferenceError, "The NajaEDA design was destroyed");
    return nullptr;
  }
  return design;
}

}  // namespace KEPLER_FORMAL
