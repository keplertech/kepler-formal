// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "KeplerNajaRuntime.h"

#include <limits>
#include <memory>

#include "KeplerPublishedNajaBuild.h"
#include "NLID.h"
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
  OwnedPyObject helper(PyImport_ImportModule("kepler_formal._published_runtime"));
  if (helper == nullptr) {
    return false;
  }
  OwnedPyObject result(PyObject_CallMethod(
      helper.get(), "validate_provider", "s", KEPLER_NAJA_PROVIDER_MANIFEST));
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
                 "Kepler Formal's published NajaEDA adapter requires %s() = %s",
                 method, expected);
    return false;
  }
  return true;
}

template <typename ID>
bool readID(PyObject* identifier, const char* method, ID& result) {
  OwnedPyObject value(PyObject_CallMethod(identifier, method, nullptr));
  if (value == nullptr) {
    return false;
  }
  if (!PyLong_Check(value.get()) || PyBool_Check(value.get())) {
    PyErr_SetString(PyExc_TypeError, "NajaEDA design IDs must be integers");
    return false;
  }
  const auto number = PyLong_AsUnsignedLongLong(value.get());
  if (PyErr_Occurred()) {
    return false;
  }
  if (number > std::numeric_limits<ID>::max()) {
    PyErr_SetString(PyExc_OverflowError, "NajaEDA design ID is out of range");
    return false;
  }
  result = static_cast<ID>(number);
  return true;
}

}  // namespace

bool validateNajaRuntime(PyObject* errorType) {
  OwnedPyObject module(PyImport_ImportModule("najaeda.naja"));
  if (module == nullptr ||
      !requireIdentityString(module.get(), "getVersion", "0.7.24", errorType) ||
      !requireIdentityString(module.get(), "getGitHash", "2263958", errorType) ||
      !validateProviderBuild()) {
    return false;
  }

  // This adapter targets the published pre-SDK runtime exclusively. A capsule
  // belongs to another contract and must not silently override these checks.
  OwnedPyObject capsule(PyObject_GetAttrString(module.get(), "_C_API"));
  if (capsule != nullptr) {
    PyErr_SetString(errorType,
                    "Incompatible NajaEDA native runtime API for the published adapter");
    return false;
  }
  if (!PyErr_ExceptionMatches(PyExc_AttributeError)) {
    return false;
  }
  PyErr_Clear();

  OwnedPyObject designType(PyObject_GetAttrString(module.get(), "SNLDesign"));
  OwnedPyObject universeType(PyObject_GetAttrString(module.get(), "NLUniverse"));
  if (designType == nullptr || universeType == nullptr) {
    return false;
  }
  // Compare live exported objects, not type names or Python wrapper layouts.
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

  // Ask the retained owner on every call. Its getter rejects destruction even
  // when a replacement design has already reused the same numerical IDs.
  OwnedPyObject identifier(PyObject_CallMethod(object, "getNLID", nullptr));
  if (identifier == nullptr) {
    if (PyErr_ExceptionMatches(PyExc_RuntimeError)) {
      PyErr_Clear();
      PyErr_SetString(PyExc_ReferenceError, "The NajaEDA design was destroyed");
    }
    return nullptr;
  }
  naja::NL::NLID::DBID dbID;
  naja::NL::NLID::LibraryID libraryID;
  naja::NL::NLID::DesignID designID;
  if (!readID(identifier.get(), "getDBID", dbID) ||
      !readID(identifier.get(), "getLibraryID", libraryID) ||
      !readID(identifier.get(), "getDesignID", designID)) {
    return nullptr;
  }
  const auto* universe = naja::NL::NLUniverse::get();
  auto* design = universe == nullptr
                     ? nullptr
                     : universe->getSNLDesign(
                           naja::NL::NLID::DesignReference(dbID, libraryID, designID));
  if (design == nullptr) {
    PyErr_SetString(PyExc_ReferenceError,
                    "The NajaEDA design is not in the active universe");
    return nullptr;
  }
  OwnedPyObject linkedDesign(PYNAJA::PySNLDesign_Link(design));
  if (linkedDesign == nullptr) {
    return nullptr;
  }
  if (linkedDesign.get() != object) {
    PyErr_SetString(PyExc_ReferenceError,
                    "The NajaEDA design no longer matches its native object");
    return nullptr;
  }
  return design;
}

}  // namespace KEPLER_FORMAL
