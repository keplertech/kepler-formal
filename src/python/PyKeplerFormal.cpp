// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "KeplerBorrowedDesigns.h"
#include "KeplerNajaRuntime.h"
#include "SNLDesign.h"

#ifndef KEPLER_FORMAL_VERSION
#define KEPLER_FORMAL_VERSION "0+unknown"
#endif

#ifndef KEPLER_FORMAL_GIT_HASH
#define KEPLER_FORMAL_GIT_HASH "unknown"
#endif

namespace {

struct PyObjectDeleter {
  void operator()(PyObject *object) const { Py_XDECREF(object); }
};

using OwnedPyObject = std::unique_ptr<PyObject, PyObjectDeleter>;

typedef struct {
  PyObject_HEAD
  PyObject *designOwner;
  PyObject *sourceOwner;
  void *design;
} PyNativeDesign;

PyTypeObject *nativeDesignType = nullptr;

bool ensureGilEnabled() {
#ifdef Py_GIL_DISABLED
  // Importing this GIL-requiring module normally enables the GIL, but an
  // explicit PYTHON_GIL=0 / -X gil=0 overrides that protection.
  OwnedPyObject sys(PyImport_ImportModule("sys"));
  if (sys == nullptr) {
    return false;
  }
  OwnedPyObject isGilEnabled(PyObject_GetAttrString(sys.get(), "_is_gil_enabled"));
  if (isGilEnabled == nullptr) {
    return false;
  }
  OwnedPyObject gilEnabled(PyObject_CallNoArgs(isGilEnabled.get()));
  if (gilEnabled == nullptr) {
    return false;
  }
  const int hasGil = PyObject_IsTrue(gilEnabled.get());
  if (hasGil < 0) {
    return false;
  }
  if (!hasGil) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Kepler Formal requires Python's GIL; use PYTHON_GIL=1 "
                    "or -X gil=1 for verification");
    return false;
  }
#endif
  return true;
}

void nativeDesignDealloc(PyNativeDesign *self) {
  PyObject_GC_UnTrack(self);
  Py_CLEAR(self->designOwner);
  Py_CLEAR(self->sourceOwner);
  PyTypeObject *type = Py_TYPE(self);
  type->tp_free(reinterpret_cast<PyObject *>(self));
  Py_DECREF(type);
}

int nativeDesignTraverse(PyNativeDesign *self, visitproc visit, void *arg) {
  Py_VISIT(self->designOwner);
  Py_VISIT(self->sourceOwner);
  return 0;
}

int nativeDesignClear(PyNativeDesign *self) {
  Py_CLEAR(self->designOwner);
  Py_CLEAR(self->sourceOwner);
  self->design = nullptr;
  return 0;
}

PyObject *nativeDesignNew(PyTypeObject *, PyObject *, PyObject *) {
  PyErr_SetString(PyExc_TypeError,
                  "NativeDesign objects are created by from_najaeda()");
  return nullptr;
}

PyObject *nativeDesignRepr(PyNativeDesign *) {
  return PyUnicode_FromString("<kepler_formal.NativeDesign: borrowed NajaEDA design>");
}

PyObject *nativeDesignSource(PyNativeDesign *self, void *) {
  if (self->sourceOwner == nullptr) {
    PyErr_SetString(PyExc_ReferenceError, "NativeDesign is no longer valid");
    return nullptr;
  }
  return Py_NewRef(self->sourceOwner);
}

PyObject *nativeDesignNajaedaDesign(PyNativeDesign *self, void *) {
  if (self->designOwner == nullptr) {
    PyErr_SetString(PyExc_ReferenceError, "NativeDesign is no longer valid");
    return nullptr;
  }
  return Py_NewRef(self->designOwner);
}

PyGetSetDef nativeDesignGetSet[] = {
    {const_cast<char *>("source"),
     reinterpret_cast<getter>(nativeDesignSource), nullptr,
     const_cast<char *>("Original object captured by from_najaeda()."), nullptr},
    {const_cast<char *>("najaeda_design"),
     reinterpret_cast<getter>(nativeDesignNajaedaDesign), nullptr,
     const_cast<char *>("Resolved raw najaeda.naja.SNLDesign."), nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr},
};

PyType_Slot nativeDesignSlots[] = {
    {Py_tp_dealloc, reinterpret_cast<void *>(nativeDesignDealloc)},
    {Py_tp_traverse, reinterpret_cast<void *>(nativeDesignTraverse)},
    {Py_tp_clear, reinterpret_cast<void *>(nativeDesignClear)},
    {Py_tp_new, reinterpret_cast<void *>(nativeDesignNew)},
    {Py_tp_repr, reinterpret_cast<void *>(nativeDesignRepr)},
    {Py_tp_getset, reinterpret_cast<void *>(nativeDesignGetSet)},
    {Py_tp_doc,
     const_cast<char *>(
         "Stable borrowed-design handle retaining its NajaEDA Python wrappers.")},
    {0, nullptr},
};

PyType_Spec nativeDesignSpec = {
    "kepler_formal.NativeDesign",
    sizeof(PyNativeDesign),
    0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    nativeDesignSlots,
};

int setDictItem(PyObject *dictionary, const char *key, PyObject *value) {
  if (value == nullptr) {
    return -1;
  }
  const int status = PyDict_SetItemString(dictionary, key, value);
  Py_DECREF(value);
  return status;
}

PyObject *optionalString(const std::string &value) {
  if (value.empty()) {
    Py_RETURN_NONE;
  }
  return PyUnicode_DecodeUTF8(value.data(), value.size(), "surrogateescape");
}

PyObject *stringList(const std::vector<std::string> &values) {
  PyObject *list = PyList_New(static_cast<Py_ssize_t>(values.size()));
  if (list == nullptr) {
    return nullptr;
  }
  for (size_t i = 0; i < values.size(); ++i) {
    PyObject *value = PyUnicode_DecodeUTF8(values[i].data(), values[i].size(),
                                           "surrogateescape");
    if (value == nullptr) {
      Py_DECREF(list);
      return nullptr;
    }
    PyList_SET_ITEM(list, static_cast<Py_ssize_t>(i), value);
  }
  return list;
}

PyObject *resultToDictionary(const KEPLER_FORMAL::RunResult &result) {
  PyObject *dictionary = PyDict_New();
  if (dictionary == nullptr) {
    return nullptr;
  }

  const auto set = [&](const char *key, PyObject *value) {
    if (setDictItem(dictionary, key, value) < 0) {
      Py_DECREF(dictionary);
      return false;
    }
    return true;
  };

  if (!set("status",
           PyUnicode_FromString(KEPLER_FORMAL::runStatusName(result.status))) ||
      !set("exit_code", PyLong_FromLong(result.exitCode)) ||
      !set("input_format", optionalString(result.inputFormat)) ||
      !set("verification", optionalString(result.verification)) ||
      !set("log_file", optionalString(result.logFile)) ||
      !set("bound", PyLong_FromSize_t(result.bound)) ||
      !set("reason", optionalString(result.reason)) ||
      !set("covered_outputs", PyLong_FromSize_t(result.coveredOutputs)) ||
      !set("total_outputs", PyLong_FromSize_t(result.totalOutputs)) ||
      !set("proven_outputs", PyLong_FromSize_t(result.provenOutputs)) ||
      !set("unproven_outputs", stringList(result.unprovenOutputs)) ||
      !set("skipped_observed_outputs",
           stringList(result.skippedObservedOutputs))) {
    return nullptr;
  }
  return dictionary;
}


bool dictionaryString(PyObject *dictionary, const char *key,
                      std::string &value, bool allowNone = false) {
  PyObject *item = PyDict_GetItemString(dictionary, key);
  if (item == nullptr || (allowNone && item == Py_None)) {
    return true;
  }
  if (!PyUnicode_Check(item)) {
    PyErr_Format(PyExc_TypeError, "%s must be a string%s", key,
                 allowNone ? " or None" : "");
    return false;
  }
  Py_ssize_t size = 0;
  const char *data = PyUnicode_AsUTF8AndSize(item, &size);
  if (data == nullptr) {
    return false;
  }
  value.assign(data, static_cast<size_t>(size));
  if (value.find('\0') != std::string::npos) {
    PyErr_Format(PyExc_ValueError, "%s cannot contain NUL bytes", key);
    return false;
  }
  return true;
}

bool dictionaryBoolean(PyObject *dictionary, const char *key, bool &value) {
  PyObject *item = PyDict_GetItemString(dictionary, key);
  if (item == nullptr) {
    return true;
  }
  if (!PyBool_Check(item)) {
    PyErr_Format(PyExc_TypeError, "%s must be a bool", key);
    return false;
  }
  value = item == Py_True;
  return true;
}

bool dictionarySize(PyObject *dictionary, const char *key, size_t &value) {
  PyObject *item = PyDict_GetItemString(dictionary, key);
  if (item == nullptr || item == Py_None) {
    return true;
  }
  if (!PyLong_Check(item) || PyBool_Check(item)) {
    PyErr_Format(PyExc_TypeError, "%s must be a non-negative integer", key);
    return false;
  }
  const size_t parsed = PyLong_AsSize_t(item);
  if (parsed == static_cast<size_t>(-1) && PyErr_Occurred()) {
    return false;
  }
  value = parsed;
  return true;
}

bool dictionaryBoundaryPairs(PyObject *dictionary, const char *key,
                             KEPLER_FORMAL::BoundaryPairs &value) {
  PyObject *item = PyDict_GetItemString(dictionary, key);
  if (item == nullptr) {
    return true;
  }
  if (!PyList_Check(item) && !PyTuple_Check(item)) {
    PyErr_Format(PyExc_TypeError,
                 "%s must be a list or tuple of path pairs", key);
    return false;
  }
  OwnedPyObject pairs(
      PySequence_Fast(item, "set_as_boundary must be a list or tuple"));
  if (pairs == nullptr) {
    return false;
  }
  KEPLER_FORMAL::BoundaryPairs parsed;
  const Py_ssize_t pairCount = PySequence_Fast_GET_SIZE(pairs.get());
  parsed.reserve(static_cast<size_t>(pairCount));
  for (Py_ssize_t index = 0; index < pairCount; ++index) {
    PyObject *pair = PySequence_Fast_GET_ITEM(pairs.get(), index);
    if (!PyList_Check(pair) && !PyTuple_Check(pair)) {
      PyErr_Format(PyExc_TypeError,
                   "%s[%zd] must be a list or tuple of two paths", key,
                   index);
      return false;
    }
    OwnedPyObject pathPair(
        PySequence_Fast(pair, "boundary path pair must be a list or tuple"));
    if (pathPair == nullptr) {
      return false;
    }
    if (PySequence_Fast_GET_SIZE(pathPair.get()) != 2) {
      PyErr_Format(PyExc_ValueError,
                   "%s[%zd] must contain exactly two paths", key, index);
      return false;
    }
    std::string paths[2];
    for (Py_ssize_t side = 0; side < 2; ++side) {
      PyObject *path = PySequence_Fast_GET_ITEM(pathPair.get(), side);
      if (!PyUnicode_Check(path)) {
        PyErr_Format(PyExc_TypeError, "%s[%zd][%zd] must be a string", key,
                     index, side);
        return false;
      }
      Py_ssize_t size = 0;
      const char *data = PyUnicode_AsUTF8AndSize(path, &size);
      if (data == nullptr) {
        return false;
      }
      paths[side].assign(data, static_cast<size_t>(size));
      if (paths[side].empty()) {
        PyErr_Format(PyExc_ValueError, "%s[%zd][%zd] must not be empty", key,
                     index, side);
        return false;
      }
      if (paths[side].find('\0') != std::string::npos) {
        PyErr_Format(PyExc_ValueError,
                     "%s[%zd][%zd] cannot contain NUL bytes", key, index,
                     side);
        return false;
      }
    }
    parsed.emplace_back(std::move(paths[0]), std::move(paths[1]));
  }
  value = std::move(parsed);
  return true;
}

bool parseBorrowedOptions(PyObject *object,
                          KEPLER_FORMAL::BorrowedDesignOptions &options) {
  if (!PyDict_Check(object)) {
    PyErr_SetString(PyExc_TypeError,
                    "verify_designs() native options must be a dict");
    return false;
  }
  static const std::unordered_set<std::string_view> allowedKeys = {
      "mode",          "solver",          "max_k",
      "sec_engine",    "sec_encoding",    "allow_boundary_mismatch",
      "set_as_boundary", "report_skipped_outputs", "log_file", "log_level"};
  Py_ssize_t position = 0;
  PyObject *key = nullptr;
  PyObject *value = nullptr;
  while (PyDict_Next(object, &position, &key, &value)) {
    if (!PyUnicode_Check(key)) {
      PyErr_SetString(PyExc_TypeError,
                      "verify_designs() option names must be strings");
      return false;
    }
    const char *name = PyUnicode_AsUTF8(key);
    if (name == nullptr) {
      return false;
    }
    if (!allowedKeys.contains(name)) {
      PyErr_Format(PyExc_TypeError,
                   "unknown verify_designs() native option: %s", name);
      return false;
    }
  }

  std::string mode = "lec";
  std::string solver = "kissat";
  std::string secEngine = "pdr";
  std::string secEncoding = "dual_rail_steady";
  if (!dictionaryString(object, "mode", mode) ||
      !dictionaryString(object, "solver", solver) ||
      !dictionaryString(object, "sec_engine", secEngine) ||
      !dictionaryString(object, "sec_encoding", secEncoding) ||
      !dictionaryString(object, "log_file", options.logFile, true) ||
      !dictionaryString(object, "log_level", options.logLevel, true) ||
      !dictionarySize(object, "max_k", options.maxK) ||
      !dictionaryBoundaryPairs(object, "set_as_boundary",
                               options.setAsBoundary) ||
      !dictionaryBoolean(object, "allow_boundary_mismatch",
                         options.allowBoundaryMismatch) ||
      !dictionaryBoolean(object, "report_skipped_outputs",
                         options.reportSkippedOutputs)) {
    return false;
  }

  if (mode == "lec") {
    options.mode = KEPLER_FORMAL::BorrowedVerificationMode::LEC;
  } else if (mode == "sec") {
    options.mode = KEPLER_FORMAL::BorrowedVerificationMode::SEC;
  } else {
    PyErr_SetString(PyExc_ValueError, "mode must be 'lec' or 'sec'");
    return false;
  }
  if (solver == "kissat") {
    options.solver = KEPLER_FORMAL::Config::SolverType::KISSAT;
  } else if (solver == "cadical") {
    options.solver = KEPLER_FORMAL::Config::SolverType::CADICAL;
  } else if (solver == "glucose") {
    options.solver = KEPLER_FORMAL::Config::SolverType::GLUCOSE;
  } else {
    PyErr_SetString(PyExc_ValueError,
                    "solver must be 'kissat', 'cadical', or 'glucose'");
    return false;
  }
  if (secEngine == "pdr") {
    options.secEngine = KEPLER_FORMAL::SEC::SecEngine::Pdr;
  } else if (secEngine == "k_induction") {
    options.secEngine = KEPLER_FORMAL::SEC::SecEngine::KInduction;
  } else if (secEngine == "imc") {
    options.secEngine = KEPLER_FORMAL::SEC::SecEngine::Imc;
  } else {
    PyErr_SetString(
        PyExc_ValueError,
        "sec_engine must be 'pdr', 'k_induction', or 'imc'");
    return false;
  }
  if (secEncoding == "dual_rail_steady") {
    options.secEncoding = KEPLER_FORMAL::SEC::SecEncoding::DualRailSteady;
  } else if (secEncoding == "binary") {
    options.secEncoding = KEPLER_FORMAL::SEC::SecEncoding::Binary;
  } else {
    PyErr_SetString(PyExc_ValueError,
                    "sec_encoding must be 'dual_rail_steady' or 'binary'");
    return false;
  }
  if (!options.logLevel.empty() && options.logLevel != "debug" &&
      options.logLevel != "info") {
    PyErr_SetString(PyExc_ValueError,
                    "log_level must be 'debug', 'info', or None");
    return false;
  }
  return true;
}

PyObject *fromNajaeda(PyObject *, PyObject *args) {
  if (!ensureGilEnabled()) {
    return nullptr;
  }
  PyObject *designOwner = nullptr;
  PyObject *sourceOwner = nullptr;
  if (!PyArg_ParseTuple(args, "O|O:from_najaeda", &designOwner,
                        &sourceOwner)) {
    return nullptr;
  }
  if (sourceOwner == nullptr) {
    sourceOwner = designOwner;
  }
  if (!KEPLER_FORMAL::validateNajaRuntime(PyExc_RuntimeError)) {
    return nullptr;
  }
  void *design = KEPLER_FORMAL::unwrapNajaDesign(designOwner);
  if (design == nullptr) {
    return nullptr;
  }
  auto *handle = reinterpret_cast<PyNativeDesign *>(
      nativeDesignType->tp_alloc(nativeDesignType, 0));
  if (handle == nullptr) {
    return nullptr;
  }
  handle->designOwner = Py_NewRef(designOwner);
  handle->sourceOwner = Py_NewRef(sourceOwner);
  handle->design = design;
  return reinterpret_cast<PyObject *>(handle);
}

naja::NL::SNLDesign *unwrapNativeDesign(PyObject *object, const char *label) {
  if (nativeDesignType == nullptr ||
      !PyObject_TypeCheck(object, nativeDesignType)) {
    PyErr_Format(PyExc_TypeError, "%s must be a NativeDesign", label);
    return nullptr;
  }
  auto *handle = reinterpret_cast<PyNativeDesign *>(object);
  void *current = KEPLER_FORMAL::unwrapNajaDesign(handle->designOwner);
  if (current == nullptr) {
    return nullptr;
  }
  if (current != handle->design) {
    PyErr_Format(PyExc_ReferenceError,
                 "%s no longer refers to the captured NajaEDA design", label);
    return nullptr;
  }
  return static_cast<naja::NL::SNLDesign *>(current);
}

PyObject *verifyDesigns(PyObject *, PyObject *args) {
  if (!ensureGilEnabled()) {
    return nullptr;
  }
  PyObject *firstObject = nullptr;
  PyObject *secondObject = nullptr;
  PyObject *optionObject = nullptr;
  if (!PyArg_ParseTuple(args, "OOO:verify_designs", &firstObject,
                        &secondObject, &optionObject)) {
    return nullptr;
  }
  KEPLER_FORMAL::BorrowedDesignOptions options;
  if (!parseBorrowedOptions(optionObject, options)) {
    return nullptr;
  }
  if (!KEPLER_FORMAL::validateNajaRuntime(PyExc_RuntimeError)) {
    return nullptr;
  }
  auto *first = unwrapNativeDesign(firstObject, "design1");
  if (first == nullptr) {
    return nullptr;
  }
  auto *second = unwrapNativeDesign(secondObject, "design2");
  if (second == nullptr) {
    return nullptr;
  }
  try {
    KEPLER_FORMAL::RunResult result;
    KEPLER_FORMAL::verifyBorrowedDesigns(first, second, options, result);
    return resultToDictionary(result);
  } catch (const std::exception &error) {
    PyErr_SetString(PyExc_RuntimeError, error.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError,
                    "unknown native borrowed-design verification failure");
    return nullptr;
  }
}


PyObject *version(PyObject *, PyObject *) {
  return PyUnicode_FromString(KEPLER_FORMAL_VERSION);
}

PyObject *gitHash(PyObject *, PyObject *) {
  return PyUnicode_FromString(KEPLER_FORMAL_GIT_HASH);
}

PyMethodDef methods[] = {
    {"from_najaeda", fromNajaeda, METH_VARARGS,
     "Capture a live NajaEDA SNLDesign without copying its native netlist."},
    {"verify_designs", verifyDesigns, METH_VARARGS,
     "Verify two captured NajaEDA designs without taking ownership."},
    {"get_version", version, METH_NOARGS, "Return the Kepler Formal version."},
    {"get_git_hash", gitHash, METH_NOARGS, "Return the build git hash."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_native",
    "Native Kepler Formal verification API.",
    -1,
    methods,
};

} // namespace

PyMODINIT_FUNC PyInit__native() {
  if (!KEPLER_FORMAL::validateNajaRuntime(PyExc_ImportError)) {
    return nullptr;
  }

  OwnedPyObject nativeDesign(PyType_FromSpec(&nativeDesignSpec));
  if (nativeDesign == nullptr) {
    return nullptr;
  }
  nativeDesignType = reinterpret_cast<PyTypeObject *>(nativeDesign.get());

  OwnedPyObject result(PyModule_Create(&module));
  if (result == nullptr) {
    nativeDesignType = nullptr;
    return nullptr;
  }
  if (PyModule_AddObjectRef(result.get(), "NativeDesign", nativeDesign.get()) < 0) {
    nativeDesignType = nullptr;
    return nullptr;
  }
#ifdef KEPLER_USE_PUBLISHED_NAJAEDA
  constexpr const char* providerMode = "published";
#else
  constexpr const char* providerMode = "sdk";
#endif
  if (PyModule_AddStringConstant(result.get(), "_provider_mode", providerMode) < 0) {
    nativeDesignType = nullptr;
    return nullptr;
  }
  return result.release();
}
