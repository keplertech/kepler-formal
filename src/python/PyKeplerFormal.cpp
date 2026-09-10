// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <climits>
#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "KeplerBorrowedDesigns.h"
#include "KeplerFormalDriver.h"
#include "NajaPythonRuntimeAPI.h"
#include "NajaRuntimeBuild.h"
#include "NLUniverse.h"
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

bool getRuntimeAPI(const NajaPythonRuntimeAPI *&api,
                   PyObject *errorType = PyExc_RuntimeError) {
  api = NajaPythonRuntime_Import();
  if (api == nullptr) {
    return false;
  }
  if (api->build_id == nullptr || api->runtime_identity == nullptr ||
      api->get_universe == nullptr || api->unwrap_design == nullptr) {
    PyErr_SetString(errorType, "Incomplete NajaEDA native runtime API");
    return false;
  }
  if (std::strcmp(api->build_id, NAJA_RUNTIME_BUILD_ID) != 0) {
    PyErr_Format(errorType,
                 "NajaEDA native build mismatch: provider %.80s, Kepler %.80s",
                 api->build_id, NAJA_RUNTIME_BUILD_ID);
    return false;
  }
  const auto localIdentity = naja::NL::NLUniverse::getRuntimeIdentity();
  if (api->runtime_identity != localIdentity) {
    PyErr_SetString(
        errorType,
        "NajaEDA and Kepler Formal did not load the same native runtime");
    return false;
  }
  if (api->get_universe() != naja::NL::NLUniverse::get()) {
    PyErr_SetString(errorType,
                    "NajaEDA and Kepler Formal disagree on the active universe");
    return false;
  }
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

bool appendArgument(PyObject *value, std::vector<std::string> &arguments) {
  OwnedPyObject path(PyOS_FSPath(value));
  if (path == nullptr) {
    return false;
  }

  PyObject *bytes = nullptr;
  if (!PyUnicode_FSConverter(path.get(), &bytes)) {
    return false;
  }
  OwnedPyObject ownedBytes(bytes);

  char *data = nullptr;
  Py_ssize_t size = 0;
  if (PyBytes_AsStringAndSize(ownedBytes.get(), &data, &size) < 0) {
    return false;
  }
  std::string argument(data, static_cast<size_t>(size));
  if (argument.find('\0') != std::string::npos) {
    PyErr_SetString(PyExc_ValueError,
                    "Kepler Formal arguments cannot contain NUL bytes");
    return false;
  }
  arguments.push_back(std::move(argument));
  return true;
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
      "report_skipped_outputs", "log_file", "log_level"};
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
  const NajaPythonRuntimeAPI *api = nullptr;
  if (!getRuntimeAPI(api)) {
    return nullptr;
  }
  void *design = api->unwrap_design(designOwner);
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

naja::NL::SNLDesign *unwrapNativeDesign(PyObject *object, const char *label,
                                       const NajaPythonRuntimeAPI *api) {
  if (nativeDesignType == nullptr ||
      !PyObject_TypeCheck(object, nativeDesignType)) {
    PyErr_Format(PyExc_TypeError, "%s must be a NativeDesign", label);
    return nullptr;
  }
  auto *handle = reinterpret_cast<PyNativeDesign *>(object);
  void *current = api->unwrap_design(handle->designOwner);
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
  const NajaPythonRuntimeAPI *api = nullptr;
  if (!getRuntimeAPI(api)) {
    return nullptr;
  }
  auto *first = unwrapNativeDesign(firstObject, "design1", api);
  if (first == nullptr) {
    return nullptr;
  }
  auto *second = unwrapNativeDesign(secondObject, "design2", api);
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

PyObject *run(PyObject *, PyObject *args) {
  if (!ensureGilEnabled()) {
    return nullptr;
  }
  const NajaPythonRuntimeAPI *api = nullptr;
  if (!getRuntimeAPI(api)) {
    return nullptr;
  }

  PyObject *suppliedArguments = nullptr;
  if (!PyArg_ParseTuple(args, "O:run", &suppliedArguments)) {
    return nullptr;
  }

  try {
    if (PyUnicode_Check(suppliedArguments) ||
        PyBytes_Check(suppliedArguments) ||
        PyByteArray_Check(suppliedArguments)) {
      PyErr_SetString(
          PyExc_TypeError,
          "run() expects a sequence of arguments, not a single path");
      return nullptr;
    }
    OwnedPyObject sequence(
        PySequence_Fast(suppliedArguments,
                        "run() expects a sequence of command-line arguments"));
    if (sequence == nullptr) {
      return nullptr;
    }
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(sequence.get());
    if (count >= INT_MAX) {
      PyErr_SetString(PyExc_OverflowError, "too many Kepler Formal arguments");
      return nullptr;
    }

    std::vector<std::string> argumentStorage;
    argumentStorage.reserve(static_cast<size_t>(count) + 1);
    argumentStorage.emplace_back("kepler-formal-python");
    PyObject **items = PySequence_Fast_ITEMS(sequence.get());
    for (Py_ssize_t i = 0; i < count; ++i) {
      if (!appendArgument(items[i], argumentStorage)) {
        return nullptr;
      }
    }

    std::vector<char *> argv;
    argv.reserve(argumentStorage.size());
    for (auto &argument : argumentStorage) {
      argv.push_back(argument.data());
    }

    KEPLER_FORMAL::RunResult result;
    KEPLER_FORMAL::Config::ScopedVerificationContext verificationContext;
    KEPLER_FORMAL::runKeplerFormal(static_cast<int>(argv.size()), argv.data(),
                                   result);
    return resultToDictionary(result);
  } catch (const std::exception &error) {
    PyErr_SetString(PyExc_RuntimeError, error.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "unknown native Kepler Formal failure");
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
    {"run", run, METH_VARARGS,
     "Run Kepler Formal in process and return an owning result dictionary."},
    {"from_najaeda", fromNajaeda, METH_VARARGS,
     "Capture a live NajaEDA SNLDesign without copying its native netlist."},
    {"verify_designs", verifyDesigns, METH_VARARGS,
     "Verify two captured NajaEDA designs without serialization or copying."},
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
  const NajaPythonRuntimeAPI *api = nullptr;
  if (!getRuntimeAPI(api, PyExc_ImportError)) {
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
  return result.release();
}
