// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <climits>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "KeplerFormalDriver.h"

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

PyObject *run(PyObject *, PyObject *args) {
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

PyMODINIT_FUNC PyInit__native() { return PyModule_Create(&module); }
