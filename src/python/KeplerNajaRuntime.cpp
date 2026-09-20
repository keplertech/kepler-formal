// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "KeplerNajaRuntime.h"

#include <cstring>

#include "NajaPythonRuntimeAPI.h"
#include "NajaRuntimeBuild.h"
#include "NLUniverse.h"

namespace KEPLER_FORMAL {
namespace {

const NajaPythonRuntimeAPI* getRuntimeAPI(PyObject* errorType) {
  const auto* api = NajaPythonRuntime_Import();
  if (api == nullptr) {
    return nullptr;
  }
  if (api->build_id == nullptr || api->runtime_identity == nullptr ||
      api->get_universe == nullptr || api->unwrap_design == nullptr) {
    PyErr_SetString(errorType, "Incomplete NajaEDA native runtime API");
    return nullptr;
  }
  if (std::strcmp(api->build_id, NAJA_RUNTIME_BUILD_ID) != 0) {
    PyErr_Format(errorType,
                 "NajaEDA native build mismatch: provider %.80s, Kepler %.80s",
                 api->build_id, NAJA_RUNTIME_BUILD_ID);
    return nullptr;
  }
  if (api->runtime_identity != naja::NL::NLUniverse::getRuntimeIdentity()) {
    PyErr_SetString(
        errorType,
        "NajaEDA and Kepler Formal did not load the same native runtime");
    return nullptr;
  }
  if (api->get_universe() != naja::NL::NLUniverse::get()) {
    PyErr_SetString(errorType,
                    "NajaEDA and Kepler Formal disagree on the active universe");
    return nullptr;
  }
  return api;
}

}  // namespace

bool validateNajaRuntime(PyObject* errorType) {
  return getRuntimeAPI(errorType) != nullptr;
}

naja::NL::SNLDesign* unwrapNajaDesign(PyObject* object) {
  const auto* api = getRuntimeAPI(PyExc_RuntimeError);
  return api == nullptr
             ? nullptr
             : static_cast<naja::NL::SNLDesign*>(api->unwrap_design(object));
}

}  // namespace KEPLER_FORMAL
