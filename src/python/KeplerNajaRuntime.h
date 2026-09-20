// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Python.h>

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL {

// All calls require the GIL. The returned design remains owned by NajaEDA.
bool validateNajaRuntime(PyObject* errorType = PyExc_RuntimeError);
naja::NL::SNLDesign* unwrapNajaDesign(PyObject* object);

}  // namespace KEPLER_FORMAL
