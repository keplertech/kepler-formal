/* Copyright 2024-2026 keplertech.io
 * SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#ifdef _WIN32
/* Glucose's serial SolverTypes.h includes pthread.h but uses no pthread
 * declarations. This is deliberately not a pthread emulation: Kepler does
 * not build the POSIX-only Glucose parallel library on Windows. */
#else
#include_next <pthread.h>
#endif
