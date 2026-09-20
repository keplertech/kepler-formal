// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <stddef.h>

#include "clause.h"
#include "watch.h"

// Check the actual solver headers with the same flags as its implementation.
// Unlike Kissat's runtime assertions, these remain active in release wheels.
_Static_assert(sizeof(watch) == sizeof(unsigned),
               "Kissat watches must fit in their unsigned raw storage");
_Static_assert(offsetof(clause, searched) == sizeof(unsigned),
               "Kissat clause bitfields must occupy one unsigned word");
