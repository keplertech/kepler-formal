// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

// XLS tracks newer Clang headers where TypeBase.h is split out. Ubuntu 22.04
// ships Clang 14, where those declarations still live in Type.h. Keep this
// compatibility shim in Kepler-owned include space instead of patching XLS.
#if defined(__has_include_next)
#  if __has_include_next("clang/AST/TypeBase.h")
#    include_next "clang/AST/TypeBase.h"
#  elif __has_include("clang/AST/Type.h")
#    include "clang/AST/Type.h"
#  else
#    error "Missing Clang AST type headers"
#  endif
#elif defined(__has_include)
#  if __has_include("clang/AST/Type.h")
#    include "clang/AST/Type.h"
#  else
#    error "Missing Clang AST type headers"
#  endif
#else
#  include "clang/AST/Type.h"
#endif
