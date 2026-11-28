// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "BoolExpr.h"

// Boolean connectives that fold constants. BoolExpr nodes are hash-consed, so
// the two constants compare by pointer.
namespace KEPLER_FORMAL::SEC::PropertyLogic {

inline BoolExpr* constTrue() { return BoolExpr::createTrue(); }
inline BoolExpr* constFalse() { return BoolExpr::createFalse(); }

inline BoolExpr* notExpr(BoolExpr* a) {
  if (a == constTrue()) {
    return constFalse();
  }
  if (a == constFalse()) {
    return constTrue();
  }
  return BoolExpr::Not(a);
}

inline BoolExpr* andExpr(BoolExpr* a, BoolExpr* b) {
  if (a == constFalse() || b == constFalse()) {
    return constFalse();
  }
  if (a == constTrue()) {
    return b;
  }
  if (b == constTrue()) {
    return a;
  }
  return BoolExpr::And(a, b);
}

inline BoolExpr* orExpr(BoolExpr* a, BoolExpr* b) {
  if (a == constTrue() || b == constTrue()) {
    return constTrue();
  }
  if (a == constFalse()) {
    return b;
  }
  if (b == constFalse()) {
    return a;
  }
  return BoolExpr::Or(a, b);
}

inline BoolExpr* xorExpr(BoolExpr* a, BoolExpr* b) {
  if (a == constFalse()) {
    return b;
  }
  if (b == constFalse()) {
    return a;
  }
  if (a == constTrue()) {
    return notExpr(b);
  }
  if (b == constTrue()) {
    return notExpr(a);
  }
  return BoolExpr::Xor(a, b);
}

}  // namespace KEPLER_FORMAL::SEC::PropertyLogic
