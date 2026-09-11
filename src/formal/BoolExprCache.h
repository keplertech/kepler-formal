// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <atomic>
#include <cstddef>
#include <memory>

namespace KEPLER_FORMAL {

class BoolExpr;

enum class Op { VAR, AND, OR, NOT, XOR, NONE };

// Minimal POD representing the cache query. Use raw pointers for children to
// avoid inclusion cycles.
struct BoolExprCacheKey {
  Op op = Op::NONE;
  // init with max size_t for invalid 
  size_t varId = (size_t)-1;  // only used if op == VAR; otherwise ignored
  BoolExpr* l  = nullptr;
  BoolExpr* r = nullptr;
};

class BoolExprCache {
 public:
  using Key = BoolExprCacheKey;

  // An isolated cache for a synchronous embedding operation. Existing users
  // retain the original global-cache behavior unless they opt into this scope.
  // Callers must serialize all cache use for the lifetime of the scope.
  class ScopedContext {
   public:
    ScopedContext();
    ~ScopedContext();
    ScopedContext(const ScopedContext&) = delete;
    ScopedContext& operator=(const ScopedContext&) = delete;

   private:
    struct State;
    std::unique_ptr<State> state_;
  };

  // Lookup-or-create API
  static BoolExpr* getExpression(Key const& k);
  static void destroy();

 private:
  struct Impl;
  static Impl& impl();
  // destructor that will delete all stored std::shared_ptr<BoolExpr>
  static std::atomic<size_t> lastID_;
  static size_t numQuaries_;
  static size_t numMiss_;
  static size_t numHit_;
};

}  // namespace KEPLER_FORMAL
