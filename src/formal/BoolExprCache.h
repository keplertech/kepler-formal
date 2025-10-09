#pragma once

#include <memory>
#include <atomic>
#include <cstddef>

namespace KEPLER_FORMAL {

class BoolExpr;

enum class Op { VAR, AND, OR, NOT, XOR, NONE };

// Minimal POD representing the cache query. Use raw pointers for children to avoid inclusion cycles.
struct BoolExprCacheKey {
  Op op;
  size_t varId;
  BoolExpr* l;   // raw pointer — not owning; use index/ptr identity for the key
  BoolExpr* r;   // raw pointer
};

class BoolExprCache {
public:
  using Key = BoolExprCacheKey;

  // Lookup-or-create API
  static BoolExpr* getExpression(Key const& k);
  static void destroy();
private:
  struct Impl;
  static Impl& impl();
  // destructor that will delete all stored BoolExpr*
  static std::atomic<size_t> lastID_;
};

} // namespace KEPLER_FORMAL
