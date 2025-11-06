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
  Op op;
  size_t varId;
  std::shared_ptr<BoolExpr>
      l;  // raw pointer — not owning; use index/ptr identity for the key
  std::shared_ptr<BoolExpr> r;  // raw pointer
};

class BoolExprCache {
 public:
  using Key = BoolExprCacheKey;

  // Lookup-or-create API
  static std::shared_ptr<BoolExpr> getExpression(Key const& k);
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
