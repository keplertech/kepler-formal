#pragma once
#include <cstddef>
#include <memory>

namespace KEPLER_FORMAL {

class BoolExpr;  // forward declaration
enum class Op;   // forward declaration

class BoolExprCache {
 public:
  struct Key {
    Op op;
    size_t varId;
    BoolExpr* l;
    BoolExpr* r;
  };

  static std::shared_ptr<BoolExpr> getExpression(Key const& k);

 private:
  // declare opaque alias names only; full types defined in cpp
  struct Impl;
  static Impl& impl();
};

}  // namespace KEPLER_FORMAL
