// MiterStrategyWithGlucose.cpp

#include "MiterStrategy.h"
#include "BuildPrimaryOutputClauses.h"
#include "SNLDesignModeling.h"
#include "SNLTruthTable2BoolExpr.h"
#include "SNLDesignTruthTable.h"
#include "SNLLogicCloud.h"
#include "NLUniverse.h"
#include "BoolExpr.h"

// include Glucose headers (adjust path to your checkout)
#include "core/Solver.h"
#include "simp/SimpSolver.h"

#include <unordered_map>
#include <string>

using namespace KEPLER_FORMAL;

namespace {
  //
  // A tiny Tseitin‐translator from BoolExpr→Glucose CNF.
  //
  // Returns a Glucose::Lit that stands for `e`, and adds
  // all necessary clauses to S so that Lit ↔ (e) holds.
  //
  // node2var      caches each subformula’s fresh variable index.
  // varName2idx   coalesces all inputs of the same name to one var.
  //
  Glucose::Lit tseitinEncode(
      Glucose::SimpSolver&                       S,
      const std::shared_ptr<BoolExpr>&          e,
      std::unordered_map<const BoolExpr*, int>& node2var,
      std::unordered_map<std::string,int>&      varName2idx)
  {
    // If we’ve already encoded this node, reuse its var
    auto it = node2var.find(e.get());
    if (it != node2var.end())
      return Glucose::mkLit(it->second);

    // Handle leaf‐vars specially to coalesce by name
    if (e->getOp() == Op::VAR) {
      const std::string& name = e->getName();
      int v;
      auto it2 = varName2idx.find(name);
      if (it2 == varName2idx.end()) {
        v = S.newVar();
        varName2idx[name] = v;
      } else {
        v = it2->second;
      }
      node2var[e.get()] = v;
      return Glucose::mkLit(v);
    }

    // Otherwise allocate a fresh var for this gate
    int v = S.newVar();
    Glucose::Lit lit_v = Glucose::mkLit(v);
    node2var[e.get()] = v;

    // Recursively encode children
    auto encodeChild = [&](const std::shared_ptr<BoolExpr>& c){
      return tseitinEncode(S, c, node2var, varName2idx);
    };

    switch (e->getOp()) {
      case Op::NOT: {
        auto a = encodeChild(e->getLeft());
        // v ↔ ¬a  gives: (¬v ∨ ¬a) ∧ (v ∨ a)
        S.addClause(~lit_v, ~a);
        S.addClause( lit_v,  a);
        break;
      }

      case Op::AND: {
        auto a = encodeChild(e->getLeft());
        auto b = encodeChild(e->getRight());
        // v ↔ (a ∧ b)  ===  (¬v ∨ a) ∧ (¬v ∨ b) ∧ (v ∨ ¬a ∨ ¬b)
        S.addClause(~lit_v,  a);
        S.addClause(~lit_v,  b);
        S.addClause( lit_v, ~a, ~b);
        break;
      }

      case Op::OR: {
        auto a = encodeChild(e->getLeft());
        auto b = encodeChild(e->getRight());
        // v ↔ (a ∨ b)  ===  (¬a ∨ v) ∧ (¬b ∨ v) ∧ (¬v ∨ a ∨ b)
        S.addClause(~a,  lit_v);
        S.addClause(~b,  lit_v);
        S.addClause(~lit_v, a, b);
        break;
      }

      case Op::XOR: {
        auto a = encodeChild(e->getLeft());
        auto b = encodeChild(e->getRight());
        // v ↔ (a XOR b)
        S.addClause(~lit_v, ~a, ~b);
        S.addClause(~lit_v,  a,  b);
        S.addClause( lit_v, ~a,  b);
        S.addClause( lit_v,  a, ~b);
        break;
      }

      default:
        // unreachable
        break;
    }

    return lit_v;
  }
}

bool MiterStrategy::run() {
  BuildPrimaryOutputClauses builder;

  // build both sets of POs
  SNLDesign* topInit = NLUniverse::get()->getTopDesign();
  NLUniverse* univ = NLUniverse::get();

  univ->setTopDesign(top0_);
  builder.build();
  auto POs0 = builder.getPOs();

  univ->setTopDesign(top1_);
  builder.build();
  auto POs1 = builder.getPOs();

  univ->setTopDesign(topInit);

  if (POs0.empty() || POs1.empty())
    return false;

  // build the Boolean‐miter expression
  auto miter = buildMiter(POs0, POs1);

  // Now SAT check via Glucose
  Glucose::SimpSolver solver;

  // mappings for Tseitin encoding
  std::unordered_map<const BoolExpr*,int> node2var;
  std::unordered_map<std::string,int>      varName2idx;

  // Tseitin‐encode & get the literal for the root
  Glucose::Lit rootLit = tseitinEncode(
      solver, miter, node2var, varName2idx);

  // Assert root == true
  solver.addClause(rootLit);

  // solve with no assumptions
  bool sat = solver.solve();

  // if UNSAT → miter can never be true → outputs identical
  return !sat;
}

std::shared_ptr<BoolExpr> MiterStrategy::buildMiter(
    const std::vector<std::shared_ptr<BoolExpr>>& A,
    const std::vector<std::shared_ptr<BoolExpr>>& B) const
{
  if (A.size() != B.size()) {
    printf("Miter inputs must match in length: %zu vs %zu\n", A.size(), B.size());
    assert(false && "Miter inputs must match in length");
  }

  // Empty miter = always‐false (no outputs to compare)
  if (A.empty())
    return BoolExpr::createFalse();

  // Start with the first XOR
  auto miter = BoolExpr::Xor(A[0], B[0]);

  // OR in the rest
  for (size_t i = 1, n = A.size(); i < n; ++i) {
    auto diff = BoolExpr::Xor(A[i], B[i]);
    miter = BoolExpr::Or(miter, diff);
  }
  return miter;
}
