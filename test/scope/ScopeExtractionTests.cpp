// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: Apache-2.0
//
// Unit tests for ScopeExtraction::collectVerificationScopes().
//
// The two designs of a comparison live in separate libraries, as the CLI
// loads them, so modules, instances and nets of the two sides carry the same
// names. Primitive models are either shared (one primitives library) or
// duplicated per side to exercise name-based leaf comparison.

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "ScopeExtraction.h"

#include "NLLibraryTruthTables.h"
#include "NLUniverse.h"
#include "SNLBusNet.h"
#include "SNLBusNetBit.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstParameter.h"
#include "SNLInstance.h"
#include "SNLParameter.h"
#include "SNLSVConstructor.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"

using namespace naja;
using namespace naja::NL;

// Test helper derived from ScopeExtraction to expose internals for assertions.
class TestScopeExtraction : public ScopeExtraction {
 public:
  TestScopeExtraction() : ScopeExtraction() {}
  void setTop0(naja::NL::SNLDesign* d) { top0_ = d; }
  void setTop1(naja::NL::SNLDesign* d) { top1_ = d; }

  const std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>>&
  getDesignsToVerify() const {
    return designsToVerify_;
  }

  void runCollect() { collectVerificationScopes(); }
};

namespace {

using ScopePair = std::pair<SNLDesign*, SNLDesign*>;

// One side of a comparison: a standard library for the modules and, when the
// test wants per-side primitives, a primitives library of its own.
struct Side {
  NLDB* db = nullptr;
  NLLibrary* designs = nullptr;
  NLLibrary* primitives = nullptr;
};

Side createSide(NLUniverse* universe, const std::string& suffix, bool ownPrimitives) {
  Side side;
  side.db = NLDB::create(universe);
  side.designs = NLLibrary::create(side.db, NLLibrary::Type::Standard,
                                   NLName("designs" + suffix));
  if (ownPrimitives) {
    side.primitives = NLLibrary::create(side.db, NLLibrary::Type::Primitives,
                                        NLName("primitives" + suffix));
  }
  return side;
}

// A one-output primitive with `inputs` inputs and the given truth table.
SNLDesign* createPrimitive(NLLibrary* primitives,
                           const std::string& name,
                           size_t inputs,
                           uint64_t table) {
  auto* model = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName(name));
  for (size_t i = 0; i < inputs; ++i) {
    SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("i" + std::to_string(i)));
  }
  SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("o"));
  SNLDesignModeling::setTruthTable(
      model, SNLTruthTable(static_cast<uint32_t>(inputs), table,
                           SNLTruthTable::fullDependencies(static_cast<uint32_t>(inputs))));
  return model;
}

// Shared primitives: a 2-input AND, a 2-input OR and an inverter.
struct Primitives {
  SNLDesign* and2 = nullptr;
  SNLDesign* or2 = nullptr;
  SNLDesign* inv = nullptr;
};

Primitives createPrimitives(NLLibrary* primitives) {
  Primitives p;
  p.and2 = createPrimitive(primitives, "and2", 2, 0b1000);
  p.or2 = createPrimitive(primitives, "or2", 2, 0b1110);
  p.inv = createPrimitive(primitives, "inv", 1, 0b01);
  NLLibraryTruthTables::construct(primitives);
  return p;
}

SNLScalarTerm* input(SNLDesign* design, const std::string& name) {
  return SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName(name));
}

SNLScalarTerm* output(SNLDesign* design, const std::string& name) {
  return SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName(name));
}

SNLScalarNet* net(SNLDesign* design, const std::string& name) {
  return SNLScalarNet::create(design, NLName(name));
}

void connect(SNLInstance* instance, const std::string& pin, SNLBitNet* bitNet) {
  auto* term = instance->getModel()->getScalarTerm(NLName(pin));
  ASSERT_NE(term, nullptr) << pin;
  instance->getInstTerm(term)->setNet(bitNet);
}

// A gate module: o = gate(a, b) with the gate chosen by the caller.
//   module gate_mod(input a, input b, output o);
//     <gate> u_g(.i0(a), .i1(b), .o(o));
SNLDesign* createGateModule(NLLibrary* designs, const std::string& name, SNLDesign* gate) {
  auto* design = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName(name));
  auto* a = input(design, "a");
  auto* b = input(design, "b");
  auto* o = output(design, "o");
  auto* na = net(design, "a");
  auto* nb = net(design, "b");
  auto* no = net(design, "o");
  a->setNet(na);
  b->setNet(nb);
  o->setNet(no);
  auto* g = SNLInstance::create(design, gate, NLName("u_g"));
  connect(g, "i0", na);
  connect(g, "i1", nb);
  connect(g, "o", no);
  return design;
}

// A top with two instances of `child` and an inverter:
//   module top(input x, input y, output p, output q);
//     child u0(.a(x), .b(y), .o(p));
//     child u1(.a(y), .b(x), .o(t));
//     inv   u_inv(.i0(t), .o(q));
SNLDesign* createTop(NLLibrary* designs,
                     const std::string& name,
                     SNLDesign* child,
                     SNLDesign* inv,
                     bool swapSecondInstance = false) {
  auto* design = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName(name));
  auto* x = input(design, "x");
  auto* y = input(design, "y");
  auto* p = output(design, "p");
  auto* q = output(design, "q");
  auto* nx = net(design, "x");
  auto* ny = net(design, "y");
  auto* np = net(design, "p");
  auto* nq = net(design, "q");
  auto* nt = net(design, "t");
  x->setNet(nx);
  y->setNet(ny);
  p->setNet(np);
  q->setNet(nq);
  auto* u0 = SNLInstance::create(design, child, NLName("u0"));
  connect(u0, "a", nx);
  connect(u0, "b", ny);
  connect(u0, "o", np);
  auto* u1 = SNLInstance::create(design, child, NLName("u1"));
  connect(u1, "a", swapSecondInstance ? nx : ny);
  connect(u1, "b", swapSecondInstance ? ny : nx);
  connect(u1, "o", nt);
  auto* uInv = SNLInstance::create(design, inv, NLName("u_inv"));
  connect(uInv, "i0", nt);
  connect(uInv, "o", nq);
  return design;
}

std::set<ScopePair> collect(SNLDesign* top0, SNLDesign* top1, size_t* identical = nullptr) {
  ScopeExtraction extraction(top0, top1);
  extraction.collectVerificationScopes();
  if (identical != nullptr) {
    *identical = extraction.getIdenticalModelPairCount();
  }
  return extraction.getScopesToVerify();
}

SNLDesign* loadSystemVerilog(NLLibrary* designs,
                             const std::string& tag,
                             const std::string& topName,
                             const std::string& sourceText) {
  const auto dir = std::filesystem::temp_directory_path() / ("scope_sv_" + tag);
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  const auto path = dir / (topName + ".sv");
  {
    std::ofstream file(path);
    file << sourceText;
  }
  SNLSVConstructor constructor(designs);
  SNLSVConstructor::ConstructOptions options;
  options.diagnosticsReportPath.reset();
  constructor.construct(path, options);
  auto* top = designs->getSNLDesign(NLName(topName));
  if (top == nullptr) {
    throw std::runtime_error("SystemVerilog top `" + topName + "` was not constructed");
  }
  return top;
}

// Local creation-order IDs of the instances and nets of a module, which the
// lowering must reproduce for an unchanged module.
std::vector<std::string> localShape(SNLDesign* design) {
  std::vector<std::string> shape;
  for (auto* instance : design->getInstances()) {
    shape.push_back("i" + std::to_string(instance->getID()) + ":" +
                    instance->getModel()->getName().getString() + ":" +
                    instance->getName().getString());
  }
  for (auto* bitNet : design->getBitNets()) {
    shape.push_back("n" + std::to_string(bitNet->getID()) + ":" +
                    bitNet->getName().getString());
  }
  return shape;
}

}  // namespace

class ScopeExtractionUnitTests : public ::testing::Test {
 protected:
  void SetUp() override {
    universe_ = NLUniverse::create();
    // Each side keeps its own primitives, as the CLI's per-design databases
    // do: an instance must not outlive the database of its model.
    side0_ = createSide(universe_, "0", true);
    side1_ = createSide(universe_, "1", true);
    prims0_ = createPrimitives(side0_.primitives);
    prims1_ = createPrimitives(side1_.primitives);
  }

  void TearDown() override {
    NLUniverse::get()->destroy();
    universe_ = nullptr;
  }

  NLUniverse* universe_ = nullptr;
  Primitives prims0_;
  Primitives prims1_;
  Side side0_;
  Side side1_;
};

// Identical two-level hierarchies: nothing to verify, every pair identical.
TEST_F(ScopeExtractionUnitTests, IdenticalHierarchies_NoScopes) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.and2);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  size_t identical = 0;
  const auto scopes = collect(top0, top1, &identical);
  EXPECT_TRUE(scopes.empty());
  // top and child (compared once although instantiated twice); shared leaf
  // models are decided inside their parent and not counted as pairs.
  EXPECT_EQ(identical, 2u);
}

// The legacy entry point used by the CLI: a default-constructed extraction
// with the tops set afterwards behaves the same.
TEST_F(ScopeExtractionUnitTests, TestSubclassEntryPoint) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.and2);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  TestScopeExtraction se;
  se.setTop0(top0);
  se.setTop1(top1);
  se.runCollect();
  EXPECT_TRUE(se.getDesignsToVerify().empty());
}

// A different instance count at the top makes the top pair the scope.
TEST_F(ScopeExtractionUnitTests, DifferentInstanceCount_TopIsScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.and2);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);
  // One more inverter on side 1, loadless but present.
  SNLInstance::create(top1, prims1_.inv, NLName("u_extra"));

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(top0, top1));
}

// Same instances, swapped wiring of one instance: the top's own
// implementation differs, so the top is the scope and the child is not.
TEST_F(ScopeExtractionUnitTests, SwappedWiring_TopIsScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.and2);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv, /*swap=*/true);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(top0, top1));
}

// Only the child's gate changed (same interface, same wiring): the child
// pair is the single scope although it is instantiated twice, and the top,
// whose own wiring is identical, is not a scope.
TEST_F(ScopeExtractionUnitTests, ChangedLeafModel_ChildIsTheOnlyScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.or2);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(child0, child1));
}

// A child with an extra internal net and gate is the scope; the parent stays
// out because its own connections are unchanged.
TEST_F(ScopeExtractionUnitTests, ChangedChildBody_ChildIsTheOnlyScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.and2);
  // Side 1: o = and(a, b) still, plus an unused inverter of a.
  auto* extra = net(child1, "a_n");
  auto* uInv = SNLInstance::create(child1, prims1_.inv, NLName("u_extra"));
  connect(uInv, "i0", child1->getScalarNet(NLName("a")));
  connect(uInv, "o", extra);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(child0, child1));
}

// When the child's interface changes, the pair cannot be verified alone and
// the parent, whose connections to it differ, becomes the scope.
TEST_F(ScopeExtractionUnitTests, ChangedChildInterface_ParentIsScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.and2);
  auto* enable = input(child1, "en");
  enable->setNet(net(child1, "en"));
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(top0, top1));
}

// Parent wiring and child both differ: the parent is verified as a whole,
// the child is not listed separately.
TEST_F(ScopeExtractionUnitTests, OwnAndChildDiffer_OnlyParentIsScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = createGateModule(side1_.designs, "child", prims1_.or2);
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv, /*swap=*/true);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(top0, top1));
}

// Leaf models from separate primitive libraries pair by name, interface and
// parameters, never by model ID.
TEST_F(ScopeExtractionUnitTests, LeafModelsPairByNameAcrossLibraries) {
  Side own0 = createSide(universe_, "p0", true);
  Side own1 = createSide(universe_, "p1", true);
  // Create the primitives in different orders so their IDs differ.
  auto* inv0 = createPrimitive(own0.primitives, "inv", 1, 0b01);
  auto* and0 = createPrimitive(own0.primitives, "and2", 2, 0b1000);
  NLLibraryTruthTables::construct(own0.primitives);
  auto* and1 = createPrimitive(own1.primitives, "and2", 2, 0b1000);
  auto* inv1 = createPrimitive(own1.primitives, "inv", 1, 0b01);
  NLLibraryTruthTables::construct(own1.primitives);
  ASSERT_NE(and0->getID(), and1->getID());

  auto* child0 = createGateModule(own0.designs, "child", and0);
  auto* child1 = createGateModule(own1.designs, "child", and1);
  auto* top0 = createTop(own0.designs, "top", child0, inv0);
  auto* top1 = createTop(own1.designs, "top", child1, inv1);

  EXPECT_TRUE(collect(top0, top1).empty());

  // A parameter on one leaf model breaks the identity.
  SNLParameter::create(and1, NLName("DRIVE"), SNLParameter::Type::Decimal, "2");
  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(child0, child1));
}

// An instance parameter value is part of the module's own implementation.
TEST_F(ScopeExtractionUnitTests, InstanceParameterDiffers_ModuleIsScope) {
  const auto createLut = [](NLLibrary* primitives) {
    auto* lut = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName("lut1"));
    input(lut, "i0");
    input(lut, "i1");
    output(lut, "o");
    SNLParameter::create(lut, NLName("INIT"), SNLParameter::Type::Binary, "2'h2");
    return lut;
  };
  auto* lut0 = createLut(side0_.primitives);
  auto* lut1 = createLut(side1_.primitives);

  auto* child0 = createGateModule(side0_.designs, "child", lut0);
  auto* child1 = createGateModule(side1_.designs, "child", lut1);
  SNLInstParameter::create(child0->getInstance(NLName("u_g")),
                           lut0->getParameter(NLName("INIT")), "2'h2");
  SNLInstParameter::create(child1->getInstance(NLName("u_g")),
                           lut1->getParameter(NLName("INIT")), "2'h1");
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(child0, child1));
}

// A renamed instance is a different implementation: identity is strict, a
// false "identical" would be unsound while a false "different" only costs a
// proof.
TEST_F(ScopeExtractionUnitTests, RenamedInstance_ModuleIsScope) {
  auto* child0 = createGateModule(side0_.designs, "child", prims0_.and2);
  auto* child1 = SNLDesign::create(side1_.designs, SNLDesign::Type::Standard, NLName("child"));
  {
    auto* a = input(child1, "a");
    auto* b = input(child1, "b");
    auto* o = output(child1, "o");
    auto* na = net(child1, "a");
    auto* nb = net(child1, "b");
    auto* no = net(child1, "o");
    a->setNet(na);
    b->setNet(nb);
    o->setNet(no);
    auto* g = SNLInstance::create(child1, prims1_.and2, NLName("u_gate"));
    connect(g, "i0", na);
    connect(g, "i1", nb);
    connect(g, "o", no);
  }
  auto* top0 = createTop(side0_.designs, "top", child0, prims0_.inv);
  auto* top1 = createTop(side1_.designs, "top", child1, prims1_.inv);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(child0, child1));
}

// Bus ports and nets compare bit by bit.
TEST_F(ScopeExtractionUnitTests, BusWiring_SwappedBitsIsScope) {
  const auto build = [&](NLLibrary* designs, SNLDesign* gate, bool swap) {
    auto* design = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName("bus_top"));
    auto* in = SNLBusTerm::create(design, SNLTerm::Direction::Input, 1, 0, NLName("d"));
    auto* o = output(design, "o");
    auto* nd = SNLBusNet::create(design, 1, 0, NLName("d"));
    auto* no = net(design, "o");
    in->getBit(0)->setNet(nd->getBit(0));
    in->getBit(1)->setNet(nd->getBit(1));
    o->setNet(no);
    auto* g = SNLInstance::create(design, gate, NLName("u_g"));
    connect(g, "i0", nd->getBit(swap ? 1 : 0));
    connect(g, "i1", nd->getBit(swap ? 0 : 1));
    connect(g, "o", no);
    return design;
  };
  auto* same0 = build(side0_.designs, prims0_.and2, false);
  auto* same1 = build(side1_.designs, prims1_.and2, false);
  EXPECT_TRUE(collect(same0, same1).empty());

  Side extra0 = createSide(universe_, "b0", true);
  Side extra1 = createSide(universe_, "b1", true);
  auto* top0 = build(extra0.designs, createPrimitives(extra0.primitives).and2, false);
  auto* top1 = build(extra1.designs, createPrimitives(extra1.primitives).and2, true);
  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(top0, top1));
}

// Elaborated SystemVerilog: anonymous lowered gates and nets of an unchanged
// module get the same creation-order IDs on both sides, so only the module
// whose body changed is a scope, and the parents are not.
TEST_F(ScopeExtractionUnitTests, SystemVerilogChangedSubmodule_IsTheOnlyScope) {
  const auto source = [](const std::string& backBody) {
    return std::string(R"(
module sv_front(
  input  logic       clk,
  input  logic       rst_n,
  input  logic [3:0] k,
  output logic       armed
);
  logic [2:0] c0;
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      armed <= 1'b0;
      c0    <= 3'd0;
    end else begin
      if (c0 != 3'd6) c0 <= c0 + 3'd1;
      if (c0 == 3'd6 && k == 4'h5) armed <= 1'b1;
    end
  end
endmodule
module sv_back(
  input  logic       clk,
  input  logic       rst_n,
  input  logic [3:0] k,
  input  logic       armed,
  output logic       fire
);
  logic [2:0] c1;
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      c1   <= 3'd0;
      fire <= 1'b0;
    end else begin
      if (armed && c1 != 3'd6) c1 <= c1 + 3'd1;
      )") + backBody + R"(
    end
  end
endmodule
module sv_mid(
  input  logic       clk,
  input  logic       rst_n,
  input  logic [3:0] k,
  output logic       fire
);
  logic armed;
  sv_front u_front(.clk(clk), .rst_n(rst_n), .k(k), .armed(armed));
  sv_back  u_back(.clk(clk), .rst_n(rst_n), .k(k), .armed(armed), .fire(fire));
endmodule
module sv_top(
  input  logic       clk,
  input  logic       rst_n,
  input  logic [3:0] k,
  output logic       q0,
  output logic       q1
);
  sv_mid u_mid0(.clk(clk), .rst_n(rst_n), .k(k), .fire(q0));
  sv_mid u_mid1(.clk(clk), .rst_n(rst_n), .k(~k), .fire(q1));
endmodule
)";
  };
  auto* top0 = loadSystemVerilog(side0_.designs, "before", "sv_top",
                                 source("fire <= fire | (c1 == 3'd6 && k == 4'h3);"));
  auto* top1 = loadSystemVerilog(
      side1_.designs, "after", "sv_top",
      source("fire <= fire | (k == 4'h3 && (c1 == 3'd6 || c1 == 3'd5));"));

  auto* front0 = side0_.designs->getSNLDesign(NLName("sv_front"));
  auto* front1 = side1_.designs->getSNLDesign(NLName("sv_front"));
  auto* back0 = side0_.designs->getSNLDesign(NLName("sv_back"));
  auto* back1 = side1_.designs->getSNLDesign(NLName("sv_back"));
  auto* mid0 = side0_.designs->getSNLDesign(NLName("sv_mid"));
  auto* mid1 = side1_.designs->getSNLDesign(NLName("sv_mid"));
  ASSERT_NE(front0, nullptr);
  ASSERT_NE(back1, nullptr);
  ASSERT_NE(mid1, nullptr);

  // Determinism of the lowering for the unchanged modules.
  EXPECT_EQ(localShape(front0), localShape(front1));
  EXPECT_EQ(localShape(mid0), localShape(mid1));
  EXPECT_EQ(localShape(top0), localShape(top1));
  EXPECT_NE(localShape(back0), localShape(back1));
  EXPECT_GT(front0->getInstances().size(), 2u);

  const auto scopes = collect(top0, top1);
  ASSERT_EQ(scopes.size(), 1u);
  EXPECT_EQ(*scopes.begin(), ScopePair(back0, back1));
}

// Required main function for Google Test
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
