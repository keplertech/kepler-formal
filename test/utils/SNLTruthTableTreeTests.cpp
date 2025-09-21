// SNLTruthTableTree_test.cpp

#include "SNLTruthTableTree.h"
#include "SNLTruthTable.h"

#include <gtest/gtest.h>
#include <bitset>
#include <memory>
#include <vector>
#include <stdexcept>

using namespace naja::NL;
using namespace KEPLER_FORMAL;

// alias the unified Node type
using Node = SNLTruthTableTree::Node;

//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

// Build a small truth‐table via the "mask" constructor (valid for size ≤ 6).
//   size  mask
//  -----  ------------
//    1    b1 b0   (bit0=LSB)
//    2    b3 b2 b1 b0
//   etc.
static SNLTruthTable makeMaskTable(uint32_t size, uint64_t mask) {
  return SNLTruthTable(size, mask);
}

//------------------------------------------------------------------------------
// Leaf (formerly InputNode) tests
//------------------------------------------------------------------------------

TEST(InputNodeTest, ReturnsCorrectValue) {
  std::vector<bool> inputs{false, true, false};
  Node leaf(1);  // inputIndex = 1

  EXPECT_TRUE (leaf.eval(inputs));
  inputs[1] = false;
  EXPECT_FALSE(leaf.eval(inputs));
}

TEST(InputNodeTest, ThrowsIfIndexOutOfRange) {
  std::vector<bool> inputs{true, false};
  Node leaf(2);  // inputIndex = 2
  EXPECT_THROW(leaf.eval(inputs), std::out_of_range);
}

//------------------------------------------------------------------------------
// Table (formerly TableNode) tests
//------------------------------------------------------------------------------

TEST(TableNodeTest, AndGateLogic) {
  auto andTable = makeMaskTable(/*size=*/2, /*mask=*/0b1000);
  Node tbl(andTable);
  tbl.addChild(std::make_unique<Node>(0));
  tbl.addChild(std::make_unique<Node>(1));

  EXPECT_FALSE(tbl.eval({false, false}));
  EXPECT_FALSE(tbl.eval({false, true }));
  EXPECT_FALSE(tbl.eval({true,  false}));
  EXPECT_TRUE (tbl.eval({true,  true }));
}

TEST(TableNodeTest, NotGateLogic) {
  auto notTable = makeMaskTable(/*size=*/1, /*mask=*/0b01);
  Node tbl(notTable);
  tbl.addChild(std::make_unique<Node>(0));

  EXPECT_TRUE (tbl.eval({false}));
  EXPECT_FALSE(tbl.eval({true }));
}

TEST(TableNodeTest, ThrowsOnTableSizeMismatch) {
  auto tinyTable = makeMaskTable(/*size=*/1, /*mask=*/0b01);
  Node tbl(tinyTable);
  tbl.addChild(std::make_unique<Node>(0));
  tbl.addChild(std::make_unique<Node>(1));

  EXPECT_THROW(tbl.eval({false, true}), std::logic_error);
}

//------------------------------------------------------------------------------
// SNLTruthTableTree tests
//------------------------------------------------------------------------------

// Compose (AND → NOT) to get NAND
TEST(SNLTruthTableTreeTest, ComposeAndNotIsNand) {
  // AND subtree
  auto andTbl = makeMaskTable(2, /*0b1000*/8);
  auto andN = std::make_unique<Node>(andTbl);
  andN->addChild(std::make_unique<Node>(0));
  andN->addChild(std::make_unique<Node>(1));

  // NOT wrapper
  auto notTbl = makeMaskTable(1, /*0b01*/1);
  auto notN   = std::make_unique<Node>(notTbl);
  notN->addChild(std::move(andN));

  // full tree expects 2 external inputs
  SNLTruthTableTree tree(std::move(notN), /*numExternalInputs=*/2);

  struct Case { bool a, b, out; };
  std::vector<Case> cases = {
    {false, false, true},
    {false, true,  true},
    {true,  false, true},
    {true,  true,  false},
  };

  for (auto c : cases) {
    EXPECT_EQ(tree.eval({c.a, c.b}), c.out)
      << "a=" << c.a << " b=" << c.b;
  }
}

TEST(SNLTruthTableTreeTest, ThrowsOnWrongExternalSize) {
  auto bufTbl = makeMaskTable(1, /*0b10*/2);
  auto bufN   = std::make_unique<Node>(bufTbl);
  bufN->addChild(std::make_unique<Node>(0));

  SNLTruthTableTree tree(std::move(bufN), /*numExternalInputs=*/1);

  EXPECT_THROW(tree.eval({}),             std::invalid_argument);
  EXPECT_THROW(tree.eval({true, false}),  std::invalid_argument);
}

//------------------------------------------------------------------------------
// Dynamic child-addition tests
//------------------------------------------------------------------------------

TEST(TableNodeDynamicTest, ThreeInputOrLogic) {
  auto or3Table = makeMaskTable(/*size=*/3, /*mask=*/0b11111110);
  Node tbl(or3Table);

  tbl.addChild(std::make_unique<Node>(0));
  tbl.addChild(std::make_unique<Node>(1));
  tbl.addChild(std::make_unique<Node>(2));

  for (uint32_t i = 0; i < (1u << 3); ++i) {
    bool a = (i >> 0) & 1;
    bool b = (i >> 1) & 1;
    bool c = (i >> 2) & 1;
    bool expected = a || b || c;
    EXPECT_EQ(tbl.eval({a, b, c}), expected)
        << "failed OR3 for bits=" << std::bitset<3>(i);
  }
}

TEST(TableNodeDynamicTest, TwoOfThreeThresholdLogic) {
  auto thrTable = makeMaskTable(/*size=*/3, /*mask=*/0b11101000);
  Node tbl(thrTable);

  tbl.addChild(std::make_unique<Node>(0));
  tbl.addChild(std::make_unique<Node>(1));
  tbl.addChild(std::make_unique<Node>(2));

  for (uint32_t i = 0; i < (1u << 3); ++i) {
    int count = ((i>>0)&1) + ((i>>1)&1) + ((i>>2)&1);
    bool expected = (count >= 2);
    bool a = (i >> 0) & 1;
    bool b = (i >> 1) & 1;
    bool c = (i >> 2) & 1;
    EXPECT_EQ(tbl.eval({a, b, c}), expected)
        << "failed threshold2/3 for bits=" << std::bitset<3>(i);
  }
}

//------------------------------------------------------------------------------
// Pyramid‐of‐And‐Gates test (8→4→2→1)
//------------------------------------------------------------------------------

TEST(TableNodePyramidTest, EightInputAndPyramid) {
  auto andTbl = makeMaskTable(/*size=*/2, /*mask=*/0b1000);

  auto and01 = std::make_unique<Node>(andTbl);
  and01->addChild(std::make_unique<Node>(0));
  and01->addChild(std::make_unique<Node>(1));

  auto and23 = std::make_unique<Node>(andTbl);
  and23->addChild(std::make_unique<Node>(2));
  and23->addChild(std::make_unique<Node>(3));

  auto and45 = std::make_unique<Node>(andTbl);
  and45->addChild(std::make_unique<Node>(4));
  and45->addChild(std::make_unique<Node>(5));

  auto and67 = std::make_unique<Node>(andTbl);
  and67->addChild(std::make_unique<Node>(6));
  and67->addChild(std::make_unique<Node>(7));

  auto and0123 = std::make_unique<Node>(andTbl);
  and0123->addChild(std::move(and01));
  and0123->addChild(std::move(and23));

  auto and4567 = std::make_unique<Node>(andTbl);
  and4567->addChild(std::move(and45));
  and4567->addChild(std::move(and67));

  auto root = std::make_unique<Node>(andTbl);
  root->addChild(std::move(and0123));
  root->addChild(std::move(and4567));

  SNLTruthTableTree tree(std::move(root), /*numExternalInputs=*/8);

  for (uint32_t mask = 0; mask < (1u << 8); ++mask) {
    std::vector<bool> in(8);
    for (int i = 0; i < 8; ++i) {
      in[i] = ((mask >> i) & 1) != 0;
    }
    bool expected = (mask == 0xFF);
    EXPECT_EQ(tree.eval(in), expected)
        << "mask=" << std::bitset<8>(mask);
  }
}

// Build a tree via concat() operations
TEST(SNLTruthTableTreeTest, Pyramid4Ands2Ors1NorViaConcat) {
  auto andTbl = makeMaskTable(2, 0b1000);
  auto  orTbl = makeMaskTable(2, 0b1110);
  auto norTbl = makeMaskTable(2, 0b0001);

  SNLTruthTableTree tree(
    std::make_unique<Node>(0), // start from leaf 0
    /*numExternalInputs=*/1
  );

  // Step 1: wrap A in a 2-input NOR → now expects A,B
  tree.concat(/*borderIndex=*/0, norTbl);

  // graft left OR
  tree.concat(0, orTbl);
  // graft right OR
  tree.concat(1, orTbl);

  // graft four ANDs
  for (size_t slot = 0; slot < 4; ++slot)
    tree.concat(slot, andTbl);

  for (uint32_t mask = 0; mask < (1u << 8); ++mask) {
    std::vector<bool> in(8);
    for (int i = 0; i < 8; ++i)
      in[i] = ((mask >> i) & 1) != 0;

    bool a0 = in[0] && in[4];
    bool a1 = in[1] && in[5];
    bool a2 = in[2] && in[6];
    bool a3 = in[3] && in[7];
    bool o0 = a0 || a2;
    bool o1 = a1 || a3;
    bool expected = !(o0 || o1);

    EXPECT_EQ(tree.eval(in), expected)
        << "mask=" << std::bitset<8>(mask);
  }
}

// Test concatFull
TEST(SNLTruthTableTreeTest, ConcatFullInvertsBothLeaves) {
  auto orTbl  = makeMaskTable(2, /*mask=*/0b1110);
  auto notTbl = makeMaskTable(1, /*mask=*/0b01);

  auto root = std::make_unique<Node>(orTbl);
  root->addChild(std::make_unique<Node>(0));
  root->addChild(std::make_unique<Node>(1));

  SNLTruthTableTree tree(std::move(root), /*numExternalInputs=*/2);

  tree.concatFull({ notTbl, notTbl });

  struct Case { bool a, b, out; };
  std::vector<Case> cases = {
    {false, false, true},
    {false, true,  true},
    {true,  false, true},
    {true,  true,  false},
  };

  for (auto c : cases) {
    EXPECT_EQ(tree.eval({c.a, c.b}), c.out)
        << "a=" << c.a << " b=" << c.b;
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
