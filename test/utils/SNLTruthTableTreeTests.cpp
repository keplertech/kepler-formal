// SNLTruthTableTree_test.cpp

#include "SNLTruthTableTree.h"
#include <gtest/gtest.h>
#include <bitset>
#include <memory>
#include <vector>
#include <stdexcept>

using namespace naja::NL;
using namespace KEPLER_FORMAL;

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
// InputNode tests
//------------------------------------------------------------------------------

TEST(InputNodeTest, ReturnsCorrectValue) {
  std::vector<bool> inputs{false, true, false};
  InputNode node(1);

  EXPECT_TRUE (node.eval(inputs));
  inputs[1] = false;
  EXPECT_FALSE(node.eval(inputs));
}

TEST(InputNodeTest, ThrowsIfIndexOutOfRange) {
  std::vector<bool> inputs{true, false};
  InputNode node(2);
  EXPECT_THROW(node.eval(inputs), std::out_of_range);
}

//------------------------------------------------------------------------------
// TableNode tests
//------------------------------------------------------------------------------

// AND(a,b): 00→0, 01→0, 10→0, 11→1  mask = 0b1000 = 8
TEST(TableNodeTest, AndGateLogic) {
  auto andTable = makeMaskTable(/*size=*/2, /*mask=*/0b1000);
  TableNode andNode(andTable);
  andNode.addChild(std::make_unique<InputNode>(0));
  andNode.addChild(std::make_unique<InputNode>(1));

  EXPECT_FALSE(andNode.eval({false, false}));
  EXPECT_FALSE(andNode.eval({false, true }));
  EXPECT_FALSE(andNode.eval({true,  false}));
  EXPECT_TRUE (andNode.eval({true,  true }));
}

// NOT(x): 0→1, 1→0  mask = 0b01 = 1
TEST(TableNodeTest, NotGateLogic) {
  auto notTable = makeMaskTable(/*size=*/1, /*mask=*/0b01);
  TableNode notNode(notTable);
  notNode.addChild(std::make_unique<InputNode>(0));

  EXPECT_TRUE (notNode.eval({false}));
  EXPECT_FALSE(notNode.eval({true }));
}

// wiring‐size mismatch → bits.size() != (1<<children.size()) → exception
TEST(TableNodeTest, ThrowsOnTableSizeMismatch) {
  // table size=1 ⇒ bitsSize=2, but we wire 2 children ⇒ expect logic_error
  auto tinyTable = makeMaskTable(/*size=*/1, /*mask=*/0b01);
  TableNode node(tinyTable);
  node.addChild(std::make_unique<InputNode>(0));
  node.addChild(std::make_unique<InputNode>(1));

  EXPECT_THROW(node.eval({false, true}), std::logic_error);
}

//------------------------------------------------------------------------------
// SNLTruthTableTree tests
//------------------------------------------------------------------------------

// Compose (AND → NOT) to get NAND
TEST(SNLTruthTableTreeTest, ComposeAndNotIsNand) {
  // AND subtree
  auto andTbl = makeMaskTable(2, /*0b1000*/8);
  auto andN   = std::make_unique<TableNode>(andTbl);
  andN->addChild(std::make_unique<InputNode>(0));
  andN->addChild(std::make_unique<InputNode>(1));

  // NOT wrapper
  auto notTbl = makeMaskTable(1, /*0b01*/1);
  auto notN   = std::make_unique<TableNode>(notTbl);
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
  // identity (buf) table: size=1, mask=0b10, but we only ever need size=1 tests
  auto bufTbl = makeMaskTable(1, /*0b10*/2);
  auto bufN   = std::make_unique<TableNode>(bufTbl);
  bufN->addChild(std::make_unique<InputNode>(0));

  SNLTruthTableTree tree(std::move(bufN), /*numExternalInputs=*/1);

  EXPECT_THROW(tree.eval({}),       std::invalid_argument);
  EXPECT_THROW(tree.eval({true, false}), std::invalid_argument);
}

//------------------------------------------------------------------------------
// Dynamic child-addition tests
//------------------------------------------------------------------------------

TEST(TableNodeDynamicTest, ThreeInputOrLogic) {
  // 3-input OR: size=3, mask = 0b11111110 (any input=1 → output=1)
  auto or3Table = makeMaskTable(/*size=*/3, /*mask=*/0b11111110);
  TableNode or3Node(or3Table);

  // dynamically add all three leaves
  or3Node.addChild(std::make_unique<InputNode>(0));
  or3Node.addChild(std::make_unique<InputNode>(1));
  or3Node.addChild(std::make_unique<InputNode>(2));

  // exhaustively verify
  for (uint32_t i = 0; i < (1u << 3); ++i) {
    bool a = (i >> 0) & 1;
    bool b = (i >> 1) & 1;
    bool c = (i >> 2) & 1;
    bool expected = a || b || c;
    EXPECT_EQ(or3Node.eval({a, b, c}), expected)
        << "failed OR3 for bits=" << std::bitset<3>(i);
  }
}


TEST(TableNodeDynamicTest, TwoOfThreeThresholdLogic) {
  // threshold-2-of-3: size=3, mask = 0b11101000 (≥2 ones → output=1)
  auto thrTable = makeMaskTable(/*size=*/3, /*mask=*/0b11101000);
  TableNode thrNode(thrTable);

  // dynamically add three inputs
  thrNode.addChild(std::make_unique<InputNode>(0));
  thrNode.addChild(std::make_unique<InputNode>(1));
  thrNode.addChild(std::make_unique<InputNode>(2));

  // exhaustively verify
  for (uint32_t i = 0; i < (1u << 3); ++i) {
    int count = ((i>>0)&1) + ((i>>1)&1) + ((i>>2)&1);
    bool expected = (count >= 2);
    bool a = (i >> 0) & 1;
    bool b = (i >> 1) & 1;
    bool c = (i >> 2) & 1;
    EXPECT_EQ(thrNode.eval({a, b, c}), expected)
        << "failed threshold2/3 for bits=" << std::bitset<3>(i);
  }
}

//------------------------------------------------------------------------------
// Pyramid‐of‐And‐Gates test (8→4→2→1)
//------------------------------------------------------------------------------

TEST(TableNodePyramidTest, EightInputAndPyramid) {
    // Re-use the 2-input AND truth-table: mask = 0b1000 
    auto andTbl = makeMaskTable(/*size=*/2, /*mask=*/0b1000);

    // Level1: 4 AND nodes
    auto and01 = std::make_unique<TableNode>(andTbl);
    and01->addChild(std::make_unique<InputNode>(0));
    and01->addChild(std::make_unique<InputNode>(1));

    auto and23 = std::make_unique<TableNode>(andTbl);
    and23->addChild(std::make_unique<InputNode>(2));
    and23->addChild(std::make_unique<InputNode>(3));

    auto and45 = std::make_unique<TableNode>(andTbl);
    and45->addChild(std::make_unique<InputNode>(4));
    and45->addChild(std::make_unique<InputNode>(5));

    auto and67 = std::make_unique<TableNode>(andTbl);
    and67->addChild(std::make_unique<InputNode>(6));
    and67->addChild(std::make_unique<InputNode>(7));

    // Level2: two ANDs of Level1 results
    auto and0123 = std::make_unique<TableNode>(andTbl);
    and0123->addChild(std::move(and01));
    and0123->addChild(std::move(and23));

    auto and4567 = std::make_unique<TableNode>(andTbl);
    and4567->addChild(std::move(and45));
    and4567->addChild(std::move(and67));

    // Level3: root AND of the two Level2 nodes
    auto root = std::make_unique<TableNode>(andTbl);
    root->addChild(std::move(and0123));
    root->addChild(std::move(and4567));

    // Build a tree expecting 8 external inputs
    SNLTruthTableTree tree(std::move(root), /*numExternalInputs=*/8);

    // Test all 256 input patterns: only all-ones → true, else false
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

#include "SNLTruthTableTree.h"
#include "SNLTruthTable.h"
#include <gtest/gtest.h>
#include <bitset>
#include <bitset>

using namespace naja::NL;

// Build a tree that implements:
//   f(a0…a7) = NOR( OR(AND(a0,a4), AND(a2,a6)),
//                  OR(AND(a1,a5), AND(a3,a7)) )
//
// That is: 4 AND gates → 2 OR gates → final NOR.
// We construct it entirely via concat(), then verify on all 256 inputs.

TEST(SNLTruthTableTreeTest, Pyramid4Ands2Ors1NorViaConcat) {
    // 2-input tables: AND, OR, NOR
    auto andTbl = SNLTruthTable(2, 0b1000);  // only bit index 0b11 = 1
    auto  orTbl = SNLTruthTable(2, 0b1110);  // bits {01,10,11} = 1
    auto norTbl = SNLTruthTable(2, 0b0001);  // only bit index 0b00 = 1

    // Start from one external input A
    SNLTruthTableTree tree(
        std::make_unique<InputNode>(0),
        /*numExternalInputs=*/1
    );
    EXPECT_EQ(tree.size(), 1u);

    // Step 1: wrap A in a 2-input NOR → now expects A,B
    tree.concat(/*borderIndex=*/0, norTbl);
    EXPECT_EQ(tree.size(), 2u);

    // Step 2: graft left OR onto B-order slot 0 (the NOR’s first input)
    tree.concat(/*borderIndex=*/0, orTbl);
    EXPECT_EQ(tree.size(), 3u);

    // Step 3: graft right OR onto slot 1 (the NOR’s second input)
    tree.concat(/*borderIndex=*/1, orTbl);
    EXPECT_EQ(tree.size(), 4u);

    // Now leaves are 4 free inputs: indices 0,1,2,3
    // Step 4–7: graft the four AND gates at slots 0–3
    for (size_t slot = 0; slot < 4; ++slot) {
        tree.concat(slot, andTbl);
    }
    EXPECT_EQ(tree.size(), 8u);

    // Verify truth table on all 256 patterns
    for (uint32_t mask = 0; mask < (1u << 8); ++mask) {
        std::vector<bool> in(8);
        for (int i = 0; i < 8; ++i) {
            in[i] = ((mask >> i) & 1) != 0;
        }

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

// Add this to your SNLTruthTableTree_test.cpp

TEST(SNLTruthTableTreeTest, ConcatFullInvertsBothLeaves) {
    // Build a 2-input OR gate: f(a,b) = a OR b
    auto orTbl = makeMaskTable(2, /*mask=*/0b1110);
    auto root  = std::make_unique<TableNode>(orTbl);
    root->addChild(std::make_unique<InputNode>(0));
    root->addChild(std::make_unique<InputNode>(1));
    SNLTruthTableTree tree(std::move(root), /*numExternalInputs=*/2);

    // Prepare a NOT table to invert each leaf
    auto notTbl = makeMaskTable(1, /*mask=*/0b01);

    // Graft a NOT onto each of the two border leaves
    tree.concatFull({ notTbl, notTbl });

    // OR(not(a), not(b)) == NAND(a,b)
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
