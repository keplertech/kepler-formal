// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "Tree2BoolExpr.h"
#include "BoolExpr.h"
#include "DNL.h"
#include "NLDB0.h"
#include "SNLBusTerm.h"
#include "SNLTruthTable.h"
#include "SNLTruthTableTree.h"
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb_allocator.h>
#include <bitset>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace naja::NL;
using namespace KEPLER_FORMAL;

// Global concurrent map from DNL iso IDs to already-created BoolExpr pointers.
// This allows sharing of BoolExpr objects across conversions and threads.
tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*> Tree2BoolExpr::iso2boolExpr_ =
    tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*>();

// Helper typedefs for thread-local containers. Each pair stores a vector
// allocated with TBB allocator and a size_t representing the logical size.
// These thread-local containers avoid repeated allocations during traversal.

typedef std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>> TermsPair;

// Thread-local storage for DNF terms while building expressions for a node.
thread_local TermsPair termsETS;

TermsPair& getTErmsETS() {
  return termsETS;
}

size_t sizeOfTermsETS() {
  return getTErmsETS().size();
}

void clearTermsETS() {
  auto& termsLocal = getTErmsETS();
  termsLocal.clear();
}

void pushBackTermsETS(BoolExpr* term) {
  getTErmsETS().emplace_back(term);
}

bool emptyTermsETS() {
  return getTErmsETS().empty();
}

// Relevant inputs bitset per node: which child inputs actually affect the table.
// Stored as a vector<uint8_t> (bool-like) with a logical size.
typedef std::pair<std::vector<uint8_t, tbb::tbb_allocator<uint8_t>>, size_t> RelevantPair;

thread_local RelevantPair relevantETS;

RelevantPair& getRelevantETS() {
  return relevantETS;
}

// LCOV_EXCL_START
size_t sizeOfRelevantETS() {  // LCOV_EXCL_LINE
  return getRelevantETS().second;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

void clearRelevantETS() {
  auto& relevantLocal = getRelevantETS();
  relevantLocal.first.clear();
}

void setRelevantETS(size_t i, bool b) {
  auto& relevantLocal = getRelevantETS();
  if (i >= relevantLocal.second) {
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    assert(false && "setRelevantETS: index out of range");
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }
  relevantLocal.first[i] = b;
}

bool getRelevantETS(size_t i) {
  auto& relevantLocal = getRelevantETS();
  if (i >= relevantLocal.second) {
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    throw std::out_of_range("getRelevantETS: index out of range");
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }
  return relevantLocal.first[i];
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

// Ensure the relevant vector has at least n entries and initialize them to false.
// The logical size is stored in the second element of the pair.
void reserveRelevantETSwithFalse(size_t n) {
  auto& relevantLocal = getRelevantETS();
  auto& vec = relevantLocal.first;
  auto& sz = relevantLocal.second;
  if (vec.size() >= n) {
    // LCOV_EXCL_START
    vec.assign(n, false);  // LCOV_EXCL_LINE
    sz = n;  // LCOV_EXCL_LINE
    return;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  size_t oldSize = vec.size();
  vec.resize(n, false);
  vec.assign(n, false);
  sz = n;
}

// Memoization table (thread-local) mapping node IDs to BoolExpr* results.
// The vector is preallocated and indexed by node ID for O(1) lookup.
typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> MemoPair;

thread_local MemoPair memoETS;

MemoPair& getMemoETS() {
  return memoETS;
}

void clearMemoETS() {
  auto& memoLocal = getMemoETS();
  memoLocal.first.clear();
}

// Reserve and initialize the memo table to a given logical size.
// The vector elements are set to nullptr to indicate "not computed".
void reserveMemoETS(size_t n) {
  auto& memoLocal = getMemoETS();
  auto& vec = memoLocal.first;
  auto& sz = memoLocal.second;
  if (vec.size() >= n) {
    // LCOV_EXCL_START
    sz = n;  // LCOV_EXCL_LINE
    vec.assign(n, nullptr);  // LCOV_EXCL_LINE
    return;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  vec.resize(n);
  sz = n;
  vec.assign(n, nullptr);
}

void setMemoETS(size_t i, BoolExpr* expr) {
  auto& memoLocal = getMemoETS();
  assert(i < memoLocal.second && "setMemoETS: index out of range");
  memoLocal.first[i] = expr;
}

BoolExpr* getMemoETS(size_t i) {
  auto& memoLocal = getMemoETS();
  assert(i < memoLocal.second && "getMemoETS: index out of range");
  return memoLocal.first[i];
}

// Temporary storage for child BoolExpr pointers while processing a table node.
typedef std::pair<std::vector<BoolExpr*, tbb::tbb_allocator<BoolExpr*>>, size_t> ChildFETSPair;

thread_local ChildFETSPair childFETS;

ChildFETSPair& getChildFETS() {
  return childFETS;
}

void clearChildFETS() {
  auto& childLocal = getChildFETS();
  childLocal.second = 0;
}

// Reserve child storage and initialize entries to nullptr.
void reserveChildFETS(size_t n) {
  auto& childLocal = getChildFETS();
  auto& vec = childLocal.first;
  auto& sz = childLocal.second;
  if (vec.size() >= n) {
    sz = n;
    vec.assign(n, nullptr);
    return;
  }
  vec.resize(n);
  sz = n;
  vec.assign(n, nullptr);
}

BoolExpr* getChildFETS(size_t i) {
  auto& childLocal = getChildFETS();
  assert(i < childLocal.second && "getChildFETS: index out of range");
  return childLocal.first[i];
}

void setChildFETS(size_t i, BoolExpr* expr) {
  auto& childLocal = getChildFETS();
  assert(i < childLocal.second && "setChildFETS: index out of range");
  childLocal.first[i] = expr;
}

BoolExpr* buildTableSelectTruthTableExprRecursive(
    uint32_t addressSize,
    uint32_t depth,
    uint32_t addressIndex,
    uint64_t prefix) {
  const uint32_t remaining = addressSize - addressIndex;
  if (remaining < 63) {
    if ((prefix << remaining) >= depth) {
      return BoolExpr::createFalse();
    }
  } else if (prefix != 0) {
    return BoolExpr::createFalse();
  }

  if (addressIndex == addressSize) {
    assert(prefix < depth);
    return getChildFETS(addressSize + static_cast<uint32_t>(prefix));
  }

  BoolExpr* low = buildTableSelectTruthTableExprRecursive(
      addressSize, depth, addressIndex + 1, prefix << 1);
  BoolExpr* high = buildTableSelectTruthTableExprRecursive(
      addressSize, depth, addressIndex + 1, (prefix << 1) | 1);
  if (low == high) {
    return low;
  }

  BoolExpr* addressBit = getChildFETS(addressIndex);
  return BoolExpr::Or(
      BoolExpr::And(BoolExpr::Not(addressBit), low),
      BoolExpr::And(addressBit, high));
}

BoolExpr* buildTableSelectTruthTableExpr(const SNLTruthTable& tbl,
                                         uint32_t k) {
  const uint32_t addressSize = tbl.getTableSelectAddressSize();
  const uint32_t depth = tbl.getTableSelectDepth();
  if (k != addressSize + depth) {
    throw std::runtime_error("TABLE_SELECT truth table arity mismatch");
  }
  return buildTableSelectTruthTableExprRecursive(
      addressSize, depth, 0, 0);
}  // LCOV_EXCL_LINE

namespace {

using BoolExprBits = std::vector<BoolExpr*>;

struct DivModExprs {
  BoolExprBits quotient;
  BoolExprBits remainder;
};

struct SubtractExprs {
  BoolExprBits difference;
  BoolExpr* borrow = nullptr;
};

std::mutex divModBuildMutex;

BoolExpr* makeIte(
    BoolExpr* condition,
    BoolExpr* whenTrue,
    BoolExpr* whenFalse) {
  if (whenTrue == whenFalse) {
    return whenTrue;
  }
  return BoolExpr::Or(
      BoolExpr::And(condition, whenTrue),
      BoolExpr::And(BoolExpr::Not(condition), whenFalse));
}

SubtractExprs subtractBits(
    const BoolExprBits& lhs,
    const BoolExprBits& rhs) {
  if (lhs.size() != rhs.size()) {
    throw std::runtime_error("DIVMOD subtract width mismatch");  // LCOV_EXCL_LINE
  }

  SubtractExprs result;
  result.difference.resize(lhs.size());
  BoolExpr* borrow = BoolExpr::createFalse();
  for (size_t bit = 0; bit < lhs.size(); ++bit) {
    const auto lhsXorRhs = BoolExpr::Xor(lhs[bit], rhs[bit]);
    result.difference[bit] = BoolExpr::Xor(lhsXorRhs, borrow);
    borrow = BoolExpr::Or(
        BoolExpr::And(
            BoolExpr::Not(lhs[bit]),
            BoolExpr::Or(rhs[bit], borrow)),
        BoolExpr::And(rhs[bit], borrow));
  }
  result.borrow = borrow;
  return result;
}

BoolExprBits conditionalNegate(
    const BoolExprBits& bits,
    BoolExpr* negate) {
  BoolExprBits result(bits.size());
  BoolExpr* carry = negate;
  for (size_t bit = 0; bit < bits.size(); ++bit) {
    BoolExpr* inverted = BoolExpr::Xor(bits[bit], negate);
    result[bit] = BoolExpr::Xor(inverted, carry);
    carry = BoolExpr::And(inverted, carry);
  }
  return result;
}

DivModExprs buildUnsignedDivMod(
    const BoolExprBits& dividend,
    const BoolExprBits& divisor) {
  if (dividend.empty() || dividend.size() != divisor.size()) {
    throw std::runtime_error("DIVMOD operand width mismatch");  // LCOV_EXCL_LINE
  }

  const size_t width = dividend.size();
  DivModExprs result;
  result.quotient.assign(width, BoolExpr::createFalse());

  BoolExprBits remainder(width + 1, BoolExpr::createFalse());
  BoolExprBits extendedDivisor = divisor;
  extendedDivisor.push_back(BoolExpr::createFalse());

  for (size_t dividendBit = width; dividendBit-- > 0;) {
    BoolExprBits shifted(width + 1);
    shifted[0] = dividend[dividendBit];
    for (size_t bit = 1; bit <= width; ++bit) {
      shifted[bit] = remainder[bit - 1];
    }

    const auto subtraction = subtractBits(shifted, extendedDivisor);
    BoolExpr* takeDifference = BoolExpr::Not(subtraction.borrow);
    result.quotient[dividendBit] = takeDifference;
    for (size_t bit = 0; bit <= width; ++bit) {
      remainder[bit] = makeIte(
          takeDifference, subtraction.difference[bit], shifted[bit]);
    }
  }

  result.remainder.assign(remainder.begin(), remainder.begin() + width);
  return result;
}

DivModExprs buildDivMod(
    const BoolExprBits& dividend,
    const BoolExprBits& divisor,
    bool isSigned) {
  if (!isSigned) {
    return buildUnsignedDivMod(dividend, divisor);
  }

  BoolExpr* dividendSign = dividend.back();
  BoolExpr* divisorSign = divisor.back();
  auto magnitudeResult = buildUnsignedDivMod(
      conditionalNegate(dividend, dividendSign),
      conditionalNegate(divisor, divisorSign));
  magnitudeResult.quotient = conditionalNegate(
      magnitudeResult.quotient,
      BoolExpr::Xor(dividendSign, divisorSign));
  magnitudeResult.remainder = conditionalNegate(
      magnitudeResult.remainder, dividendSign);
  return magnitudeResult;
}

BoolExpr* buildDivModTruthTableExpr(
    const SNLTruthTable& table,
    uint32_t arity,
    const SNLTruthTableTree::Node* node,
    naja::DNL::DNLID isoID) {
  const uint32_t width = table.getDivModWidth();
  if (arity != width * 2 || node == nullptr) {
    throw std::runtime_error("DIVMOD truth table metadata mismatch");  // LCOV_EXCL_LINE
  }

  std::lock_guard<std::mutex> lock(divModBuildMutex);
  if (isoID != naja::DNL::DNLID_MAX) {
    const auto cached = Tree2BoolExpr::iso2boolExpr_.find(isoID);
    if (cached != Tree2BoolExpr::iso2boolExpr_.end() &&
        cached->second != nullptr && cached->second->isValid()) {
      return cached->second;
    }
  }

  // Naja's packed [WIDTH-1:0] operand terms are flattened MSB first.
  // The divider helpers use conventional LSB-first vectors.
  BoolExprBits dividend(width);
  BoolExprBits divisor(width);
  for (uint32_t bit = 0; bit < width; ++bit) {
    dividend[bit] = getChildFETS(width - 1 - bit);
    divisor[bit] = getChildFETS(width * 2 - 1 - bit);
  }

  const auto expressions =
      buildDivMod(dividend, divisor, table.isDivModSigned());
  const auto& currentTerm =
      naja::DNL::get()->getDNLTerminalFromID(node->data.termid);
  const auto& instance = currentTerm.getDNLInstance();
  const auto* model = instance.getSNLModel();
  if (!NLDB0::isDivMod(model)) {
    throw std::runtime_error("DIVMOD truth table used by a non-divmod model");  // LCOV_EXCL_LINE
  }

  const auto* quotientTerm = NLDB0::getDivModQuotient(model);
  const auto* remainderTerm = NLDB0::getDivModRemainder(model);
  for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX &&
       termID <= instance.getTermIndexes().second;
       ++termID) {
    const auto& term = naja::DNL::get()->getDNLTerminalFromID(termID);
    const auto* bitTerm = term.getSnlBitTerm();
    if (bitTerm->getDirection() != SNLBitTerm::Direction::Output) {
      continue;
    }

    const auto bit = static_cast<size_t>(bitTerm->getBit());
    if (bit >= width) {
      throw std::runtime_error("DIVMOD output bit is out of range");  // LCOV_EXCL_LINE
    }
    BoolExpr* expression = nullptr;
    if (bitTerm->getID() == quotientTerm->getID()) {
      expression = expressions.quotient[bit];
    } else if (bitTerm->getID() == remainderTerm->getID()) {
      expression = expressions.remainder[bit];
    } else {
      continue;  // LCOV_EXCL_LINE
    }

    const auto outputIsoID = term.getIsoID();
    if (outputIsoID != naja::DNL::DNLID_MAX) {
      Tree2BoolExpr::iso2boolExpr_.insert({outputIsoID, expression});
    }
  }

  const auto outputBit = table.getDivModOutputBit();
  BoolExpr* selected =
      table.getDivModResult() == SNLTruthTable::DivModResult::QUOTIENT
          ? expressions.quotient[outputBit]
          : expressions.remainder[outputBit];
  if (isoID != naja::DNL::DNLID_MAX) {
    const auto cached = Tree2BoolExpr::iso2boolExpr_.find(isoID);
    if (cached != Tree2BoolExpr::iso2boolExpr_.end()) {
      selected = cached->second;
    }
  }
  return selected;
}

}  // namespace

BoolExpr* buildGenericTruthTableExpr(
    const SNLTruthTable& tbl,
    uint32_t k,
    const SNLTruthTableTree::Node* node,
    naja::DNL::DNLID isoID) {
  assert(k > 0);

  BoolExpr* expr = getChildFETS(0);
  switch (tbl.getGenericType()) {
    case SNLTruthTable::GenericType::AND:
    case SNLTruthTable::GenericType::NAND:
      for (uint32_t j = 1; j < k; ++j) {
        expr = BoolExpr::And(expr, getChildFETS(j));
      }
      if (tbl.getGenericType() == SNLTruthTable::GenericType::NAND) {
        expr = BoolExpr::Not(expr);
      }
      return expr;
    case SNLTruthTable::GenericType::OR:
    case SNLTruthTable::GenericType::NOR:
      for (uint32_t j = 1; j < k; ++j) {
        expr = BoolExpr::Or(expr, getChildFETS(j));
      }
      if (tbl.getGenericType() == SNLTruthTable::GenericType::NOR) {
        expr = BoolExpr::Not(expr);
      }
      return expr;
    case SNLTruthTable::GenericType::XOR:
    case SNLTruthTable::GenericType::XNOR:
      for (uint32_t j = 1; j < k; ++j) {
        expr = BoolExpr::Xor(expr, getChildFETS(j));
      }
      if (tbl.getGenericType() == SNLTruthTable::GenericType::XNOR) {
        expr = BoolExpr::Not(expr);
      }
      return expr;
    case SNLTruthTable::GenericType::TABLE_SELECT:
      return buildTableSelectTruthTableExpr(tbl, k);
    case SNLTruthTable::GenericType::DIVMOD:
      return buildDivModTruthTableExpr(tbl, k, node, isoID);
    case SNLTruthTable::GenericType::NONE:
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      break;
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  // LCOV_DISABLED_START
  throw std::runtime_error("Unsupported generic truth table type");
  // LCOV_DISABLED_STOP
  // LCOV_EXCL_STOP
}

BoolExpr* buildGenericTruthTableExpr(const SNLTruthTable& tbl, uint32_t k) {
  return buildGenericTruthTableExpr(
      tbl, k, nullptr, naja::DNL::DNLID_MAX);
}

// Frame type used for explicit stack-based post-order traversal.
// Each frame holds a pointer to a node and a boolean indicating whether
// the node has been visited (post-visit) or not (pre-visit).
using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
thread_local std::vector<Frame, tbb::tbb_allocator<Frame>> stackETS;

std::vector<Frame, tbb::tbb_allocator<Frame>>& getStackETS() {
  return stackETS;
}

// Main conversion routine: converts an SNLTruthTableTree into a BoolExpr.
// The varNames vector maps SNL variable indices to desired variable IDs
// (or special markers like 0/1 for constants).
BoolExpr* Tree2BoolExpr::convert(
  const SNLTruthTableTree& tree, const std::vector<size_t>& varNames) {

  const auto root = tree.getRoot();
  if (!root) {
    // LCOV_EXCL_START
    return nullptr; // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // Determine maximum node ID to size memoization structures.
  size_t maxID = tree.getMaxID();

  // 2) memo table: clear and reserve memoization storage for all node IDs.
  clearMemoETS();
  reserveMemoETS(maxID + 1);

  // 3) post-order build using an explicit stack to avoid recursion.
  auto & stack = getStackETS();
  stack.clear();
  stack.emplace_back(root.get(), false);

  while (!stack.empty()) {
    Frame f = stack.back();
    stack.pop_back();
    const SNLTruthTableTree::Node* node = f.first;

    naja::DNL::DNLID isoID = naja::DNL::DNLID_MAX;
    if (node->type != SNLTruthTableTree::Node::Type::Input) {
      isoID = naja::DNL::get()->getDNLTerminalFromID(node->data.termid).getIsoID();
    }

    bool visited = f.second;
    size_t id = node->nodeID;

    if (!visited) {
      // Pre-visit: attempt to reuse an existing BoolExpr from iso2boolExpr_
      // if the node corresponds to a DNL terminal that was already converted.
      if (isoID != naja::DNL::DNLID_MAX) {
        auto it = iso2boolExpr_.find(isoID);
        if (it != iso2boolExpr_.end() && it->second != nullptr &&
            it->second->isValid()) {
          setMemoETS(id, it->second);
        }
      }
      // If memo already contains an expression for this node, skip processing.
      if (getMemoETS(id) != nullptr) {
        continue;
      }

      // If node is a Table or P node, push it back as visited and push children.
      if (node->type == SNLTruthTableTree::Node::Type::Table || node->type == SNLTruthTableTree::Node::Type::P) {
        stack.emplace_back(node, true);
        for (const auto& c : node->childrenIds) stack.emplace_back(node->tree->nodeFromId(c).get(), false);
      } else {
        // Input node handling: map input nodes to variables or constants.
        assert(node->type == SNLTruthTableTree::Node::Type::Input);
        if (node->parentIds.size() > 1) {
          #ifdef DEBUG_PRINTS
          // Debug logging for inputs with multiple parents (should be rare).
          for (const auto& pid : node->parentIds) {
            DEBUG_LOG("%s\n", naja::DNL::get()->getDNLTerminalFromID(tree.nodeFromId(pid)->data.termid)
                     .getSnlBitTerm()->getString().c_str());
            DEBUG_LOG("of model %s\n", naja::DNL::get()->getDNLTerminalFromID(tree.nodeFromId(pid)->data.termid)
                   .getDNLInstance().getSNLModel()->getString().c_str());
          }
          #endif
        // LCOV_EXCL_START
        }  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        if (node->parentIds.empty()) { 
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Input node has no parent"); 
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        assert(node->parentIds.size() == 1);
        SNLTruthTableTree::Node* const parent = node->tree->nodeFromId(node->parentIds[0]).get();
        assert(parent && parent->type == SNLTruthTableTree::Node::Type::P);
        if (isoID != naja::DNL::DNLID_MAX) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error(
          // LCOV_DISABLED_STOP
              "Input node unexpectedly has a cacheable iso ID");
          // LCOV_EXCL_STOP
        }
        assert(parent->data.termid < varNames.size());
        const auto& name = varNames[parent->data.termid];
        if (name == (size_t)-1) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Input variable index is SIZE_MAX");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        // Special handling for constant mappings: 0 -> false, 1 -> true.
        if (name == 0) {
           BoolExpr* expr = BoolExpr::createFalse();
           setMemoETS(id, expr);
        } else if (name == 1) {
           BoolExpr* expr = BoolExpr::createTrue();
           setMemoETS(id, expr);
        } else {
          // Normal variable mapping.
          BoolExpr* expr = BoolExpr::Var(name);
          setMemoETS(id, expr);
        }
      }
    } else {
      // Post-visit: node is a Table or P node and all children have been processed.
      const SNLTruthTable& tbl = node->getTruthTable();
      DEBUG_LOG("Processing node ID %zu with table:\n%s\n", id, tbl.toString().c_str());
      uint32_t k = tbl.size();
      DEBUG_LOG("Node ID %zu has %u inputs\n", id, k);

      if (tbl.all0() || tbl.all1()) {
        BoolExpr* expr =
            tbl.all1() ? BoolExpr::createTrue() : BoolExpr::createFalse();
        if (isoID != naja::DNL::DNLID_MAX) {
          auto result = iso2boolExpr_.insert({isoID, expr});
          if (!result.second) {
            // LCOV_EXCL_START
            expr = result.first->second;  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        setMemoETS(id, expr);
        continue;
      }

      {
        // Gather child BoolExpr pointers into a temporary array for quick access.
        clearChildFETS();
        reserveChildFETS(k);
        for (uint32_t i = 0; i < k; ++i) {
          size_t cid = node->tree->nodeFromId(node->childrenIds[i])->nodeID;
          setChildFETS(i, getMemoETS(cid));
        }

        if (tbl.isGeneric()) {
          BoolExpr* expr = buildGenericTruthTableExpr(tbl, k, node, isoID);
          if (isoID != naja::DNL::DNLID_MAX) {
            // LCOV_EXCL_START
            // The duplicate insert path only happens through concurrent reuse of
            // the same isoID during conversion.
            auto result = iso2boolExpr_.insert({isoID, expr});
            if (!result.second) {
              // LCOV_DISABLED_START
              expr = result.first->second;
            }
            // LCOV_DISABLED_STOP
            // LCOV_EXCL_STOP
          }
          setMemoETS(id, expr);
          continue;
        }

        const uint64_t rows = uint64_t{1} << k;
        DEBUG_LOG("Node ID %zu has %llu rows\n", id, rows);

        // Determine which inputs actually matter for this truth table.
        // For each input j, check if flipping bit j changes the table output.
        clearRelevantETS();
        reserveRelevantETSwithFalse(k);
        for (uint32_t j = 0; j < k; ++j) {
          for (uint64_t m = 0; m < rows; ++m) {
            bool b0 = tbl.bits().bit(m);
            bool b1 = tbl.bits().bit(m ^ (uint64_t{1} << j));
            if (b0 != b1) { setRelevantETS(j, true); break; }
          }
        }

        // Count how many inputs are relevant.
        size_t numRelIdx = 0;
        for (uint32_t j = 0; j < k; ++j) {
          if (getRelevantETS(j)) {
            numRelIdx++;
          }
        }

        // The algorithm expects at least one relevant input for a PI node.
        assert(numRelIdx > 0 && "No relevant inputs for node");
        {
          // Build DNF terms by iterating over rows where the table output is 1.
          // For each such row, create a conjunction of literals for relevant inputs.
          clearTermsETS();
          for (uint64_t m = 0; m < rows; ++m) {
            if (!tbl.bits().bit(m)) {
              continue;
            }
            BoolExpr* term = nullptr;
            bool firstLit = true;
            BoolExpr* lit = nullptr;
            // For each relevant input, pick the literal (child or its negation)
            // according to the bit value in row m.
            for (uint32_t j = 0; j < k; ++j) { 
              if (!getRelevantETS(j)) {
                // LCOV_EXCL_START
                continue;  // LCOV_EXCL_LINE
                // LCOV_EXCL_STOP
              }
              bool bit1 = ((m >> j) & 1) != 0;
              lit = bit1 ? getChildFETS(j) : BoolExpr::Not(getChildFETS(j));
              if (firstLit) { term = lit; firstLit = false; }
              else { assert(term != nullptr); assert(lit != nullptr); term = BoolExpr::And(term, lit); }
            }
            // Only push the term if at least one literal was included.
            if (term) { pushBackTermsETS(std::move(term)); }
          }

          // There must be at least one term for a PI node.
          assert(!emptyTermsETS()); // Should be a PI

          {
            // Fold the list of terms into a single expression by OR-ing them.
            BoolExpr* expr = getTErmsETS()[0];
            DEBUG_LOG("number of rows for node ID %zu: %zu\n", id, sizeOfTermsETS());
            for (size_t t = 1; t < sizeOfTermsETS(); ++t) {
              expr = BoolExpr::Or(expr, getTErmsETS()[t]);
              DEBUG_LOG("Intermediate OR expr for node ID %zu: %s\n", id, expr->toString().c_str());
            }
            // Store the resulting expression in the memo table and in the iso map.
            if (isoID != naja::DNL::DNLID_MAX) {
                auto result = iso2boolExpr_.insert({isoID, expr});
                if (!result.second) {
                    // Another thread inserted concurrently.
                    // Reuse canonical instance (do NOT delete expr; ownership may not be raw).
                    expr = result.first->second;
                }
            }
            setMemoETS(id, expr);
            DEBUG_LOG("Bool expression for node ID %zu: %s\n", id, expr->toString().c_str());
          }
        }
      }
    }
  }

  // 4) return root expression from memo table.
  return getMemoETS(root->nodeID);
}
