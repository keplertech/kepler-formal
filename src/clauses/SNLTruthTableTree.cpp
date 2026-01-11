// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "SNLTruthTableTree.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdio>
#include <limits>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <tbb/enumerable_thread_specific.h>
#include "tbb/concurrent_vector.h"


using namespace KEPLER_FORMAL;

// #define DEBUG_CHECKS

// #define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

typedef std::pair<std::vector<KEPLER_FORMAL::SNLTruthTableTree::BorderLeaf, tbb::tbb_allocator<KEPLER_FORMAL::SNLTruthTableTree::BorderLeaf>>, size_t> BorderLeavesPair;

// ets for std::vector<BorderLeaf, tbb::tbb_allocator<BorderLeaf>> newBorderLeaves;
// tbb::enumerable_thread_specific<BorderLeavesPair> newBorderLeavesETS;

// tbb::concurrent_vector<BorderLeavesPair*>
//     newBorderLeavesETSVector =
//         tbb::concurrent_vector<BorderLeavesPair*>(40, nullptr);

// void initNewBorderLeavesETS() {
//  size_t idx = tbb::this_task_arena::current_thread_index() >= 0 ? tbb::this_task_arena::current_thread_index() : 0;
//  if (newBorderLeavesETSVector.size() <= idx) {
//     for (size_t i = newBorderLeavesETSVector.size(); i <= idx; i++) {
//       newBorderLeavesETSVector.emplace_back(nullptr);
//     }
//   }
//   if (newBorderLeavesETSVector
//           [idx] == nullptr) {
//     newBorderLeavesETSVector[idx] =
//         &newBorderLeavesETS.local();
//   }
// }

thread_local BorderLeavesPair newBorderLeavesETS;

BorderLeavesPair& getNewBorderLeavesETS() {
  //size_t idx = tbb::this_task_arena::current_thread_index() >= 0 ? tbb::this_task_arena::current_thread_index() : 0;
  //return *newBorderLeavesETSVector[idx];
  return newBorderLeavesETS;
}

size_t getSizeOfNewBorderLeavesETS() {
  return getNewBorderLeavesETS().first.size();
}

void pushBackNewBorderLeavesETS(const KEPLER_FORMAL::SNLTruthTableTree::BorderLeaf& leaf) {
  // auto& pair = getNewBorderLeavesETS();
  // auto& vec = pair.first;
  // auto& sz = pair.second;
  // if (vec.size() > sz) {
  //   vec[sz] = leaf;
  //   sz++;
  //   return;
  // }
  // vec.emplace_back(leaf);
  // sz++;
  getNewBorderLeavesETS().first.emplace_back(leaf);
}

void clearNewBorderLeavesETS() {
  auto& pair = getNewBorderLeavesETS();
  pair.first.clear();
}

// same for
// std::vector<std::vector<std::shared_ptr<Node>, tbb::tbb_allocator<std::shared_ptr<Node>>>,
//     tbb::tbb_allocator<std::vector<std::shared_ptr<Node>, tbb::tbb_allocator<std::shared_ptr<Node>>>>> resolvedChildren(
//       nodes_.size());

// prefer clear aliases so allocator types are consistent
using NodeRaw = SNLTruthTableTree::Node*;
using NodeVecRaw = std::vector<NodeRaw, tbb::tbb_allocator<NodeRaw>>;
using NodeVecVecRaw = std::vector<NodeVecRaw, tbb::tbb_allocator<NodeVecRaw>>;
using ResolvedChildrenPairRaw = std::pair<NodeVecVecRaw, size_t>;
// tbb::enumerable_thread_specific<ResolvedChildrenPairRaw> resolvedChildrenETS;

// tbb::concurrent_vector<ResolvedChildrenPairRaw*>
//     resolvedChildrenETSVector =
//         tbb::concurrent_vector<ResolvedChildrenPairRaw*>(40, nullptr);

// void initResolvedChildrenETS() {
//   size_t idx = tbb::this_task_arena::current_thread_index() >= 0 ? tbb::this_task_arena::current_thread_index() : 0;
//   if (resolvedChildrenETSVector.size() <=
//       idx) {
//     for (size_t i = resolvedChildrenETSVector.size();
//          i <= idx; i++) {
//       resolvedChildrenETSVector.emplace_back(nullptr);
//     }
//   }
//   if (resolvedChildrenETSVector
//           [idx] == nullptr) {
//     resolvedChildrenETSVector[idx] =
//         &resolvedChildrenETS.local();
//   }
// }

thread_local ResolvedChildrenPairRaw resolvedChildrenETS;

ResolvedChildrenPairRaw& getResolvedChildrenETS() {
  //size_t idx = tbb::this_task_arena::current_thread_index() >= 0 ? tbb::this_task_arena::current_thread_index() : 0;
  //return *resolvedChildrenETSVector
  //    [idx];
  return resolvedChildrenETS;
}

void reserveResolvedChildrenETS(size_t n) {
  auto& pair = getResolvedChildrenETS();
  for (size_t i = 0; i < pair.first.size(); ++i) {
    pair.first[i].clear();
  }
  if (pair.first.size() < n) {
    for (size_t i = pair.first.size(); i < n; ++i) {
      pair.first.emplace_back(
          NodeVecRaw());
    }
  }
}

// Init Ptable holder
const SNLTruthTable SNLTruthTableTree::PtableHolder_ = SNLTruthTable(1, 2);

// diagnostic global
static std::atomic<size_t> g_live_nodes{0};

// NodeLifetimeCounter impl
// SNLTruthTableTree::Node::NodeLifetimeCounter::NodeLifetimeCounter()  {
// g_live_nodes.fetch_add(1, std::memory_order_relaxed); }
// SNLTruthTableTree::Node::NodeLifetimeCounter::~NodeLifetimeCounter() {
// g_live_nodes.fetch_sub(1, std::memory_order_relaxed); }

//----------------------------------------------------------------------
// Node ctors / dtor
//----------------------------------------------------------------------
SNLTruthTableTree::Node::Node(uint32_t idx, SNLTruthTableTree* t)
    : type(Type::Input),
      /*nodeID(0),*/ nodeID(SNLTruthTableTree::kInvalidId),
      tree(t) {
  data.inputIndex = idx;
  if (tree && tree->lastID_ == std::numeric_limits<unsigned>::max()) {
    // LCOV_EXCL_START
    throw std::overflow_error("Node ID overflow");
    // LCOV_EXCL_STOP
  }
  if (tree)
    nodeID = (uint32_t)tree->lastID_++;
}

SNLTruthTableTree::Node::Node(SNLTruthTableTree* t,
                              naja::DNL::DNLID instid,
                              naja::DNL::DNLID term,
                              Type type_)
    : type(type_),
      /*nodeID(0),*/ nodeID(SNLTruthTableTree::kInvalidId),
      tree(t) {
  data.termid = term;
  if (tree && tree->lastID_ == std::numeric_limits<unsigned>::max()) {
    // LCOV_EXCL_START
    throw std::overflow_error("Node ID overflow");
    // LCOV_EXCL_STOP
  }
  if (tree)
    nodeID = (uint32_t)tree->lastID_++;
  if (type == Type::Table) {
    truthTable = SNLDesignModeling::getTruthTable(naja::DNL::get()
                                             ->getDNLTerminalFromID(data.termid)
                                             .getDNLInstance()
                                             .getSNLModel(), naja::DNL::get()
                                    ->getDNLTerminalFromID(data.termid)
                                    .getSnlBitTerm()
                                    ->getOrderID());
  }
}

SNLTruthTableTree::Node::~Node() {
  childrenIds.clear();
  parentIds.clear();
  tree = nullptr;
}

//----------------------------------------------------------------------
// Node::getTruthTable
//----------------------------------------------------------------------
const SNLTruthTable& SNLTruthTableTree::Node::getTruthTable() const {
  if (type == Type::Table) {
    if (!truthTable.isInitialized()) {
      // LCOV_EXCL_START
      throw std::logic_error("getTruthTable: uninitialized Table node");
      // LCOV_EXCL_STOP
    }
    return truthTable;
  } else if (type == Type::P || type == Type::Input) {
    return PtableHolder_;
  }
  // LCOV_EXCL_START
  throw std::logic_error("getTruthTable: not a Table/P node");
  // LCOV_EXCL_STOP
}

static std::shared_ptr<SNLTruthTableTree::Node> nullNodePtr = nullptr;

//----------------------------------------------------------------------
// nodeFromId helper
//----------------------------------------------------------------------
const std::shared_ptr<SNLTruthTableTree::Node>& SNLTruthTableTree::nodeFromId(
    uint32_t id) const {
  assert((id != kInvalidId));
  assert((id >= kIdOffset));
  size_t idx = (size_t)(id - kIdOffset);
  assert(idx < nodes_.size());
  const auto& sp = nodes_[idx];
  assert(sp.get() != nullptr);
  assert(sp->nodeID == id);
  return sp;
}

//----------------------------------------------------------------------
// Node::eval (resolves children via ids)
//----------------------------------------------------------------------
bool SNLTruthTableTree::Node::eval(const std::vector<bool>& extInputs) const {
  if (type != Type::Table && type != Type::P && type != Type::Input) {
    // LCOV_EXCL_START
    throw std::logic_error("eval: node not Table/P/Input");
    // LCOV_EXCL_STOP
  }
  const auto& tbl = getTruthTable();
  auto arity = tbl.size();
  if (childrenIds.size() != arity) {
    // LCOV_EXCL_START
    throw std::logic_error("TableNode: children count mismatch");
    // LCOV_EXCL_STOP
  }

  uint32_t idx = 0;
  for (uint32_t i = 0; i < arity; ++i) {
    bool bit = false;
    uint32_t cid = childrenIds[i];
    if (cid == kInvalidId) {
      // LCOV_EXCL_START
      throw std::logic_error("Invalid child id");
      // LCOV_EXCL_STOP
    }
    Node* const childSp = tree->nodeFromId(cid).get();
    if (!childSp) {
      // LCOV_EXCL_START
      throw std::logic_error("Null child node");
      // LCOV_EXCL_STOP
    }
    if (childSp->type == Type::Input) {
      size_t inx = childSp->data.inputIndex;
      if (inx >= extInputs.size()) {
        // LCOV_EXCL_START
        throw std::out_of_range("Input index out of range");
        // LCOV_EXCL_STOP
      }
      bit = extInputs[inx];
    } else {
      bit = childSp->eval(extInputs);
    }
    if (bit)
      idx |= (1u << i);
  }
  return tbl.bits().bit(idx);
}

//----------------------------------------------------------------------
// addChildId: set parent/child relationship via ids
//----------------------------------------------------------------------
void SNLTruthTableTree::Node::addChildId(uint32_t childId) {
  if (childId == kInvalidId) {
    // LCOV_EXCL_START
    throw std::invalid_argument("addChildId: invalid id");
    // LCOV_EXCL_STOP
  }

#ifdef DEBUG_CHECKS
  uint32_t cur = this->parentId;
  while (cur != SNLTruthTableTree::kInvalidId) {
    if (cur == childId)
      throw std::invalid_argument("addChildId: cycle detected");
    auto p = tree->nodeFromId(cur);
    if (!p)
      break;
    cur = p->parentId;
  }
#endif

  childrenIds.emplace_back(childId);

  auto childSp = tree->nodeFromId(childId);
  if (childSp)
    childSp->parentIds.emplace_back(this->nodeID);
}

//----------------------------------------------------------------------
// allocateNode helper - assigns id before publishing into nodes_
//----------------------------------------------------------------------
uint32_t SNLTruthTableTree::allocateNode(std::shared_ptr<Node>& np) {
  if (!np) {
    // LCOV_EXCL_START
    throw std::invalid_argument("allocateNode: null");
    // LCOV_EXCL_STOP
  }
  auto iter = termid2nodeid_.find(np->data.termid);
  if (np->type == Node::Type::Table && iter != termid2nodeid_.end()) {
    np = nodeFromId(iter->second);
    return iter->second;
  }
  uint32_t id = static_cast<uint32_t>(nodes_.size()) + kIdOffset;
  np->nodeID = id;
  np->tree = this;
  nodes_.emplace_back(np);
  if (np->type == Node::Type::Table) {
    termid2nodeid_[np->data.termid] = id;
  }
  return id;
}

//----------------------------------------------------------------------
// updateBorderLeaves
//----------------------------------------------------------------------
void SNLTruthTableTree::updateBorderLeaves() {
  borderLeaves_.clear();
  size_t externalIndex = 0;
  if (rootId_ == kInvalidId)
    return;
  std::vector<uint32_t, tbb::tbb_allocator<uint32_t>> stk;
  stk.reserve(64);
  stk.emplace_back(rootId_);
  // replace:
  // std::set<uint32_t, std::less<uint32_t>, tbb::tbb_allocator<uint32_t>> visited;

  // with:
  // std::unordered_set<uint32_t,
  //                   std::hash<uint32_t>,
  //                   std::equal_to<uint32_t>,
  //                   tbb::tbb_allocator<uint32_t>> visited;
  // visited.reserve(nodes_.size() * 2);        // avoid rehashes; tune factor to expected size
  // visited.max_load_factor(0.7f);             // optional: control bucket density
  markAllUnvisited();
  while (!stk.empty()) {
    uint32_t nid = stk.back();
    stk.pop_back();
    //auto [itr, inserted] = visited.insert(nid);
    //if (!inserted)
    //  continue;
    const auto& nsp = nodeFromId(nid).get();
    if (nsp->visited)
      continue;
    nsp->visited = true;
    if (!nsp)
      assert(false && "updateBorderLeaves: null node in tree");
    assert(nsp->childrenIds.size() > 0);
    for (size_t i = 0; i < nsp->childrenIds.size(); ++i) {
      uint32_t cid = nsp->childrenIds[i];
      const auto& ch = nodeFromId(cid).get();
      if (!ch)
        assert(false && "updateBorderLeaves: null child node in tree");
      if (ch->type == Node::Type::Input || ch->type == Node::Type::P) {
        BorderLeaf bl;
        if (ch->type == Node::Type::P) {
          bl.parentId = cid;
          bl.childPos = 0;
        } else {
          bl.parentId = (nid);
          bl.childPos = i;
        }
        bl.extIndex = externalIndex;
        DEBUG_LOG(
            "updateBorderLeaves: found border leaf parentId=%u childPos=%zu "
            "extIndex=%zu\n",
            bl.parentId, bl.childPos, bl.extIndex);
        externalIndex++;
        borderLeaves_.emplace_back(bl);
      } else {
        stk.emplace_back(cid);
      }
    }
  }
  if (borderLeaves_.size() != numExternalInputs_) {
    DEBUG_LOG(
        "updateBorderLeaves: mismatch in border leaves count %zu vs "
        "numExternalInputs %zu\n",
        borderLeaves_.size(), numExternalInputs_);
    assert(false && "border leaves count mismatch");
  }
  std::sort(
      borderLeaves_.begin(), borderLeaves_.end(),
      [](auto const& a, auto const& b) { return a.extIndex < b.extIndex; });
}

//----------------------------------------------------------------------
// Constructors for tree
//----------------------------------------------------------------------
SNLTruthTableTree::SNLTruthTableTree()
    : rootId_(kInvalidId), numExternalInputs_(0) {
      //initNewBorderLeavesETS();
      //initResolvedChildrenETS();
    }

SNLTruthTableTree::SNLTruthTableTree(naja::DNL::DNLID instid,
                                     naja::DNL::DNLID termid,
                                     Node::Type type) {
  //initNewBorderLeavesETS();
  //initResolvedChildrenETS();
  auto rootNode = std::make_shared<Node>(this, instid, termid, type);
  uint32_t id = allocateNode(rootNode);
  rootId_ = id;

  if (type == Node::Type::P || type == Node::Type::Input) {
    auto inNode = std::make_shared<Node>(0u, this);
    uint32_t inId = allocateNode(inNode);
    rootNode->childrenIds.emplace_back(inId);
    inNode->parentIds.emplace_back(rootId_);
    assert(inNode->parentIds.size() == 1);
    numExternalInputs_ = 1;
    updateBorderLeaves();
    return;
  }

  const auto& table = rootNode->getTruthTable();

  auto arity = table.size();
  for (uint32_t i = 0; i < arity; ++i) {
    auto inNode = std::make_shared<Node>(i, this);
    uint32_t inId = allocateNode(inNode);
    rootNode->childrenIds.emplace_back(inId);
    inNode->parentIds.emplace_back(rootId_);
    assert(inNode->parentIds.size() == 1);
  }
  numExternalInputs_ = arity;
  updateBorderLeaves();
}

//----------------------------------------------------------------------
// size / eval
//----------------------------------------------------------------------
size_t SNLTruthTableTree::size() const {
  return numExternalInputs_;
}

bool SNLTruthTableTree::eval(const std::vector<bool>& extInputs) const {
  if (rootId_ == kInvalidId || extInputs.size() != numExternalInputs_) {
    // LCOV_EXCL_START
    throw std::invalid_argument("wrong input size or uninitialized tree");
    // LCOV_EXCL_STOP
  }
  Node* const rootSp = nodeFromId(rootId_).get();
  if (!rootSp) {
    // LCOV_EXCL_START
    throw std::logic_error("Missing root");
    // LCOV_EXCL_STOP
  }
  return rootSp->eval(extInputs);
}

//----------------------------------------------------------------------
// concatBody
//----------------------------------------------------------------------
const SNLTruthTableTree::Node& SNLTruthTableTree::concatBody(
    size_t borderIndex,
    naja::DNL::DNLID instid,
    naja::DNL::DNLID termid) {
  if (borderIndex >= borderLeaves_.size()) {
    // LCOV_EXCL_START
    throw std::out_of_range("concat: leafIndex out of range");
    // LCOV_EXCL_STOP
  }
  const auto& leaf = borderLeaves_[borderIndex];

  uint32_t parentId = (leaf.parentId);
  auto parentSp = nodeFromId(parentId).get();
  if (!parentSp) {
    // LCOV_EXCL_START
    throw std::logic_error("concat: null parent");
    // LCOV_EXCL_STOP
  }

  uint32_t oldChildId = parentSp->childrenIds[leaf.childPos];

  uint32_t arity = 1;
  std::shared_ptr<Node> newNodeSp;
  if (instid != naja::DNL::DNLID_MAX) {
    newNodeSp = std::make_shared<Node>(this, instid, termid, Node::Type::Table);
    const auto& tbl = newNodeSp->getTruthTable();
    arity = tbl.size();
    auto iter = termid2nodeid_.find(termid);
    if (iter != termid2nodeid_.end()) {
      DEBUG_LOG(
          "###@@@@concat: node for termid %zu %s %s already exists, reusing\n",
          termid,
          naja::DNL::get()
              ->getDNLTerminalFromID(termid)
              .getSnlBitTerm()
              ->getName()
              .getString()
              .c_str(),
          naja::DNL::get()
              ->getDNLTerminalFromID(termid)
              .getDNLInstance()
              .getSNLModel()
              ->getName()
              .getString()
              .c_str());
      // node exist, just connect it to the new parent, but leave the child
      // connections intact
      newNodeSp = nodeFromId(iter->second);
      assert(newNodeSp->type == Node::Type::Table);
      newNodeSp->parentIds.emplace_back(parentId);
      parentSp->childrenIds[leaf.childPos] = newNodeSp->nodeID;
      // assert at least one child for newNodeSp
      if (newNodeSp->childrenIds.size() == 0) {
        // LCOV_EXCL_START
        throw std::logic_error("concat: existing node has no children");
        // LCOV_EXCL_STOP
      }
      return *newNodeSp;
    }
  } else {
    arity = 1;
    newNodeSp = std::make_shared<Node>(this, instid, termid, Node::Type::P);
  }

  uint32_t newNodeId = allocateNode(newNodeSp);

  // Connecting children, skipped if node already existed

  newNodeSp->childrenIds.emplace_back(oldChildId);
  auto oldChildSp = nodeFromId(oldChildId).get();
  if (oldChildSp) {
    assert(oldChildSp->type == Node::Type::Input);
    assert(oldChildSp->parentIds.size() == 1);
    oldChildSp->parentIds[0] = (newNodeId);
    oldChildSp->data.inputIndex = numExternalInputs_;
    numExternalInputs_++;
    DEBUG_LOG("concating with inputIndex %zu\n", oldChildSp->data.inputIndex);
  } else {
    // LCOV_EXCL_START
    throw std::logic_error("concat: null old child");
    // LCOV_EXCL_STOP
  }

  if (newNodeSp->type == Node::Type::Table) {
    for (uint32_t i = 1; i < arity; ++i) {
      auto inNode = std::make_shared<Node>(numExternalInputs_, this);
      numExternalInputs_++;
      uint32_t inId = allocateNode(inNode);
      newNodeSp->childrenIds.emplace_back(inId);
      inNode->parentIds.emplace_back(newNodeId);
      assert(inNode->parentIds.size() == 1);
    }
  }

  parentSp->childrenIds[leaf.childPos] = newNodeId;
  newNodeSp->parentIds.emplace_back(parentId);
  if (!(newNodeSp->parentIds.size() == 1 ||
        newNodeSp->type == Node::Type::Table)) {
    DEBUG_LOG("concat: new node parent count %zu\n",
              newNodeSp->parentIds.size());
    DEBUG_LOG("concat: new node type %s\n",
              newNodeSp->type == Node::Type::Table ? "Table"
              : newNodeSp->type == Node::Type::P   ? "P"
                                                   : "Input");
    assert(newNodeSp->parentIds.size() == 1 ||
           newNodeSp->type == Node::Type::Table);
  }

  return *newNodeSp;
}

void SNLTruthTableTree::concatFull(
    const std::vector<
        std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
        tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>&
        tables, size_t size) {
#ifdef DEBUG_CHECKS
  // print tables
  DEBUG_LOG("Tables in concatFull:\n");
  for (size_t i = 0; i < tables.size(); ++i) {
    DEBUG_LOG("  table %zu termid %zu %s %s\n", i, tables[i].second,
              naja::DNL::get()
                  ->getDNLTerminalFromID(tables[i].second)
                  .getSnlBitTerm()
                  ->getName()
                  .getString()
                  .c_str(),
              naja::DNL::get()
                  ->getDNLTerminalFromID(tables[i].second)
                  .getDNLInstance()
                  .getSNLModel()
                  ->getName()
                  .getString()
                  .c_str());
  }
  // print border leaves
  DEBUG_LOG("Border leaves in concatFull:\n");
  for (const auto& bl : borderLeaves_) {
    auto parentPtr = nodeFromId(bl.parentId);
    if (!parentPtr)
      continue;
    naja::DNL::DNLTerminalFull term =
        naja::DNL::get()->getDNLTerminalFromID(parentPtr->data.termid);
    naja::DNL::DNLInstanceFull inst =
        naja::DNL::get()
            ->getDNLTerminalFromID(parentPtr->data.termid)
            .getDNLInstance();
    DEBUG_LOG("  border leaf instance %s %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              inst.getSNLModel()->getName().getString().c_str());
  }
  // int newInputs = (int)numExternalInputs_;
  std::set<naja::DNL::DNLID> BorderLeafInstances;
  std::set<naja::DNL::DNLID> BorderPIs;
  for (const auto& bl : borderLeaves_) {
    auto parentPtr = nodeFromId(bl.parentId);
    if (parentPtr->type == Node::Type::P) {
      BorderPIs.insert(parentPtr->data.termid);
    }
    if (!parentPtr)
      assert(false);
    if (parentPtr->type == Node::Type::P) {
      // PI table, skip check
      continue;
    }
    naja::DNL::DNLInstanceFull inst =
        naja::DNL::get()
            ->getDNLTerminalFromID(parentPtr->data.termid)
            .getDNLInstance();
    BorderLeafInstances.insert(inst.getSNLInstance()->getID());
  }
  for (auto table : tables) {
    if (table.first == naja::DNL::DNLID_MAX) {
      // PI table, skip check
      continue;
    }
    auto iso = naja::DNL::get()->getDNLIsoDB().getIsoFromIsoIDconst(
        naja::DNL::get()->getDNLTerminalFromID(table.second).getIsoID());
    auto readers = iso.getReaders();
    bool drivingBorderLeaf = false;
    for (auto reader : readers) {
      if (naja::DNL::get()
              ->getDNLTerminalFromID(reader)
              .getDNLInstance()
              .getSNLInstance() == nullptr) {
        continue;
      }
      if (BorderLeafInstances.find(naja::DNL::get()
                                       ->getDNLTerminalFromID(reader)
                                       .getDNLInstance()
                                       .getSNLInstance()
                                       ->getID()) !=
              BorderLeafInstances.end() ||
          BorderPIs.find(table.second) != BorderPIs.end()) {
        drivingBorderLeaf = true;
      }
    }
    if (!drivingBorderLeaf) {
      // print border leaves
      for (const auto& bl : borderLeaves_) {
        auto parentPtr = nodeFromId(bl.parentId);
        if (!parentPtr)
          continue;
        naja::DNL::DNLTerminalFull term =
            naja::DNL::get()->getDNLTerminalFromID(parentPtr->data.termid);
        naja::DNL::DNLInstanceFull inst =
            naja::DNL::get()
                ->getDNLTerminalFromID(parentPtr->data.termid)
                .getDNLInstance();
        DEBUG_LOG("  border leaf instance %zu %s %s\n",
                  inst.getSNLInstance()->getID(),
                  term.getSnlBitTerm()->getName().getString().c_str(),
                  inst.getSNLModel()->getName().getString().c_str());
      }
      DEBUG_LOG(
          "concatFull: table termid %zu %s %s does not drive any border leaf\n",
          table.second,
          naja::DNL::get()
              ->getDNLTerminalFromID(table.second)
              .getSnlBitTerm()
              ->getName()
              .getString()
              .c_str(),
          naja::DNL::get()
              ->getDNLTerminalFromID(table.second)
              .getDNLInstance()
              .getSNLModel()
              ->getName()
              .getString()
              .c_str());
      assert(drivingBorderLeaf &&
             "concatFull: table does not drive any border leaf");
    }
  }
  if (tables.size() > borderLeaves_.size()) {
    // print all tables
    std::set<naja::DNL::DNLID> tableTermIDs;
    for (size_t i = 0; i < tables.size(); ++i) {
      DEBUG_LOG("  table %zu termid %zu %s %s\n", i, tables[i].second,
                naja::DNL::get()
                    ->getDNLTerminalFromID(tables[i].second)
                    .getSnlBitTerm()
                    ->getName()
                    .getString()
                    .c_str(),
                naja::DNL::get()
                    ->getDNLTerminalFromID(tables[i].second)
                    .getDNLInstance()
                    .getSNLModel()
                    ->getName()
                    .getString()
                    .c_str());
      if (tableTermIDs.find(tables[i].second) != tableTermIDs.end()) {
        DEBUG_LOG("concatFull: duplicate table termid %zu %s %s\n",
                  tables[i].second,
                  naja::DNL::get()
                      ->getDNLTerminalFromID(tables[i].second)
                      .getSnlBitTerm()
                      ->getName()
                      .getString()
                      .c_str(),
                  naja::DNL::get()
                      ->getDNLTerminalFromID(tables[i].second)
                      .getDNLInstance()
                      .getSNLModel()
                      ->getName()
                      .getString()
                      .c_str());
      }
      tableTermIDs.insert(tables[i].second);
    }
    DEBUG_LOG("  tableTermIDs %zu\n", tableTermIDs.size());
    DEBUG_LOG("  tables %zu\n", tables.size());
    DEBUG_LOG("  borderLeaves_ %zu\n", borderLeaves_.size());
    assert(tables.size() == tableTermIDs.size() &&
           "concatFull: duplicate tables in input");
    DEBUG_LOG("concatFull: too many tables %zu > %zu\n", tables.size(),
              borderLeaves_.size());
    throw std::invalid_argument("too many tables in concatFull");
  }
#endif
  // FUNC START

  //std::vector<BorderLeaf, tbb::tbb_allocator<BorderLeaf>> newBorderLeaves;
  auto& newBorderLeaves = getNewBorderLeavesETS();
  clearNewBorderLeavesETS();
  size_t newInputs = 0;
  size_t index = 0;
  assert(size == borderLeaves_.size());
  numExternalInputs_ = 0;
  for (size_t i = 0; i < size; ++i) {
    // For each entry in table to merge
    //assert(newBorderLeaves.size() == newInputs);
    assert(getSizeOfNewBorderLeavesETS() == newInputs);
    // Get the relevant border leaf based on order -> assuming identical order
    // between tables and border leaves
    const auto& borderLeaf = borderLeaves_[i];
    // Get parent node of current border leaf
    auto parentPtr = nodeFromId(borderLeaf.parentId).get();
    // if (!parentPtr) {
    //   // No parent so it is the root
    //   index++;
    //   newBorderLeaves.emplace_back(borderLeaf);
    //   DEBUG_LOG("--- concatBody: null parent for border leaf index %zu\n",
    //   index-1); newInputs += 1; assert(newBorderLeaves.size() == newInputs);
    //   assert(rootId_ == borderLeaf.parentId && "concatFull: null parent is
    //   not root"); continue;
    // }
    if (parentPtr->type == Node::Type::P) {
      // If it is a PI border leaf, keep the same leaf and continue, no need to
      // chain PIs
      index++;
      //newBorderLeaves.emplace_back(borderLeaf);
      pushBackNewBorderLeavesETS(borderLeaf);
      DEBUG_LOG("--- concatBody: skipping PI border leaf index %zu\n",
                index - 1);
      newInputs += 1;
      //assert(newBorderLeaves.size() == newInputs);
      assert(getSizeOfNewBorderLeavesETS() == newInputs);
      continue;
    }
    const auto& n = concatBody(index, tables[i].first, tables[i].second);
    if (n.parentIds.size() <= 1 || n.type == Node::Type::P) {
      // if new node is not reused, expand border leaves
      DEBUG_LOG("ConcatBody expanding border leaf index %zu termid %zu %s %s\n",
                index, tables[i].second,
                naja::DNL::get()
                    ->getDNLTerminalFromID(tables[i].second)
                    .getSnlBitTerm()
                    ->getName()
                    .getString()
                    .c_str(),
                naja::DNL::get()
                    ->getDNLTerminalFromID(tables[i].second)
                    .getDNLInstance()
                    .getSNLModel()
                    ->getName()
                    .getString()
                    .c_str());
      // Now we will create new border leaves for each input of the newly
      // inserted node It is in the place of the original border leaf
      uint32_t insertedId = parentPtr->childrenIds[borderLeaf.childPos];
      // assert that insertedId is an input node
      // assert(nodeFromId(insertedId)->type != Node::Type::Input &&
      //  "concatFull: inserted node is input after concatBody");
      auto insertedSp = nodeFromId(insertedId).get();
      assert(insertedSp->type != Node::Type::Input &&
             "concatFull: inserted node is input after concatBody");
      // assert the input node have only one parent
      assert(insertedSp->parentIds.size() == 1 &&
             "concatFull: inserted node has multiple parents after concatBody");
      if (!insertedSp) {
        index++;
        assert(false);
      }
      DEBUG_LOG("insertedSP %s\n",
                naja::DNL::get()
                    ->getDNLTerminalFromID(insertedSp->data.termid)
                    .getSnlBitTerm()
                    ->getName()
                    .getString()
                    .c_str());
      DEBUG_LOG("children count: %zu\n", insertedSp->childrenIds.size());
      // now next is to add border leaf on top of each input node of insertedSp
      for (size_t j = 0; j < insertedSp->childrenIds.size(); ++j) {
        uint32_t cid = insertedSp->childrenIds[j];

        auto ch = nodeFromId(cid).get();
        assert(ch);
        // assert that cid is an input node
        assert(ch->type == Node::Type::Input &&
               "concatFull: inserted node child is not input after concatBody");

        if (ch->type == Node::Type::Input) {
          // Now concat a border leaf for this input
          BorderLeaf bl;
          bl.parentId = (insertedId);
          bl.childPos = j;
          bl.extIndex = ch->data.inputIndex;  // Set correctly in concatBody
          pushBackNewBorderLeavesETS(bl);
          DEBUG_LOG(
              "--- new border leaf extIndex %zu from inserted node id %u "
              "childPos %zu\n",
              bl.extIndex, insertedId, j);
          DEBUG_LOG("--- %s %s\n",
                    naja::DNL::get()
                        ->getDNLTerminalFromID(insertedSp->data.termid)
                        .getSnlBitTerm()
                        ->getName()
                        .getString()
                        .c_str(),
                    naja::DNL::get()
                        ->getDNLTerminalFromID(insertedSp->data.termid)
                        .getDNLInstance()
                        .getSNLModel()
                        ->getName()
                        .getString()
                        .c_str());
          newInputs += 1;
          //assert(newBorderLeaves.size() == newInputs);
          assert(getSizeOfNewBorderLeavesETS() == newInputs); 
        } else {
          assert(false);
        }
      }
    } else {
    }
    index++;
  }
  numExternalInputs_ = (size_t)newInputs;
  borderLeaves_ = std::move(newBorderLeaves.first);
  DEBUG_LOG("ConcatBody done, new numExternalInputs_: %zu\n",
            numExternalInputs_);
  DEBUG_LOG("ConcatBody done, borderLeaves_ size: %zu\n", borderLeaves_.size());
  // CHECKS

#ifdef DEBUG_CHECKS
  // count all inputs and pi nodes in the tree
  std::stack<uint32_t> stk;
  stk.push(rootId_);
  std::set<uint32_t> inputs;
  while (!stk.empty()) {
    uint32_t nid = stk.top();
    stk.pop();
    auto nsp = nodeFromId(nid);
    if (!nsp)
      assert(false && "concatFull: null node in tree during input count");
    for (size_t i = 0; i < nsp->childrenIds.size(); ++i) {
      uint32_t cid = nsp->childrenIds[i];
      auto ch = nodeFromId(cid);
      if (!ch)
        assert(false &&
               "concatFull: null child node in tree during input count");
      assert(std::find(ch->parentIds.begin(), ch->parentIds.end(), nid) !=
                 ch->parentIds.end() &&
             "concatFull: child missing parent link during input count");
      assert(std::find(nsp->childrenIds.begin(), nsp->childrenIds.end(), cid) !=
                 nsp->childrenIds.end() &&
             "concatFull: parent missing child link during input count");
      if (ch->type == Node::Type::Input || ch->type == Node::Type::P) {
        inputs.insert(cid);
      } else {
        stk.push(cid);
      }
    }
  }
  DEBUG_LOG(
      "concatFull: counted inputs %zu vs numExternalInputs_ %zu after "
      "concatFull\n",
      inputs.size(), numExternalInputs_);
  assert((borderLeaves_.size() == numExternalInputs_) &&
         "concatFull: border leaves count mismatch after concatFull");
  for (const auto& bl : borderLeaves_) {
    DEBUG_LOG("1  border leaf parentId %u childPos %zu extIndex %zu\n",
              bl.parentId, bl.childPos, bl.extIndex);
  }
  assert(inputs.size() == numExternalInputs_ &&
         "concatFull: counted inputs mismatch after concatFull");
  // updateBorderLeaves(); <- not used as it changes the order and unsyncing the
  // connection of tables and border leaves order
  assert((borderLeaves_.size() == numExternalInputs_) &&
         "concatFull: border leaves count mismatch after concatFull");
  assert(inputs.size() == numExternalInputs_ &&
         "concatFull: counted inputs mismatch after concatFull");
  for (const auto& bl : borderLeaves_) {
    DEBUG_LOG("2  border leaf parentId %u childPos %zu extIndex %zu\n",
              bl.parentId, bl.childPos, bl.extIndex);
  }

  // assert all new border leaves are in table and in the right order
  size_t order = 0;
  DEBUG_LOG("@@ Border leaves size after concatFull: %zu\n",
            borderLeaves_.size());
  for (size_t i = 0; i < borderLeaves_.size(); ++i) {
    DEBUG_LOG("node id %u border leaf %zu extIndex %zu\n",
              borderLeaves_[i].parentId, i, borderLeaves_[i].extIndex);
    assert(nodeFromId(borderLeaves_[i].parentId) &&
           "concatFull: null border leaf parent after concatFull");
    // assert that node is not an input
    assert(nodeFromId(borderLeaves_[i].parentId)->type != Node::Type::Input &&
           "concatFull: border leaf parent is input after concatFull");
    naja::DNL::DNLID termid =
        naja::DNL::get()
            ->getDNLTerminalFromID(
                nodeFromId(borderLeaves_[i].parentId)->data.termid)
            .getID();
    size_t newOrder = 0;
    DEBUG_LOG("border leaf %zu termid %zu %s %s\n", i, termid,
              naja::DNL::get()
                  ->getDNLTerminalFromID(termid)
                  .getSnlBitTerm()
                  ->getName()
                  .getString()
                  .c_str(),
              naja::DNL::get()
                  ->getDNLTerminalFromID(termid)
                  .getDNLInstance()
                  .getSNLModel()
                  ->getName()
                  .getString()
                  .c_str());
    bool PI = false;
    for (size_t j = 0; j < tables.size(); ++j) {
      if (tables[j].first == naja::DNL::DNLID_MAX) {
        // PI table
        PI = true;
      }
      if (tables[j].second == termid) {
        newOrder = j;
        break;
      }
    }
    if (PI) {
      // skip PI tables in order check
      continue;
    }
    if (newOrder == 0) {
      if (nodeFromId(borderLeaves_[i].parentId)->parentIds.size() > 1) {
        // reused node, skip
        continue;
      }
    }
    DEBUG_LOG("newOrder %zu order %zu\n", newOrder, order);
    assert(newOrder >= order &&
           "concatFull: border leaves out of order after concatFull");
    if (order < newOrder) {
      order = newOrder;
    }
  }
  for (const auto& pair : tables) {
    naja::DNL::DNLID termid = pair.second;
    bool found = false;
    for (size_t i = 0; i < borderLeaves_.size(); ++i) {
      naja::DNL::DNLID btermid =
          naja::DNL::get()
              ->getDNLTerminalFromID(
                  nodeFromId(borderLeaves_[i].parentId)->data.termid)
              .getID();
      if (btermid == termid) {
        found = true;
        break;
      }
    }
    if (!found) {
      DEBUG_LOG(
          "concatFull: table termid %zu %s %s not found in border leaves after "
          "concatFull\n",
          termid,
          naja::DNL::get()
              ->getDNLTerminalFromID(termid)
              .getSnlBitTerm()
              ->getName()
              .getString()
              .c_str(),
          naja::DNL::get()
              ->getDNLTerminalFromID(termid)
              .getDNLInstance()
              .getSNLModel()
              ->getName()
              .getString()
              .c_str());
      if (termid2nodeid_.find(termid) != termid2nodeid_.end()) {
        DEBUG_LOG("  termid %zu exists in termid2nodeid_\n", termid);
      } else {
        assert(false &&
               "concatFull: table not found in border leaves after concatFull");
      }
    }
  }
#endif
}

//----------------------------------------------------------------------
// isInitialized / print
//----------------------------------------------------------------------
bool SNLTruthTableTree::isInitialized() const {
  if (rootId_ == kInvalidId)
    return false;
  std::vector<uint32_t> stk;
  stk.emplace_back(rootId_);
  while (!stk.empty()) {
    uint32_t nid = stk.back();
    stk.pop_back();
    Node* const n = nodeFromId(nid).get();
    if (!n)
      continue;
    if (n->type == Node::Type::Table) {
      if (!n->getTruthTable().isInitialized())
        return false;
    }
    for (size_t i = 0; i < n->childrenIds.size(); ++i) {
      uint32_t cid = n->childrenIds[i];
      Node* const ch = nodeFromId(cid).get();
      if (!ch)
        continue;
      if (ch->type != Node::Type::Input)
        stk.emplace_back(cid);
    }
  }
  return true;
}

// LCOV_EXCL_START
void SNLTruthTableTree::print() const {
  if (rootId_ == kInvalidId)
    return;
  std::vector<uint32_t> stk;
  stk.emplace_back(rootId_);
  while (!stk.empty()) {
    uint32_t nid = stk.back();
    stk.pop_back();
    auto n = nodeFromId(nid);
    if (!n)
      continue;
    if (n->type == Node::Type::Table) {
      printf("term: %zu nodeID=%u id=%u\n", (size_t)n->data.termid,
                n->nodeID, n->nodeID);
    } else if (n->type == Node::Type::P) {
      printf("P nodeID=%u id=%u\n", n->nodeID, n->nodeID);
    } else {
      printf("Input node index=%u nodeID=%u id=%u\n", n->data.inputIndex,
                n->nodeID, n->nodeID);
    }
    for (size_t i = 0; i < n->childrenIds.size(); ++i) {
      uint32_t cid = n->childrenIds[i];
      auto ch = nodeFromId(cid);
      if (!ch) {
        printf("  child[%zu] = null (childId=%u)\n", i, cid);
      } else if (ch->type == Node::Type::Input) {
        printf("  child[%zu] = Input(%u) id=%u\n", i, ch->data.inputIndex,
                  ch->nodeID);
      } else {
        printf("  child[%zu] = Node(id=%u)\n", i, cid);
        stk.emplace_back(cid);
      }
    }
  }
}

// LCOV_EXCL_STOP

//----------------------------------------------------------------------
// destroy
//----------------------------------------------------------------------
void SNLTruthTableTree::destroy() {
  nodes_.clear();
  rootId_ = kInvalidId;
  borderLeaves_.clear();
  numExternalInputs_ = 0;
}

thread_local std::unordered_map<uint32_t, SNLTruthTableTree::Node*, std::hash<uint32_t>, std::equal_to<uint32_t>,
   tbb::tbb_allocator<std::pair<const uint32_t, SNLTruthTableTree::Node*>>> mapById;
thread_local std::unordered_map<uint32_t, SNLTruthTableTree::Node*, std::hash<uint32_t>, std::equal_to<uint32_t>,
   tbb::tbb_allocator<std::pair<const uint32_t, SNLTruthTableTree::Node*>>> mapByNodeID;

// Build reverse map from shared_ptr pointer (address) to canonical id
using MapAlloc = tbb::tbb_allocator<std::pair<const SNLTruthTableTree::Node* const, uint32_t>>;

thread_local std::unordered_map<const SNLTruthTableTree::Node*,
                   uint32_t,
                   std::hash<const SNLTruthTableTree::Node*>,
                   std::equal_to<const SNLTruthTableTree::Node*>,
                   MapAlloc> ptrToId;

//----------------------------------------------------------------------
// finalize: repair and validation after construction
//----------------------------------------------------------------------
void SNLTruthTableTree::finalize() {
  // Build resolver maps for existing nodes based on current fields.
  // We accept that builders may have used:
  //  - correct nodeID values (index+kIdOffset)
  //  - debug nodeID values (node->nodeID)
  //  - temporary or precomputed ids (which may be wrong)
  //
  // Strategy:
  // 1) Create maps: by current nodeID (if valid), by nodeID (debug), and by
  // slot index. 2) For each node, resolve each childrenIds entry by attempting:
  //      a) match by nodeID (fast)
  //      b) match by nodeID (fallback)
  //      c) interpret as index (cid - kIdOffset) within range
  //    If resolved, record target shared_ptr.
  // 3) After all children resolved, rebuild nodes_ canonical ordering (keep
  // existing order),
  //    set node->nodeID = index + kIdOffset, set node->tree = this, and replace
  //    childrenIds with canonical ids derived from the resolved shared_ptrs.
  //
  // This repairs common builder mistakes without requiring edits in builder
  // code.

  // Step 0: quick sanity for root
  if (rootId_ == kInvalidId && nodes_.empty())
    return;
  size_t nodeSize = nodes_.size();
  // Build lookup maps
  
  //mapById.reserve(nodeSize * 2);
  mapById.clear();
  //mapByNodeID.reserve(nodeSize * 2);
  mapByNodeID.clear();

  for (size_t i = 0; i < nodeSize; ++i) {
    Node* sp = nodes_[i].get();
    if (!sp)
      continue;
    if (sp->nodeID != kInvalidId)
      mapById[sp->nodeID] = sp;
    if (sp->nodeID != 0)
      mapByNodeID[sp->nodeID] = sp;
  }
  reserveResolvedChildrenETS(nodeSize);
  for (size_t i = 0; i < nodeSize; ++i) {
    Node* sp = nodes_[i].get();
    if (!sp)
      continue;
    getResolvedChildrenETS().first[i].reserve(sp->childrenIds.size());
    for (size_t j = 0; j < sp->childrenIds.size(); ++j) {
      uint32_t cid = sp->childrenIds[j];
      Node* target = nullptr;

      // try match by exact nodeID
      auto it = mapById.find(cid);
      if (it != mapById.end())
        target = it->second;
      // fallback: match by debug nodeID
      if (!target) {
        auto it2 = mapByNodeID.find(cid);
        if (it2 != mapByNodeID.end())
          target = it2->second;
      }
      // fallback: interpret as index (cid - kIdOffset)
      if (!target) {
        if (cid >= kIdOffset) {
          size_t idx = (size_t)(cid - kIdOffset);
          if (idx < nodeSize) {
            target = nodes_[idx].get();
          }
        }
      }
      if (!target) {
        // LCOV_EXCL_START
        // cannot resolve child id: report and abort
        fprintf(stderr,
                "finalize: could not resolve child reference: parent_slot=%zu "
                "parent_assigned_id=%u childPos=%zu childId=%u nodes=%zu\n",
                i, sp->nodeID, j, cid, nodeSize);
        throw std::logic_error("finalize: unresolved child id");
        // LCOV_EXCL_STOP
      }
      //resolvedChildren[i].emplace_back(target);
      getResolvedChildrenETS().first[i].emplace_back(target);
    }
  }

  // Now assign canonical ids and remap childrenIds/parentId
  for (size_t i = 0; i < nodeSize; ++i) {
    uint32_t canonicalId = static_cast<uint32_t>(i) + kIdOffset;
    Node* sp = nodes_[i].get();
    sp->nodeID = canonicalId;
    sp->tree = this;
  }
  
  //ptrToId.reserve(nodeSize * 2);
  ptrToId.clear();
  for (size_t i = 0; i < nodeSize; ++i) {
    Node*  sp = nodes_[i].get();
    if (!sp)
      continue;
    ptrToId[sp] = static_cast<uint32_t>(i) + kIdOffset;
  }

  // Replace childrenIds with canonical ids and set parentId accordingly
  for (size_t i = 0; i < nodeSize; ++i) {
    Node* sp = nodes_[i].get();
    sp->childrenIds.clear();
    const auto& resolvedChilde = getResolvedChildrenETS().first[i];
    sp->childrenIds.reserve(resolvedChilde.size());
    for (size_t j = 0; j < resolvedChilde.size(); ++j) {
      auto targ = resolvedChilde[j];
      auto it = ptrToId.find(targ);
      if (it == ptrToId.end()) {
        // LCOV_EXCL_START
        fprintf(stderr,
                "finalize: internal error mapping ptr->id parent_slot=%zu "
                "childPos=%zu\n",
                i, j);
        throw std::logic_error("finalize: internal mapping failed");
        // LCOV_EXCL_STOP
      }
      uint32_t newCid = it->second;
      sp->childrenIds.emplace_back(newCid);
      // set child's parentId; last writer wins (ok for tree)
      auto childSp = targ;
      // replace the slot in the index that contain the value i in parentIds
      // with the sp->nodeID
      auto iter = std::find(childSp->parentIds.begin(),
                            childSp->parentIds.end(), sp->nodeID);
      if (iter != childSp->parentIds.end()) {
        *iter = sp->nodeID;
      } else {
        // LCOV_EXCL_START
        throw std::logic_error("finalize: parentIds inconsistent");
        // LCOV_EXCL_STOP
      }
    }
  }

  // Recompute rootId_: if existing rootId_ was resolvable, remap it; otherwise
  // try to keep slot 0
  if (rootId_ != kInvalidId) {
    // try to remap previous rootId_ by matching to new canonical id via
    // mapById/mapByNodeID/slot heuristic
    uint32_t newRoot = kInvalidId;
    auto itRoot = mapById.find(rootId_);
    if (itRoot != mapById.end()) {
      Node* const sp = itRoot->second;
      auto pit = ptrToId.find(sp);
      if (pit != ptrToId.end())
        newRoot = pit->second;
    }
    if (newRoot == kInvalidId) {
      // fallback: if nodes_[0] exists, use that
      if (!nodes_.empty() && nodes_[0])
        newRoot = nodes_[0]->nodeID;
    }
    rootId_ = newRoot;
  }

  // Recompute numExternalInputs_ by scanning leaves
  size_t maxInput = 0;
  numExternalInputs_ = 0;
  bool anyInput = false;
  std::vector<uint32_t, tbb::tbb_allocator<uint32_t>> stk;
  if (rootId_ != kInvalidId) {
    stk.emplace_back(rootId_);
  }
  markAllUnvisited();
  while (!stk.empty()) {
    uint32_t nid = stk.back();
    stk.pop_back();
    //auto [itr, inserted] = visited.insert(nid);
   
    Node* const n = nodeFromId(nid).get();
     if (n->visited)
      continue;
    n->visited = true;
    if (!n)
      continue;
    for (size_t k = 0; k < n->childrenIds.size(); ++k) {
      uint32_t cid = n->childrenIds[k];
      Node* const ch = nodeFromId(cid).get();
      if (!ch)
        continue;
      if (ch->type == Node::Type::Input || ch->type == Node::Type::P) {
        anyInput = true;
        // if (ch->data.inputIndex > maxInput) maxInput = ch->data.inputIndex;
        numExternalInputs_++;
      } else {
        stk.emplace_back(cid);
      }
    }
  }

  updateBorderLeaves();
}
