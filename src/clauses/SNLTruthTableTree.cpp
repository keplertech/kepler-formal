#include "SNLTruthTableTree.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <stack>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <atomic>

using namespace KEPLER_FORMAL;

// Init Ptable holder
const SNLTruthTable SNLTruthTableTree::PtableHolder_ = SNLTruthTable(1,2);

// diagnostic global
static std::atomic<size_t> g_live_nodes{0};

// NodeLifetimeCounter impl
// SNLTruthTableTree::Node::NodeLifetimeCounter::NodeLifetimeCounter()  { g_live_nodes.fetch_add(1, std::memory_order_relaxed); }
// SNLTruthTableTree::Node::NodeLifetimeCounter::~NodeLifetimeCounter() { g_live_nodes.fetch_sub(1, std::memory_order_relaxed); }

//----------------------------------------------------------------------
// Node ctors / dtor
//----------------------------------------------------------------------
SNLTruthTableTree::Node::Node(uint32_t idx, SNLTruthTableTree* t)
  : type(Type::Input), inputIndex(idx), /*nodeID(0),*/ nodeID(SNLTruthTableTree::kInvalidId),
    tree(t), termid(naja::DNL::DNLID_MAX), parentId(SNLTruthTableTree::kInvalidId)
{
  if (tree && tree->lastID_ == std::numeric_limits<unsigned>::max()) {
    throw std::overflow_error("Node ID overflow");
  }
  if (tree) nodeID = (uint32_t)tree->lastID_++;
}

SNLTruthTableTree::Node::Node(SNLTruthTableTree* t,
                              naja::DNL::DNLID instid,
                              naja::DNL::DNLID term,
                              Type type_)
  : type(type_), inputIndex(std::numeric_limits<uint32_t>::max()), /*nodeID(0),*/ nodeID(SNLTruthTableTree::kInvalidId),
    tree(t), termid(term), parentId(SNLTruthTableTree::kInvalidId)
{
  if (tree && tree->lastID_ == std::numeric_limits<unsigned>::max()) {
    throw std::overflow_error("Node ID overflow");
  }
  if (tree) nodeID = (uint32_t)tree->lastID_++;
}

SNLTruthTableTree::Node::~Node() {
  childrenIds.clear();
  parentId = SNLTruthTableTree::kInvalidId;
  tree = nullptr;
}

//----------------------------------------------------------------------
// Node::getTruthTable
//----------------------------------------------------------------------
const SNLTruthTable& SNLTruthTableTree::Node::getTruthTable() const {
  if (type == Type::Table) {
    auto* model = const_cast<SNLDesign*>(naja::DNL::get()->getDNLTerminalFromID(termid).getDNLInstance().getSNLModel());
    return model->getTruthTable(
      naja::DNL::get()->getDNLTerminalFromID(termid)
                   .getSnlBitTerm()->getOrderID());
  }
  else if (type == Type::P || type == Type::Input) {
    return PtableHolder_;
  }
  throw std::logic_error("getTruthTable: not a Table/P node");
}

static std::shared_ptr<SNLTruthTableTree::Node> nullNodePtr = nullptr;

//----------------------------------------------------------------------
// nodeFromId helper
//----------------------------------------------------------------------
const std::shared_ptr<SNLTruthTableTree::Node>& SNLTruthTableTree::nodeFromId(uint32_t id) const {
  if (id == kInvalidId) return nullNodePtr;
  if (id < kIdOffset) return nullNodePtr;
  size_t idx = (size_t)(id - kIdOffset);
  if (idx >= nodes_.size()) return nullNodePtr;
  auto& sp = nodes_[idx];
  if (!sp) return nullNodePtr;
  // sanity check: nodeID must match slot
  if (sp->nodeID != id) {
    fprintf(stderr, "nodeFromId: id mismatch requested=%u slot=%zu node->nodeID=%u\n", id, idx, sp->nodeID);
    return nullNodePtr;
  }
  return sp;
}

//----------------------------------------------------------------------
// Node::eval (resolves children via ids)
//----------------------------------------------------------------------
bool SNLTruthTableTree::Node::eval(const std::vector<bool>& extInputs) const
{
  if (type != Type::Table && type != Type::P && type != Type::Input)
    throw std::logic_error("eval: node not Table/P/Input");

  auto tbl = getTruthTable();
  auto arity = tbl.size();
  if (childrenIds.size() != arity)
    throw std::logic_error("TableNode: children count mismatch");

  uint32_t idx = 0;
  for (uint32_t i = 0; i < arity; ++i) {
    bool bit = false;
    uint32_t cid = childrenIds[i];
    if (cid == kInvalidId) throw std::logic_error("Invalid child id");
    auto childSp = tree->nodeFromId(cid);
    if (!childSp) throw std::logic_error("Null child node");
    if (childSp->type == Type::Input) {
      size_t inx = childSp->inputIndex;
      if (inx >= extInputs.size()) throw std::out_of_range("Input index out of range");
      bit = extInputs[inx];
    } else {
      bit = childSp->eval(extInputs);
    }
    if (bit) idx |= (1u << i);
  }
  return tbl.bits().bit(idx);
}

//----------------------------------------------------------------------
// addChildId: set parent/child relationship via ids
//----------------------------------------------------------------------
void SNLTruthTableTree::Node::addChildId(uint32_t childId) {
  if (childId == kInvalidId) throw std::invalid_argument("addChildId: invalid id");
  uint32_t cur = this->parentId;
  while (cur != SNLTruthTableTree::kInvalidId) {
    if (cur == childId) throw std::invalid_argument("addChildId: cycle detected");
    auto p = tree->nodeFromId(cur);
    if (!p) break;
    cur = p->parentId;
  }

  childrenIds.push_back(childId);

  auto childSp = tree->nodeFromId(childId);
  if (childSp) childSp->parentId = this->nodeID;
}

//----------------------------------------------------------------------
// allocateNode helper - assigns id before publishing into nodes_
//----------------------------------------------------------------------
uint32_t SNLTruthTableTree::allocateNode(const std::shared_ptr<Node>& np) {
  if (!np) throw std::invalid_argument("allocateNode: null");
  uint32_t id = static_cast<uint32_t>(nodes_.size()) + kIdOffset;
  np->nodeID = id;
  np->tree = this;
  nodes_.push_back(np);
  return id;
}

//----------------------------------------------------------------------
// updateBorderLeaves
//----------------------------------------------------------------------
void SNLTruthTableTree::updateBorderLeaves() {
  borderLeaves_.clear();
  if (rootId_ == kInvalidId) return;
  std::vector<uint32_t> stk;
  stk.reserve(64);
  stk.push_back(rootId_);

  while (!stk.empty()) {
    uint32_t nid = stk.back(); stk.pop_back();
    auto nsp = nodeFromId(nid);
    if (!nsp) continue;
    for (size_t i = 0; i < nsp->childrenIds.size(); ++i) {
      uint32_t cid = nsp->childrenIds[i];
      auto ch = nodeFromId(cid);
      if (!ch) continue;
      if (ch->type == Node::Type::Input) {
        BorderLeaf bl;
        bl.parentId = nid;
        bl.childPos = i;
        bl.extIndex = ch->inputIndex;
        borderLeaves_.push_back(bl);
      } else {
        stk.push_back(cid);
      }
    }
  }

  std::sort(borderLeaves_.begin(),
            borderLeaves_.end(),
            [](auto const& a, auto const& b){
              return a.extIndex < b.extIndex;
            });
}

//----------------------------------------------------------------------
// Constructors for tree
//----------------------------------------------------------------------
SNLTruthTableTree::SNLTruthTableTree()
  : rootId_(kInvalidId), numExternalInputs_(0)
{}

SNLTruthTableTree::SNLTruthTableTree(naja::DNL::DNLID instid,
                                     naja::DNL::DNLID termid, Node::Type type)
{
  auto rootNode = std::make_shared<Node>(this, instid, termid, type);
  uint32_t id = allocateNode(rootNode);
  rootId_ = id;

  if (type == Node::Type::P || type == Node::Type::Input) {
    auto inNode = std::make_shared<Node>(0u, this);
    uint32_t inId = allocateNode(inNode);
    rootNode->childrenIds.push_back(inId);
    inNode->parentId = rootId_;
    numExternalInputs_ = 1;
    updateBorderLeaves();
    return;
  }

  auto* model = const_cast<SNLDesign*>(naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel());
  const auto& table = model->getTruthTable(
    naja::DNL::get()->getDNLTerminalFromID(termid)
                 .getSnlBitTerm()->getOrderID());

  auto arity = table.size();
  for (uint32_t i = 0; i < arity; ++i) {
    auto inNode = std::make_shared<Node>(i, this);
    uint32_t inId = allocateNode(inNode);
    rootNode->childrenIds.push_back(inId);
    inNode->parentId = rootId_;
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

bool SNLTruthTableTree::eval(const std::vector<bool>& extInputs) const
{
  if (rootId_ == kInvalidId || extInputs.size() != numExternalInputs_)
    throw std::invalid_argument("wrong input size or uninitialized tree");
  auto rootSp = nodeFromId(rootId_);
  if (!rootSp) throw std::logic_error("Missing root");
  return rootSp->eval(extInputs);
}

//----------------------------------------------------------------------
// concatBody
//----------------------------------------------------------------------
const SNLTruthTableTree::Node&
SNLTruthTableTree::concatBody(size_t borderIndex,
                              naja::DNL::DNLID instid,
                              naja::DNL::DNLID termid)
{
  if (borderIndex >= borderLeaves_.size()) throw std::out_of_range("concat: leafIndex out of range");
  auto leaf = borderLeaves_[borderIndex];

  uint32_t parentId = leaf.parentId;
  auto parentSp = nodeFromId(parentId);
  if (!parentSp) throw std::logic_error("concat: null parent");

  uint32_t oldChildId = parentSp->childrenIds[leaf.childPos];

  uint32_t arity = 1;
  std::shared_ptr<Node> newNodeSp;
  if (instid != naja::DNL::DNLID_MAX) {
    auto* model = const_cast<SNLDesign*>(naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel());
    const auto& tbl = model->getTruthTable(
      naja::DNL::get()->getDNLTerminalFromID(termid)
                   .getSnlBitTerm()->getOrderID());
    arity = tbl.size();
    newNodeSp = std::make_shared<Node>(this, instid, termid, Node::Type::Table);
  } else {
    arity = 1;
    newNodeSp = std::make_shared<Node>(this, instid, termid, Node::Type::P);
  }

  uint32_t newNodeId = allocateNode(newNodeSp);

  newNodeSp->childrenIds.push_back(oldChildId);
  auto oldChildSp = nodeFromId(oldChildId);
  if (oldChildSp) {
    oldChildSp->parentId = newNodeId;
  }

  if (newNodeSp->type == Node::Type::Table) {
    for (uint32_t i = 1; i < arity; ++i) {
      auto inNode = std::make_shared<Node>(numExternalInputs_ + (i - 1), this);
      uint32_t inId = allocateNode(inNode);
      newNodeSp->childrenIds.push_back(inId);
      inNode->parentId = newNodeId;
    }
  }

  parentSp->childrenIds[leaf.childPos] = newNodeId;
  newNodeSp->parentId = parentId;

  return *newNodeSp;
}

//----------------------------------------------------------------------
// concat / concatFull
//----------------------------------------------------------------------
void SNLTruthTableTree::concat(size_t borderIndex,
                               naja::DNL::DNLID instid,
                               naja::DNL::DNLID termid)
{
  auto const& n = concatBody(borderIndex,instid,termid);
  numExternalInputs_ += (n.getTruthTable().size() - 1);
  updateBorderLeaves();
}

void SNLTruthTableTree::concatFull(
  const std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
            tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>& tables)
{
  int newInputs = (int)numExternalInputs_;
  if (tables.size() > borderLeaves_.size())
    throw std::invalid_argument("too many tables in concatFull");
  std::vector<BorderLeaf, tbb::tbb_allocator<BorderLeaf>> newBorderLeaves;
  size_t index = 0;
  for (size_t i=0; i<tables.size();++i) {
    auto borderLeaf = borderLeaves_[i];
    auto parentPtr = nodeFromId(borderLeaf.parentId);
    if (!parentPtr) {
      index++;
      newBorderLeaves.push_back(borderLeaf);
      continue;
    }
    if (parentPtr->type == Node::Type::P) {
      index++;
      newBorderLeaves.push_back(borderLeaf);
      continue;
    }

    const auto& n = concatBody(index, tables[i].first, tables[i].second);
    newInputs += (n.getTruthTable().size() - 1);

    uint32_t insertedId = parentPtr->childrenIds[borderLeaf.childPos];
    auto insertedSp = nodeFromId(insertedId);
    if (!insertedSp) { index++; continue; }

    for (size_t j = 0; j < insertedSp->childrenIds.size(); ++j) {
      uint32_t cid = insertedSp->childrenIds[j];
      auto ch = nodeFromId(cid);
      if (!ch) continue;
      if (ch->type == Node::Type::Input) {
        BorderLeaf bl;
        bl.parentId = insertedId;
        bl.childPos = j;
        bl.extIndex = ch->inputIndex;
        newBorderLeaves.push_back(bl);
      } else {
        for (size_t k = 0; k < ch->childrenIds.size(); ++k) {
          uint32_t ccid = ch->childrenIds[k];
          auto cc = nodeFromId(ccid);
          if (!cc) continue;
          if (cc->type == Node::Type::Input) {
            BorderLeaf bl2;
            bl2.parentId = cid;
            bl2.childPos = k;
            bl2.extIndex = cc->inputIndex;
            newBorderLeaves.push_back(bl2);
          }
        }
      }
    }
    index++;
  }
  numExternalInputs_ = (size_t) newInputs;
  borderLeaves_ = std::move(newBorderLeaves);
}

//----------------------------------------------------------------------
// isInitialized / print
//----------------------------------------------------------------------
bool SNLTruthTableTree::isInitialized() const {
  if (rootId_ == kInvalidId) return false;
  std::vector<uint32_t> stk;
  stk.push_back(rootId_);
  while(!stk.empty()) {
    uint32_t nid = stk.back(); stk.pop_back();
    auto n = nodeFromId(nid);
    if (!n) continue;
    if (n->type == Node::Type::Table) {
      if (!n->getTruthTable().isInitialized()) return false;
    }
    for (size_t i = 0; i < n->childrenIds.size(); ++i) {
      uint32_t cid = n->childrenIds[i];
      auto ch = nodeFromId(cid);
      if (!ch) continue;
      if (ch->type != Node::Type::Input) stk.push_back(cid);
    }
  }
  return true;
}

void SNLTruthTableTree::print() const {
  if (rootId_ == kInvalidId) return;
  std::vector<uint32_t> stk;
  stk.push_back(rootId_);
  while (!stk.empty()) {
    uint32_t nid = stk.back(); stk.pop_back();
    auto n = nodeFromId(nid);
    if (!n) continue;
    if (n->type == Node::Type::Table) {
      printf("term: %zu nodeID=%u id=%u\n", (size_t)n->termid, n->nodeID, n->nodeID);
    } else if (n->type == Node::Type::P) {
      printf("P nodeID=%u id=%u\n", n->nodeID, n->nodeID);
    } else {
      printf("Input node index=%u nodeID=%u id=%u\n", n->inputIndex, n->nodeID, n->nodeID);
    }
    for (size_t i = 0; i < n->childrenIds.size(); ++i) {
      uint32_t cid = n->childrenIds[i];
      auto ch = nodeFromId(cid);
      if (!ch) {
        printf("  child[%zu] = null (childId=%u)\n", i, cid);
      } else if (ch->type == Node::Type::Input) {
        printf("  child[%zu] = Input(%u) id=%u\n", i, ch->inputIndex, ch->nodeID);
      } else {
        printf("  child[%zu] = Node(id=%u)\n", i, cid);
        stk.push_back(cid);
      }
    }
  }
}

//----------------------------------------------------------------------
// simplify
//----------------------------------------------------------------------
void SNLTruthTableTree::simplify() {
  if (rootId_ == kInvalidId) return;

  std::vector<uint32_t> stackIds;
  stackIds.reserve(nodes_.size());
  stackIds.push_back(rootId_);

  std::vector<uint32_t> order;
  order.reserve(nodes_.size());
  std::unordered_set<uint32_t> seen;

  while (!stackIds.empty()) {
    uint32_t nid = stackIds.back(); stackIds.pop_back();
    if (seen.count(nid)) continue;
    seen.insert(nid);
    auto n = nodeFromId(nid);
    if (!n) continue;
    for (uint32_t cid : n->childrenIds) {
      auto ch = nodeFromId(cid);
      if (ch && ch->type != Node::Type::Input) stackIds.push_back(cid);
    }
    order.push_back(nid);
  }

  for (auto nid : order) {
    auto node = nodeFromId(nid);
    if (!node) continue;
    if (node->type != Node::Type::Table) continue;

    const SNLTruthTable &tbl = node->getTruthTable();
    const uint32_t arity = tbl.size();

    if (node->childrenIds.size() != arity) continue;

    if (arity == 1) {
      bool b0 = tbl.bits().bit(0);
      bool b1 = tbl.bits().bit(1);
      if (b0 == false && b1 == true) {
        uint32_t childId = node->childrenIds[0];
        if (node->parentId != kInvalidId) {
          auto parent = nodeFromId(node->parentId);
          if (parent) {
            for (size_t i = 0; i < parent->childrenIds.size(); ++i) {
              if (parent->childrenIds[i] == nid) {
                parent->childrenIds[i] = childId;
                auto child = nodeFromId(childId);
                if (child) child->parentId = parent->nodeID;
                break;
              }
            }
          }
        } else {
          rootId_ = childId;
          auto child = nodeFromId(childId);
          if (child) child->parentId = kInvalidId;
        }
        continue;
      }
    }

    bool all_equal = true;
    for (size_t i = 1; i < node->childrenIds.size(); ++i) {
      if (node->childrenIds[i] != node->childrenIds[0]) { all_equal = false; break; }
    }
    if (all_equal && arity >= 1) {
      uint32_t idx0 = 0;
      uint32_t idx1 = (1u << arity) - 1u;
      bool out0 = tbl.bits().bit(idx0);
      bool out1 = tbl.bits().bit(idx1);
      if (!out0 && out1) {
        uint32_t childId = node->childrenIds[0];
        if (node->parentId != kInvalidId) {
          auto parent = nodeFromId(node->parentId);
          if (parent) {
            for (size_t i = 0; i < parent->childrenIds.size(); ++i) {
              if (parent->childrenIds[i] == nid) {
                parent->childrenIds[i] = childId;
                auto child = nodeFromId(childId);
                if (child) child->parentId = parent->nodeID;
                break;
              }
            }
          }
        } else {
          rootId_ = childId;
          auto child = nodeFromId(childId);
          if (child) child->parentId = kInvalidId;
        }
        continue;
      }
    }
  }

  size_t maxInput = 0;
  bool anyInput = false;
  std::vector<uint32_t> stk2;
  if (rootId_ != kInvalidId) stk2.push_back(rootId_);
  while (!stk2.empty()) {
    uint32_t nid = stk2.back(); stk2.pop_back();
    auto n = nodeFromId(nid);
    if (!n) continue;
    for (size_t i = 0; i < n->childrenIds.size(); ++i) {
      uint32_t cid = n->childrenIds[i];
      auto ch = nodeFromId(cid);
      if (!ch) continue;
      if (ch->type == Node::Type::Input) {
        anyInput = true;
        if (ch->inputIndex > maxInput) maxInput = ch->inputIndex;
      } else {
        stk2.push_back(cid);
      }
    }
  }
  if (anyInput) numExternalInputs_ = maxInput + 1;
  else numExternalInputs_ = 0;

  updateBorderLeaves();
}

//----------------------------------------------------------------------
// destroy
//----------------------------------------------------------------------
void SNLTruthTableTree::destroy() {
  nodes_.clear();
  rootId_ = kInvalidId;
  borderLeaves_.clear();
  numExternalInputs_ = 0;
}

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
  // 1) Create maps: by current nodeID (if valid), by nodeID (debug), and by slot index.
  // 2) For each node, resolve each childrenIds entry by attempting:
  //      a) match by nodeID (fast)
  //      b) match by nodeID (fallback)
  //      c) interpret as index (cid - kIdOffset) within range
  //    If resolved, record target shared_ptr.
  // 3) After all children resolved, rebuild nodes_ canonical ordering (keep existing order),
  //    set node->nodeID = index + kIdOffset, set node->tree = this, and replace childrenIds with
  //    canonical ids derived from the resolved shared_ptrs.
  //
  // This repairs common builder mistakes without requiring edits in builder code.

  // Step 0: quick sanity for root
  if (rootId_ == kInvalidId && nodes_.empty()) return;

  // Build lookup maps
  std::unordered_map<uint32_t, std::shared_ptr<Node>> mapById;
  std::unordered_map<uint32_t, std::shared_ptr<Node>> mapByNodeID;
  mapById.reserve(nodes_.size() * 2);
  mapByNodeID.reserve(nodes_.size() * 2);

  for (size_t i = 0; i < nodes_.size(); ++i) {
    auto sp = nodes_[i];
    if (!sp) continue;
    if (sp->nodeID != kInvalidId) mapById[sp->nodeID] = sp;
    if (sp->nodeID != 0) mapByNodeID[sp->nodeID] = sp;
  }

  // Resolve children entries to shared_ptrs for every node
  std::vector<std::vector<std::shared_ptr<Node>>> resolvedChildren(nodes_.size());
  for (size_t i = 0; i < nodes_.size(); ++i) {
    auto sp = nodes_[i];
    if (!sp) continue;
    resolvedChildren[i].reserve(sp->childrenIds.size());
    for (size_t j = 0; j < sp->childrenIds.size(); ++j) {
      uint32_t cid = sp->childrenIds[j];
      std::shared_ptr<Node> target;

      // try match by exact nodeID
      auto it = mapById.find(cid);
      if (it != mapById.end()) target = it->second;
      // fallback: match by debug nodeID
      if (!target) {
        auto it2 = mapByNodeID.find(cid);
        if (it2 != mapByNodeID.end()) target = it2->second;
      }
      // fallback: interpret as index (cid - kIdOffset)
      if (!target) {
        if (cid >= kIdOffset) {
          size_t idx = (size_t)(cid - kIdOffset);
          if (idx < nodes_.size()) {
            target = nodes_[idx];
          }
        }
      }
      if (!target) {
        // cannot resolve child id: report and abort
        fprintf(stderr, "finalize: could not resolve child reference: parent_slot=%zu parent_assigned_id=%u childPos=%zu childId=%u nodes=%zu\n",
                i, sp->nodeID, j, cid, nodes_.size());
        throw std::logic_error("finalize: unresolved child id");
      }
      resolvedChildren[i].push_back(target);
    }
  }

  // Now assign canonical ids and remap childrenIds/parentId
  for (size_t i = 0; i < nodes_.size(); ++i) {
    uint32_t canonicalId = static_cast<uint32_t>(i) + kIdOffset;
    auto sp = nodes_[i];
    sp->nodeID = canonicalId;
    sp->tree = this;
  }

  // Build reverse map from shared_ptr pointer (address) to canonical id
  std::unordered_map<const Node*, uint32_t> ptrToId;
  ptrToId.reserve(nodes_.size()*2);
  for (size_t i = 0; i < nodes_.size(); ++i) {
    auto sp = nodes_[i];
    if (!sp) continue;
    ptrToId[sp.get()] = static_cast<uint32_t>(i) + kIdOffset;
  }

  // Replace childrenIds with canonical ids and set parentId accordingly
  for (size_t i = 0; i < nodes_.size(); ++i) {
    auto sp = nodes_[i];
    sp->childrenIds.clear();
    sp->childrenIds.reserve(resolvedChildren[i].size());
    for (size_t j = 0; j < resolvedChildren[i].size(); ++j) {
      auto targ = resolvedChildren[i][j];
      auto it = ptrToId.find(targ.get());
      if (it == ptrToId.end()) {
        fprintf(stderr, "finalize: internal error mapping ptr->id parent_slot=%zu childPos=%zu\n", i, j);
        throw std::logic_error("finalize: internal mapping failed");
      }
      uint32_t newCid = it->second;
      sp->childrenIds.push_back(newCid);
      // set child's parentId; last writer wins (ok for tree)
      auto childSp = targ;
      childSp->parentId = sp->nodeID;
    }
  }

  // Recompute rootId_: if existing rootId_ was resolvable, remap it; otherwise try to keep slot 0
  if (rootId_ != kInvalidId) {
    // try to remap previous rootId_ by matching to new canonical id via mapById/mapByNodeID/slot heuristic
    uint32_t newRoot = kInvalidId;
    auto itRoot = mapById.find(rootId_);
    if (itRoot != mapById.end()) {
      auto sp = itRoot->second;
      auto pit = ptrToId.find(sp.get());
      if (pit != ptrToId.end()) newRoot = pit->second;
    }
    if (newRoot == kInvalidId) {
      // fallback: if nodes_[0] exists, use that
      if (!nodes_.empty() && nodes_[0]) newRoot = nodes_[0]->nodeID;
    }
    rootId_ = newRoot;
  }

  // Recompute numExternalInputs_ by scanning leaves
  size_t maxInput = 0;
  bool anyInput = false;
  std::vector<uint32_t> stk;
  if (rootId_ != kInvalidId) stk.push_back(rootId_);
  while (!stk.empty()) {
    uint32_t nid = stk.back(); stk.pop_back();
    auto n = nodeFromId(nid);
    if (!n) continue;
    for (size_t k = 0; k < n->childrenIds.size(); ++k) {
      uint32_t cid = n->childrenIds[k];
      auto ch = nodeFromId(cid);
      if (!ch) continue;
      if (ch->type == Node::Type::Input) {
        anyInput = true;
        if (ch->inputIndex > maxInput) maxInput = ch->inputIndex;
      } else {
        stk.push_back(cid);
      }
    }
  }
  if (anyInput) numExternalInputs_ = maxInput + 1;
  else numExternalInputs_ = 0;

  updateBorderLeaves();
}
