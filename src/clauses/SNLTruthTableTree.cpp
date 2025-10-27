#include "SNLTruthTableTree.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>

//#define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

//#define DEBUG_CHECKS

using namespace naja::NL;
using namespace KEPLER_FORMAL;

// Init Ptable holder
const SNLTruthTable SNLTruthTableTree::PtableHolder_ = SNLTruthTable(1,2);

//----------------------------------------------------------------------
// Node::addChild
//----------------------------------------------------------------------
void SNLTruthTableTree::Node::addChild(const std::shared_ptr<Node>& child) {
  // cycle detection walking parent chain
  auto self = shared_from_this();
  #ifdef DEBUG_CHECKS
  for (auto node = self; node; node = node->parent.lock()) {
    if (node->type == Type::Table
        && node->dnlid == child->dnlid
        && node->termid == child->termid
        && child->dnlid != naja::DNL::DNLID_MAX
        && child->termid != naja::DNL::DNLID_MAX) {
      throw std::invalid_argument("addChild: cycle detected");
    }
  }
  #endif

  // attach
  children.push_back(std::move(child));
  children.back()->parent = self;
}

//----------------------------------------------------------------------
// Node::getTruthTable
//----------------------------------------------------------------------
const SNLTruthTable& SNLTruthTableTree::Node::getTruthTable() const {
  if (type == Type::Table) {
    auto* model = const_cast<SNLDesign*>(naja::DNL::get()->getDNLInstanceFromID(dnlid).getSNLModel());
    return model->getTruthTable(
      naja::DNL::get()->getDNLTerminalFromID(termid)
                   .getSnlBitTerm()->getOrderID());
  }
  else if (type == Type::P) {
    return PtableHolder_;
  }
  throw std::logic_error("getTruthTable: not a Table/P node");
}

//----------------------------------------------------------------------
// Node::eval
//----------------------------------------------------------------------
bool SNLTruthTableTree::Node::eval(
  const std::vector<bool>& extInputs) const 
{
  if (type == Type::Input) {
    if (inputIndex >= extInputs.size())
      throw std::out_of_range("InputNode: index out of range");
    return extInputs[inputIndex];
  }

  auto tbl = getTruthTable();
  auto arity = tbl.size();
  if (children.size() != arity)
    throw std::logic_error("TableNode: children count mismatch");

  uint32_t idx = 0;
  for (uint32_t i = 0; i < arity; ++i) {
    if (children[i]->eval(extInputs))
      idx |= (1u << i);
  }
  return tbl.bits().bit(idx);
}

//----------------------------------------------------------------------
// Tree ctors
//----------------------------------------------------------------------
SNLTruthTableTree::SNLTruthTableTree()
  : root_(nullptr), numExternalInputs_(0) 
{}

// SNLTruthTableTree::SNLTruthTableTree(Node::Type type) {
//   if (type != Node::Type::P)
//     throw std::invalid_argument("invalid type ctor");

//   root_ = std::make_shared<Node>(this, naja::DNL::DNLID_MAX, naja::DNL::DNLID_MAX, Node::Type::P);
//   root_->addChild(std::make_shared<Node>(0, this));
//   numExternalInputs_ = 1;
// }

SNLTruthTableTree::SNLTruthTableTree(naja::DNL::DNLID instid,
                                     naja::DNL::DNLID termid, Node::Type type)
{
  root_ = std::make_shared<Node>(this, instid, termid, type);
  if (type == Node::Type::P) {
    root_->addChild(std::make_shared<Node>(0, this));
    numExternalInputs_ = 1;
    return;
  }
  auto* model = const_cast<SNLDesign*>(naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel());
  const auto& table = model->getTruthTable(
    naja::DNL::get()->getDNLTerminalFromID(termid)
                 .getSnlBitTerm()->getOrderID());

  auto arity = table.size();
  
  for (uint32_t i = 0; i < arity; ++i)
    root_->addChild(std::make_shared<Node>(i, this));

  numExternalInputs_ = arity;
  updateBorderLeaves();
}

//----------------------------------------------------------------------
// size / eval
//----------------------------------------------------------------------
size_t SNLTruthTableTree::size() const {
  return numExternalInputs_;
}

bool SNLTruthTableTree::eval(
  const std::vector<bool>& extInputs) const 
{
  if (!root_ || extInputs.size() != numExternalInputs_)
    throw std::invalid_argument("wrong input size or uninitialized tree");
  return root_->eval(extInputs);
}

//----------------------------------------------------------------------
// updateBorderLeaves
//----------------------------------------------------------------------
void SNLTruthTableTree::updateBorderLeaves() {
  borderLeaves_.clear();
  struct Frame { Node* parent; size_t childPos; Node* node; };
  std::vector<Frame, tbb::tbb_allocator<Frame>> stk;
  stk.reserve(64);
  stk.push_back({nullptr,0,root_.get()});

  while (!stk.empty()) {
    auto f = stk.back(); stk.pop_back();
    if (f.node->type == Node::Type::Table) {
      for (size_t i = f.node->children.size(); i-- > 0;) {
        stk.push_back({f.node,i,f.node->children[i].get()});
      }
    } else if (f.node->type == Node::Type::Input) {
      borderLeaves_.push_back({f.parent,f.childPos,f.node->inputIndex});
    }
  }

  std::sort(borderLeaves_.begin(),
            borderLeaves_.end(),
            [](auto const& a, auto const& b){
              return a.extIndex < b.extIndex;
            });
}

//----------------------------------------------------------------------
// core graft logic
//----------------------------------------------------------------------
const SNLTruthTableTree::Node& 
SNLTruthTableTree::concatBody(size_t borderIndex,
                              naja::DNL::DNLID instid,
                              naja::DNL::DNLID termid)
{
  DEBUG_LOG("--- cocating index %zu with instid/termid %zu/%zu\n", borderIndex, (size_t)instid, (size_t)termid);
  DEBUG_LOG("--- cocating Instance name: %s\n",
         naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel()->getName().getString().c_str());
  DEBUG_LOG("--- cocating Terminal name: %s\n", 
         naja::DNL::get()->getDNLTerminalFromID(termid).getSnlBitTerm()->getName().getString().c_str());
  DEBUG_LOG("--- cocating  path: %s\n", 
         naja::DNL::get()->getDNLInstanceFromID(instid).getFullPath().c_str());
  // locate slot
  // size_t idx = borderIndex;
  // if (idx >= borderLeaves_.size()) {
  //   bool found=false;
  //   for (size_t i=0;i<borderLeaves_.size();++i)
  //     if (borderLeaves_[i].extIndex == borderIndex) {
  //       idx=i; found=true; break;
  //     }
  //   if (!found) throw std::out_of_range("concat: leafIndex out of range");
  // }
  auto leaf = borderLeaves_[borderIndex];
  // get iso of termid
#ifdef DEBUG_CHECKS
  const naja::DNL::DNLIso& iso = naja::DNL::get()->getDNLIsoDB().getIsoFromIsoIDconst(naja::DNL::get()->getIsoIdfromTermId(termid));
  naja::DNL::DNLID readeInst = leaf.parent ? leaf.parent->dnlid : naja::DNL::DNLID_MAX;
  bool validReader = false;
  if (readeInst != naja::DNL::DNLID_MAX) {
    for (auto reader : iso.getReaders()) {
      DEBUG_LOG("reader term id %zu inst id %zu\n", (size_t)reader, (size_t)naja::DNL::get()->getDNLTerminalFromID(reader).getDNLInstance().getID());
      DEBUG_LOG("reader path %s\n", naja::DNL::get()->getDNLTerminalFromID(reader).getDNLInstance().getFullPath().c_str());
      if (naja::DNL::get()->getDNLTerminalFromID(reader).getDNLInstance().getID() == readeInst) {
        validReader = true;
        break;
      }
    }
  }
  if (!(readeInst == naja::DNL::DNLID_MAX || validReader)) {
    auto newNode = leaf.parent;
    for (Node *node = newNode; node; node = node->parent.lock().get()) {
      DEBUG_LOG("check node with type %d: instid/termid %zu/%zu\n", node->type, (size_t)node->dnlid, (size_t)node->termid);
    }
    assert(false);
  }
  Node * parent = nullptr;
  // assert going up the parents that there is no loop, hence no node with the same instid/termid
  for (Node *node = leaf.parent; node; node = node->parent.lock().get()) {
    parent = node;
    if (node->type == Node::Type::Table
        && node->dnlid == instid
        && node->termid == termid
        && instid != naja::DNL::DNLID_MAX
        && termid != naja::DNL::DNLID_MAX) {
      // print the current tree
      //print();
      //assert termid is output
      assert(naja::DNL::get()->getDNLTerminalFromID(termid).getSnlBitTerm()->getDirection() == SNLBitTerm::Direction::Output);
      DEBUG_LOG("inst: %s\n", naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel()->getName().getString().c_str());
      DEBUG_LOG("term: %s\n", naja::DNL::get()->getDNLTerminalFromID(termid).getSnlBitTerm()->getName().getString().c_str());
      Node * parent1 = nullptr;
      for (Node *node1 = leaf.parent; node1; node1 = node1->parent.lock().get()) {
        parent1 = node1;
        DEBUG_LOG("parent node with type %d: instid/termid %zu/%zu\n", node1->type, (size_t)node1->dnlid, (size_t)node1->termid);
      } 
      assert(root_.get() == parent1);
      DEBUG_LOG("looped elemed instance %s term %s\n with path %s\n", 
             naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel()->getName().getString().c_str(),
             naja::DNL::get()->getDNLTerminalFromID(termid).getSnlBitTerm()->getName().getString().c_str(),
             naja::DNL::get()->getDNLInstanceFromID(instid).getFullPath().c_str());
      // print the path till the root1qŵ1 
      for (Node *node1 = leaf.parent; node1; node1 = node1->parent.lock().get()) {
        DEBUG_LOG("path %s \n", 
               naja::DNL::get()->getDNLInstanceFromID(node1->dnlid).getFullPath().c_str());
      }
      throw std::invalid_argument("concat: cycle detected for instid/termid" + std::to_string((size_t)instid)
                                  + "/" + std::to_string((size_t)termid) );
    }
  }
  DEBUG_LOG("root node with type %d: instid/termid %zu/%zu\n", root_.get()->type, (size_t)root_.get()->dnlid, (size_t)root_.get()->termid);
  DEBUG_LOG("parent node with type %d: instid/termid %zu/%zu\n", parent->type, (size_t)parent->dnlid, (size_t)parent->termid);
  assert(root_.get() == parent);
#endif
  // detach old leaf
  std::shared_ptr<Node> oldLeaf;
  if (leaf.parent) {
    oldLeaf = std::move(leaf.parent->children[leaf.childPos]);
  } else {
    oldLeaf = std::move(root_);
  }

  // create new table/P node
  uint32_t arity=0;
  std::shared_ptr<Node> newNode;
  const naja::DNL::DNLTerminalFull& term = naja::DNL::get()->getDNLTerminalFromID(termid);
  if (instid!=naja::DNL::DNLID_MAX) { // the condition works but need to be refactored to if P node 
    auto* model = const_cast<SNLDesign*>(naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel());
    const auto& tbl = model->getTruthTable(
      naja::DNL::get()->getDNLTerminalFromID(termid)
                   .getSnlBitTerm()->getOrderID());
    arity = tbl.size();
    newNode = std::make_shared<Node>(this, instid, termid);
  } else {
    //auto tbl = SNLTruthTable(1,2);
    arity = 1;
    newNode = std::make_shared<Node>(this, instid,  termid, Node::Type::P);
  }

  // re-use old leaf as first child, then fresh inputs
  if (arity>0) {
    newNode->addChild(std::move(oldLeaf));
    for (uint32_t i = 1; i<arity; ++i)
      newNode->addChild(
        std::make_shared<Node>(numExternalInputs_ + (i-1), this)
      );
  }

  // reattach
  if (leaf.parent) {
    leaf.parent->children[leaf.childPos] = newNode;
    newNode->parent = leaf.parent ? leaf.parent->shared_from_this() : std::weak_ptr<Node>();
  } else {
    root_ = newNode;
  }
#ifdef DEBUG_CHECKS
  // assert that newNode is or can reach root through parents
  Node* check = newNode.get();
  for (Node *node = newNode.get(); node; node = node->parent.lock().get()) {
    DEBUG_LOG("check node with type %d: instid/termid %zu/%zu\n", node->type, (size_t)node->dnlid, (size_t)node->termid);
    check = node;
  }
  DEBUG_LOG("root node with type %d: instid/termid %zu/%zu\n", root_.get()->type, (size_t)root_.get()->dnlid, (size_t)root_.get()->termid);
  DEBUG_LOG("check node with type %d: instid/termid %zu/%zu\n", check->type, (size_t)check->dnlid, (size_t)check->termid);
  assert(root_.get() == check);
#endif
  return *(leaf.parent
           ? leaf.parent->children[leaf.childPos].get()
           : root_.get());
}

//----------------------------------------------------------------------
// public concat APIs
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
    // get border leaf for i 
    auto borderLeaf = borderLeaves_[i];
    if (borderLeaf.parent->dnlid == naja::DNL::DNLID_MAX) {
      // no need to continue to expend on P type border leaves
      index++;
      newBorderLeaves.push_back(borderLeaf);
      continue;
    }
    auto const& n = concatBody(index, tables[i].first, tables[i].second);
    assert(n.getTruthTable().size() > 0 || newInputs > 0);
    newInputs += (n.getTruthTable().size() - 1);
    for (size_t i = 0; i < n.children.size(); ++i) {
      auto child = n.children[i].get();
      if (child->type == Node::Type::Input) {
        newBorderLeaves.push_back({child->parent.lock().get(), i, child->inputIndex});
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
  if (!root_) return false;
  struct Frame{Node* n;};
  std::vector<Frame, tbb::tbb_allocator<Frame>> stk{{root_.get()}};
  while(!stk.empty()) {
    auto f = stk.back(); stk.pop_back();
    if (f.n->type == Node::Type::Table) {
      if (!f.n->getTruthTable().isInitialized()) return false;
      for (auto& c : f.n->children)
        stk.push_back({c.get()});
    }
  }
  return true;
}

void SNLTruthTableTree::print() const {
  if (!root_) return;
  struct Frame{Node* n; size_t d;};
  std::vector<Frame, tbb::tbb_allocator<Frame>> stk{{root_.get(), 0}};
  while(!stk.empty()) {
    auto [n,d] = stk.back(); stk.pop_back();
    // print node content
    DEBUG_LOG("Node: Type=%d dnlid=%zu termid=%zu nodeID=%zu\n"
      , (int)n->type, (size_t)n->dnlid, (size_t)n->termid, (size_t)n->nodeID);
    // print indentation + node info here…
    if (n->type == Node::Type::Table) {
      DEBUG_LOG("instance: %s\n", naja::DNL::get()->getDNLInstanceFromID(n->dnlid).getSNLModel()->getName().getString().c_str());
      DEBUG_LOG("term: %s\n", naja::DNL::get()->getDNLTerminalFromID(n->termid).getSnlBitTerm()->getName().getString().c_str());
      for (size_t i=n->children.size(); i-->0;)
        stk.push_back({n->children[i].get(), d+1});
    }
  }
}

//----------------------------------------------------------------------
// Collapse the tree based on simple, local constant/identity rules
//----------------------------------------------------------------------
// Conservative transforms only: identity unary-table removal, duplicate-child collapse,
// and full-equal-children simplification where safe.
void SNLTruthTableTree::simplify() {
  if (!root_) return;

  // Post-order traversal collect nodes
  struct Frame { std::shared_ptr<Node> ptr; bool visited; };
  std::vector<Frame, tbb::tbb_allocator<Frame>> stack;
  stack.reserve(256);
  stack.push_back({root_, false});

  std::vector<std::shared_ptr<Node>, tbb::tbb_allocator<std::shared_ptr<Node>>> order;
  order.reserve(256);

  while (!stack.empty()) {
    auto fr = stack.back(); stack.pop_back();
    if (!fr.ptr) continue;
    if (fr.visited) {
      order.push_back(fr.ptr);
      continue;
    }
    // push again as visited then children
    stack.push_back({fr.ptr, true});
    if (fr.ptr->type == Node::Type::Table) {
      for (auto &c : fr.ptr->children)
        stack.push_back({c, false});
    }
  }

  // Helper to replace child pointer in parent
  auto replaceChildInParent = [&](const std::shared_ptr<Node>& target,
                                  const std::shared_ptr<Node>& replacement) {
    auto wp = target->parent;
    if (wp.expired()) {
      // target was root
      root_ = replacement;
      if (replacement) replacement->parent.reset();
      return;
    }
    auto parent = wp.lock();
    for (size_t i = 0; i < parent->children.size(); ++i) {
      if (parent->children[i].get() == target.get()) {
        parent->children[i] = replacement;
        if (replacement) replacement->parent = parent;
        return;
      }
    }
  };

  // Process nodes bottom-up
  for (auto &node_shared : order) {
    Node* node = node_shared.get();
    if (node->type != Node::Type::Table) continue;

    // retrieve table
    const SNLTruthTable &tbl = node->getTruthTable();
    const uint32_t arity = tbl.size();

    // safety: arity should match children count; skip otherwise
    if (node->children.size() != arity) continue;

    // 1) Unary identity removal:
    // If arity==1 and table behaves as identity on its single boolean input:
    // table.bits()[0] == 0 and table.bits()[1] == 1  (i.e., bit(0)==0, bit(1)==1)
    if (arity == 1) {
      bool b0 = tbl.bits().bit(0);
      bool b1 = tbl.bits().bit(1);
      if (b0 == false && b1 == true) {
        // identity: replace node by its single child
        auto child = node->children[0];
        replaceChildInParent(node_shared, child);
        // parent ownership updated by replaceChildInParent; no further action needed.
        continue;
      }
      // Note: if unary negation (b0==1 && b1==0) we'd need a NOT node representation;
      // skipping that because tree does not have unary NOT wrapper type.
    }

    // 2) If all children pointers are identical (all reference same subnode),
    //    try to reduce arity by reasoning on truth table as a function of single input.
    bool all_equal = true;
    for (size_t i = 1; i < node->children.size(); ++i) {
      if (node->children[i].get() != node->children[0].get()) {
        all_equal = false; break;
      }
    }
    if (all_equal && arity >= 1) {
      // When all children equal X, the table reduces to a unary function f(x) where:
      // f(0) = table.bits().bit(0), f(1) = table.bits().bit((1<<arity)-1)
      // If both outputs equal 0 => node is constant FALSE (cannot represent directly)
      // If both outputs equal 1 => node is constant TRUE (cannot represent directly)
      // If f is identity (0->0,1->1) we can replace node by child.
      // If f is negation (0->1,1->0) we cannot represent without NOT wrapper, skip.
      uint32_t idx0 = 0;
      uint32_t idx1 = (1u << arity) - 1u;
      bool out0 = tbl.bits().bit(idx0);
      bool out1 = tbl.bits().bit(idx1);
      if (!out0 && out1) {
        // behaves as identity => replace with child
        auto child = node->children[0];
        replaceChildInParent(node_shared, child);
        continue;
      }
      // if constant both 0 or both 1, we cannot create a constant leaf in this tree structure,
      // so conservatively leave node untouched.
    }

    // 3) Duplicate-child elimination when possible:
    //    If multiple children are identical and removing duplicates does not change table semantics,
    //    perform collapse. This is conservative: attempt when removing duplicate columns yields
    //    a table whose value for all remaining input combinations is consistent.
    //    Implementation: build map from child pointer to vector of positions,
    //    and if there exists at least one repeated child, test collapsing by remapping table.
    //    Restrict collapse to small arities to keep it cheap.
    if (arity > 1 && arity <= 8) { // limit complexity
      // group positions by identical child pointer
      std::unordered_map<Node*, std::vector<size_t>> posmap;
      posmap.reserve(node->children.size()*2);
      for (size_t i = 0; i < node->children.size(); ++i)
        posmap[node->children[i].get()].push_back(i);

      bool anyDup = false;
      for (auto &ent : posmap) if (ent.second.size() > 1) { anyDup = true; break; }
      if (anyDup) {
        // Attempt collapse: create canonical ordering of unique children
        std::vector<Node*> uniqueChildren;
        std::vector<std::vector<size_t>> positions;
        uniqueChildren.reserve(posmap.size());
        positions.reserve(posmap.size());
        for (auto &p : posmap) {
          uniqueChildren.push_back(p.first);
          positions.push_back(p.second);
        }

        size_t newArity = uniqueChildren.size();
        // Build map from original index -> reduced index
        std::vector<size_t> origToReduced(arity);
        for (size_t r = 0; r < positions.size(); ++r)
          for (size_t originalPos : positions[r])
            origToReduced[originalPos] = r;

        // Check that for every reduced input vector, all original input vectors that map to it
        // produce identical table output.
        bool collapseSafe = true;
        size_t origRows = 1u << arity;
        for (size_t origMask = 0; origMask < origRows && collapseSafe; ++origMask) {
          // compute reduced mask
          size_t redMask = 0;
          for (size_t bit = 0; bit < arity; ++bit) {
            if (origMask & (1u << bit))
              redMask |= (1u << origToReduced[bit]);
          }
          // collect all original masks that map to same redMask; verify their outputs equal
          bool expected = tbl.bits().bit((uint32_t)origMask);
          // For efficiency, test only other masks with same redMask by flipping duplicated positions
          // Approach: iterate over all origMask2 and compare when mapping equal
          for (size_t origMask2 = 0; origMask2 < origRows; ++origMask2) {
            // quick skip: if mapping differs, continue
            size_t redMask2 = 0;
            for (size_t bit = 0; bit < arity; ++bit)
              if (origMask2 & (1u << bit))
                redMask2 |= (1u << origToReduced[bit]);
            if (redMask2 != redMask) continue;
            if (tbl.bits().bit((uint32_t)origMask2) != expected) {
              collapseSafe = false;
              break;
            }
          }
        }

        if (collapseSafe && newArity < arity && newArity >= 1) {
          // Construct new node with uniqueChildren as children and reattach in place
          // Create a temporary new Node sharing the same dnlid/termid, then replace
          auto newNode = std::make_shared<Node>(this, node->dnlid, node->termid);
          // attach unique children in same order as uniqueChildren vector
          for (Node* cptr : uniqueChildren)
            newNode->addChild(std::shared_ptr<Node>(cptr->shared_from_this()));
          // Reattach into parent
          replaceChildInParent(node_shared, newNode);
          // Note: numExternalInputs_ may change if we removed Input leaves; we'll recompute below
          continue;
        }
      }
    }

    // No local simplification applied for this node
  } // end for order

  // Recompute border leaves and numExternalInputs_ conservatively.
  // We cannot accurately track numExternalInputs_ deltas in all cases during transformations above,
  // so recalc from root.
  // Rebuild numExternalInputs_ by scanning leaves for Type::Input and finding max index +1
  size_t maxInput = 0;
  bool anyInput = false;
  struct F2 { Node* n; };
  std::vector<F2, tbb::tbb_allocator<F2>> stk;
  if (root_) stk.push_back({root_.get()});
  while (!stk.empty()) {
    auto f = stk.back(); stk.pop_back();
    if (f.n->type == Node::Type::Table) {
      for (auto &c : f.n->children) stk.push_back({c.get()});
    } else if (f.n->type == Node::Type::Input) {
      anyInput = true;
      if (f.n->inputIndex > maxInput) maxInput = f.n->inputIndex;
    } else if (f.n->type == Node::Type::P) {
      // P nodes contain one Input child originally, nothing to do
      for (auto &c : f.n->children) stk.push_back({c.get()});
    }
  }
  if (anyInput) numExternalInputs_ = maxInput + 1;
  else numExternalInputs_ = 0;

  // update borders from new tree structure
  updateBorderLeaves();
}

