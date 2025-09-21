#include "SNLTruthTableTree.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>

using namespace naja::NL;
using namespace KEPLER_FORMAL;

//----------------------------------------------------------------------
// Node::addChild
//----------------------------------------------------------------------
void SNLTruthTableTree::Node::addChild(std::shared_ptr<Node> child) {
  // cycle detection walking parent chain
  auto self = shared_from_this();
  for (auto node = self; node; node = node->parent.lock()) {
    if (node->type == Type::Table
        && node->dnlid == child->dnlid
        && node->termid == child->termid
        && child->dnlid != naja::DNL::DNLID_MAX
        && child->termid != naja::DNL::DNLID_MAX) {
      throw std::invalid_argument("addChild: cycle detected");
    }
  }

  // attach
  children.push_back(std::move(child));
  children.back()->parent = self;
}

//----------------------------------------------------------------------
// Node::getTruthTable
//----------------------------------------------------------------------
SNLTruthTable SNLTruthTableTree::Node::getTruthTable() const {
  if (type == Type::Table) {
    return SNLDesignModeling::getTruthTable(
      naja::DNL::get()->getDNLInstanceFromID(dnlid).getSNLModel(),
      naja::DNL::get()->getDNLTerminalFromID(termid)
                   .getSnlBitTerm()->getOrderID());
  }
  else if (type == Type::P) {
    return SNLTruthTable(1,2);
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

SNLTruthTableTree::SNLTruthTableTree(Node::Type type) {
  if (type != Node::Type::P)
    throw std::invalid_argument("invalid type ctor");

  root_ = std::make_shared<Node>(this);
  root_->addChild(std::make_shared<Node>(0, this));
  numExternalInputs_ = 1;
}

SNLTruthTableTree::SNLTruthTableTree(naja::DNL::DNLID instid,
                                     naja::DNL::DNLID termid)
{
  auto table = SNLDesignModeling::getTruthTable(
    naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel(),
    naja::DNL::get()->getDNLTerminalFromID(termid)
                 .getSnlBitTerm()->getOrderID());

  auto arity = table.size();
  root_ = std::make_shared<Node>(this, instid, termid);
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
  std::vector<Frame> stk;
  stk.reserve(64);
  stk.push_back({nullptr,0,root_.get()});

  while (!stk.empty()) {
    auto f = stk.back(); stk.pop_back();
    if (f.node->type == Node::Type::Table) {
      for (size_t i = f.node->children.size(); i-- > 0;) {
        stk.push_back({f.node,i,f.node->children[i].get()});
      }
    }
    else {
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
  // locate slot
  size_t idx = borderIndex;
  if (idx >= borderLeaves_.size()) {
    bool found=false;
    for (size_t i=0;i<borderLeaves_.size();++i)
      if (borderLeaves_[i].extIndex == borderIndex) {
        idx=i; found=true; break;
      }
    if (!found) throw std::out_of_range("concat: leafIndex out of range");
  }
  auto leaf = borderLeaves_[idx];

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
  if (instid!=naja::DNL::DNLID_MAX && termid!=naja::DNL::DNLID_MAX) {
    auto tbl = SNLDesignModeling::getTruthTable(
      naja::DNL::get()->getDNLInstanceFromID(instid).getSNLModel(),
      naja::DNL::get()->getDNLTerminalFromID(termid)
                   .getSnlBitTerm()->getOrderID());
    arity = tbl.size();
    newNode = std::make_shared<Node>(this, instid, termid);
  }
  else {
    auto tbl = SNLTruthTable(1,2);
    arity = tbl.size();
    newNode = std::make_shared<Node>(this);
  }

  // re-use old leaf as first child, then fresh inputs
  if (arity>0) {
    newNode->addChild(std::move(oldLeaf));
    for (uint32_t i=1; i<arity; ++i)
      newNode->addChild(
        std::make_shared<Node>(numExternalInputs_ + (i-1), this)
      );
  }

  // reattach
  if (leaf.parent)
    leaf.parent->children[leaf.childPos] = newNode;
  else
    root_ = newNode;

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
  const std::vector<std::pair<naja::DNL::DNLID,
                              naja::DNL::DNLID>>& tables)
{
  int newInputs = (int)numExternalInputs_;
  if (tables.size() > borderLeaves_.size())
    throw std::invalid_argument("too many tables in concatFull");

  for (size_t i=0;i<tables.size();++i) {
    auto const& n = concatBody(i, tables[i].first, tables[i].second);
    newInputs += (n.getTruthTable().size() - 1);
  }
  numExternalInputs_ = (size_t)newInputs;
  updateBorderLeaves();
}

//----------------------------------------------------------------------
// isInitialized / print
//----------------------------------------------------------------------
bool SNLTruthTableTree::isInitialized() const {
  if (!root_) return false;
  struct Frame{Node* n;};
  std::vector<Frame> stk{{root_.get()}};
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
  std::vector<Frame> stk{{root_.get(), 0}};
  while(!stk.empty()) {
    auto [n,d] = stk.back(); stk.pop_back();
    // print indentation + node info here…
    if (n->type == Node::Type::Table) {
      for (size_t i=n->children.size(); i-->0;)
        stk.push_back({n->children[i].get(), d+1});
    }
  }
}
