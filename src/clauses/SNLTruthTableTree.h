#ifndef SNLTRUTHTABLETREE_H
#define SNLTRUTHTABLETREE_H

#include "SNLTruthTable.h"
#include <vector>
#include <memory>
#include <cstddef>
#include "DNL.h"
#include "SNLDesignModeling.h"

namespace KEPLER_FORMAL {

class SNLTruthTableTree {
public:
  struct Node 
    : public std::enable_shared_from_this<Node>   // enable shared_from_this
  {
    enum class Type { Input/* intermediate input of table in tree, used for concat*/, 
    Table/*table node*/, P/*primary input note, all external nodes should P in the end of cloud creation*/ } type;

    // for Input
    size_t inputIndex = (size_t)-1;

    // debug
    size_t nodeID = 0;       
    SNLTruthTableTree* tree = nullptr;

    // for Table nodes
    naja::DNL::DNLID dnlid    = naja::DNL::DNLID_MAX;
    naja::DNL::DNLID termid   = naja::DNL::DNLID_MAX;
    std::vector<std::shared_ptr<Node>> children;

    // parent pointer as weak_ptr to break cycles
    std::weak_ptr<Node> parent;     

    //--- ctors
    // explicit Node(SNLTruthTableTree* t, naja::DNL::DNLID i, 
    //      naja::DNL::DNLID term)                   // Type::P
    //   : type(Type::P), tree(t), dnlid(i), termid(term) {
    //   nodeID = tree->lastID_++;
    // }

    explicit Node(size_t idx, SNLTruthTableTree* t)       // Type::Input
      : type(Type::Input), inputIndex(idx), tree(t) {
      nodeID = tree->lastID_++;
    }

    Node(SNLTruthTableTree* t, 
         naja::DNL::DNLID i, 
         naja::DNL::DNLID term, Type type = Type::Table)                          // Type::Table
      : type(type), tree(t), dnlid(i), termid(term) 
    {
      //assert(i != naja::DNL::DNLID_MAX && term != naja::DNL::DNLID_MAX);
      nodeID = tree->lastID_++;
    }

    Node(const Node& other) = delete;

    // evaluate recursively
    bool eval(const std::vector<bool>& extInputs) const;

    // add a child, detect cycles, set child's parent
    void addChild(std::shared_ptr<Node> child);

    // get the table, only valid for Table and P nodes
    const SNLTruthTable& getTruthTable() const;
  };

  //--- public API
  SNLTruthTableTree();
  //SNLTruthTableTree(Node::Type type);
  SNLTruthTableTree(naja::DNL::DNLID instid, naja::DNL::DNLID termid, Node::Type type = Node::Type::Table);

  size_t size() const;
  bool eval(const std::vector<bool>& extInputs) const;

  void concat(size_t borderIndex,
              naja::DNL::DNLID instid,
              naja::DNL::DNLID termid);

  void concatFull(const std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
            tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>& tables);

  const Node* getRoot() const { return root_.get(); }
  bool isInitialized() const;
  void print() const;

private:
  struct BorderLeaf {
    Node*  parent;    // nullptr if root
    size_t childPos;  
    size_t extIndex;  
  };

  void updateBorderLeaves();
  const Node& concatBody(size_t borderIndex,
                         naja::DNL::DNLID instid,
                         naja::DNL::DNLID termid);

  std::shared_ptr<Node>   root_;
  size_t                  numExternalInputs_ = 0;
  std::vector<BorderLeaf> borderLeaves_;
  size_t                  lastID_ = 2;       // for debug
  static const SNLTruthTable      PtableHolder_; // dummy table for P nodes
};

} // namespace KEPLER_FORMAL

#endif // SNLTRUTHTABLETREE_H
