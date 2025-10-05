// MiterStrategyWithGlucose.cpp

#include "MiterStrategy.h"
#include "BoolExpr.h"
#include "BuildPrimaryOutputClauses.h"
#include "NLUniverse.h"
#include "SNLDesignModeling.h"
#include "SNLDesignModeling.h"
#include "SNLLogicCloud.h"
#include "SNLTruthTable2BoolExpr.h"

// include Glucose headers (adjust path to your checkout)
#include "core/Solver.h"
#include "simp/SimpSolver.h"

#include <string>
#include <unordered_map>
#include "NetlistGraph.h"
#include "SNLEquipotential.h"
#include "SNLLogicCone.h"
#include "SNLPath.h"

// For executeCommand
#include <cstdlib>
#include <stack>

using namespace naja;
using namespace naja::NL;
using namespace KEPLER_FORMAL;

namespace {

void executeCommand(const std::string& command) {
  int result = system(command.c_str());
  if (result != 0) {
    std::cerr << "Command execution failed." << std::endl;
  }
}
//
// A tiny Tseitin‐translator from BoolExpr→Glucose CNF.
//
// Returns a Glucose::Lit that stands for `e`, and adds
// all necessary clauses to S so that Lit ↔ (e) holds.
//
// node2var      caches each subformula’s fresh variable index.
// varName2idx   coalesces all inputs of the same name to one var.
//

// Assume BoolExpr, Op, getLeft(), getRight(), getOp(), getName() are defined 
// as in your original code.

Glucose::Lit tseitinEncode(
    Glucose::SimpSolver& S,
    BoolExpr* root,
    std::unordered_map<const BoolExpr*, int>& node2var,
    std::unordered_map<std::string, int>& varName2idx)
{
    auto getOrCreateVar = [&](const std::string& key) -> int {
        auto it = varName2idx.find(key);
        if (it != varName2idx.end()) return it->second;
        int v = S.newVar();
        varName2idx[key] = v;
        return v;
    };

    auto constLit = [&](bool value) -> Glucose::Lit {
        const std::string key = value ? "$__CONST_TRUE__" : "$__CONST_FALSE__";
        int v = getOrCreateVar(key);
        Glucose::Lit lv = Glucose::mkLit(v);
        S.addClause(value ? lv : ~lv);
        return lv;
    };

    struct Frame {
        BoolExpr* expr;
        bool visited = false;
        Glucose::Lit leftLit, rightLit;
    };

    std::stack<Frame> stk;
    stk.push({root, false, {}, {}});
    std::unordered_map<const BoolExpr*, Glucose::Lit> result;

    while (!stk.empty()) {
        Frame& fr = stk.top();
        BoolExpr* e = fr.expr;

        // If already encoded, reuse
        if (node2var.count(e)) {
            result[e] = Glucose::mkLit(node2var[e]);
            stk.pop();
            continue;
        }

        // Leaf VAR or CONST
        if (!fr.visited && e->getOp() == Op::VAR) {
            const std::string& name = e->getName();
            Glucose::Lit lit;
            if (name == "0" || name == "false"  || name == "False" || name == "FALSE")
                lit = constLit(false);
            else if (name == "1" || name == "true" || name == "True" || name == "TRUE")
                lit = constLit(true);
            else {
                int v = getOrCreateVar(name);
                lit = Glucose::mkLit(v);
            }
            node2var[e] = Glucose::var(lit);
            result[e] = lit;
            stk.pop();
            continue;
        }

        // First time we see this node, push children
        if (!fr.visited) {
            fr.visited = true;
            if (e->getRight()) stk.push({e->getRight(), false, {}, {}});
            if (e->getLeft())  stk.push({e->getLeft(),  false, {}, {}});
            continue;
        }

        // Children have been processed; retrieve their lits
        if (e->getLeft())  fr.leftLit  = result[e->getLeft()];
        if (e->getRight()) fr.rightLit = result[e->getRight()];

        // Create fresh var for this gate
        int v = S.newVar();
        Glucose::Lit lit_v = Glucose::mkLit(v);
        node2var[e] = v;
        result[e] = lit_v;

        // Emit Tseitin clauses
        switch (e->getOp()) {
            case Op::NOT:
                S.addClause(~lit_v, ~fr.leftLit);
                S.addClause( lit_v,  fr.leftLit);
                break;
            case Op::AND:
                S.addClause(~lit_v,  fr.leftLit);
                S.addClause(~lit_v,  fr.rightLit);
                S.addClause( lit_v, ~fr.leftLit, ~fr.rightLit);
                break;
            case Op::OR:
                S.addClause(~fr.leftLit, lit_v);
                S.addClause(~fr.rightLit, lit_v);
                S.addClause(~lit_v, fr.leftLit, fr.rightLit);
                break;
            case Op::XOR:
                S.addClause(~lit_v, ~fr.leftLit, ~fr.rightLit);
                S.addClause(~lit_v,  fr.leftLit,  fr.rightLit);
                S.addClause( lit_v, ~fr.leftLit,  fr.rightLit);
                S.addClause( lit_v,  fr.leftLit, ~fr.rightLit);
                break;
            default:
                // Handle other operators if needed
                break;
        }

        stk.pop();
    }

    return result.at(root);
}

}  // namespace

void MiterStrategy::normalizeInputs(std::vector<naja::DNL::DNLID>& inputs0,
                       std::vector<naja::DNL::DNLID>& inputs1,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& inputs0Map,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& inputs1Map) {
    // find the intersection of inputs0 and inputs1 based on the getFullPathIDs of DNLTerminal and the diffs
    std::set<std::vector<NLID::DesignObjectID>> paths0;
    std::set<std::vector<NLID::DesignObjectID>> paths1;
    std::set<std::vector<NLID::DesignObjectID>> pathsCommon;
    
    for (const auto& [path0, input0] : inputs0Map) {
      paths0.insert(path0);
    }
    for (const auto& [path1, input1] : inputs1Map) {
      paths1.insert(path1);
    }
    size_t index = 0;
    for (const auto& [path0, input0] : inputs0Map) {
      //printf("input %zu from %lu\n", index++, inputs0Map.size());
      if (paths1.find(path0) != paths1.end()) {
        pathsCommon.insert(path0);
      }
    }
    std::vector<naja::DNL::DNLID> diff0;
    for (const auto& [path0, input0] : inputs0Map) {
      if (pathsCommon.find(path0) == pathsCommon.end()) {
        diff0.push_back(input0);
        auto pathInstance = path0;
        pathInstance.pop_back(); // removing bit ID
        pathInstance.pop_back(); // removing terminal ID
        naja::NL::SNLPath p = naja::NL::SNLPath(naja::DNL::get()->getTop().getSNLModel(), pathInstance);
        printf("diff0 input: %s\n", p.getString().c_str());
      }
    }
    std::vector<naja::DNL::DNLID> diff1;
    for (const auto& [path1, input1] : inputs1Map) {
      if (pathsCommon.find(path1) == pathsCommon.end()) {
        diff1.push_back(input1);
        auto pathInstance = path1;
        pathInstance.pop_back(); // removing bit ID
        pathInstance.pop_back(); // removing terminal ID
        naja::NL::SNLPath p = naja::NL::SNLPath(naja::DNL::get()->getTop().getSNLModel(), pathInstance);
        printf("diff1 input: %s\n", p.getString().c_str());
      }
    }
    inputs0.clear();
    for (const auto& path : pathsCommon) {
      inputs0.push_back(inputs0Map.at(path));
    }
    inputs0.insert(inputs0.end(), diff0.begin(), diff0.end());
    inputs1.clear();
    for (const auto& path : pathsCommon) {
      inputs1.push_back(inputs1Map.at(path));
    }
    inputs1.insert(inputs1.end(), diff1.begin(), diff1.end());
    printf("size of common inputs: %lu\n", pathsCommon.size());
    printf("size of diff0 inputs: %lu\n", diff0.size());
    printf("size of diff1 inputs: %lu\n", diff1.size());
}

void MiterStrategy::normalizeOutputs(std::vector<naja::DNL::DNLID>& outputs0,
                        std::vector<naja::DNL::DNLID>& outputs1,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& outputs0Map,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& outputs1Map) {
    // find the intersection of outputs0 and outputs1 based on the getFullPathIDs of DNLTerminal and the diffs
    std::set<std::vector<NLID::DesignObjectID>> paths1;
    std::set<std::vector<NLID::DesignObjectID>> pathsCommon;
    for (const auto& [path1, output1] : outputs1Map) {
      paths1.insert(path1);
    }
    size_t index = 0;
    for (const auto& [path0, output0] : outputs0Map) {
      //printf("output %zu from %lu\n", index++, outputs0Map.size());
      if (paths1.find(path0) != paths1.end()) {
        pathsCommon.insert(path0);
      }
    }
    std::vector<naja::DNL::DNLID> diff0;
    for (const auto& [path0, output0] : outputs0Map) {
      if (pathsCommon.find(path0) == pathsCommon.end()) {
        diff0.push_back(output0);
      }
    }
    std::vector<naja::DNL::DNLID> diff1;
    for (const auto& [path1, output1] : outputs1Map) {
      if (pathsCommon.find(path1) == pathsCommon.end()) {
        diff1.push_back(output1);
      }
    }
    outputs0.clear();
    for (const auto& path : pathsCommon) {
      outputs0.push_back(outputs0Map.at(path));
    }
    outputs0.insert(outputs0.end(), diff0.begin(), diff0.end());
    outputs1.clear();
    for (const auto& path : pathsCommon) {
      outputs1.push_back(outputs1Map.at(path));
    }
    outputs1.insert(outputs1.end(), diff1.begin(), diff1.end());
    printf("size of common outputs: %lu\n", pathsCommon.size());
    printf("size of diff0 outputs: %lu\n", diff0.size());
    printf("size of diff1 outputs: %lu\n", diff1.size());
    if (outputs0.size() == outputs1.size()) {
      if (outputs0 != outputs1) {
        // build the paths vector for outputs0 and outputs1
        std::vector<std::vector<NLID::DesignObjectID>> paths0;
        std::vector<std::vector<NLID::DesignObjectID>> paths1;
        for (const auto& output0 : outputs0) {
          std::vector<NLID::DesignObjectID> path;
          for (const auto& [path0, output0m] : outputs0Map) {
            if (output0m == output0) {
              path = path0;
              break;
            }
          }
          paths0.push_back(path);
        }
        for (const auto& output1 : outputs1) {
          std::vector<NLID::DesignObjectID> path;
          for (const auto& [path1, output1m] : outputs1Map) {
            if (output1m == output1) {
              path = path1;
              break;
            }
          }
          paths1.push_back(path);
        }
        if (paths0 != paths1) {
          printf("Miter outputs must match in order:\n");
          for (size_t i = 0; i < paths0.size(); ++i) {
            naja::NL::SNLPath p0 = naja::NL::SNLPath(naja::DNL::get()->getTop().getSNLModel(), paths0[i]);
            naja::NL::SNLPath p1 = naja::NL::SNLPath(naja::DNL::get()->getTop().getSNLModel(), paths1[i]);
            throw std::runtime_error(
                "Output " + std::to_string(i) + " mismatch: " + p0.getString() + " vs " + p1.getString());

          }
          assert(false && "Miter outputs must match in order");
        }
      }
    }
}

bool MiterStrategy::run() {
  // build both sets of POs
  topInit_ = NLUniverse::get()->getTopDesign();
  NLUniverse* univ = NLUniverse::get();
  naja::DNL::destroy();
  univ->setTopDesign(top0_);
  BuildPrimaryOutputClauses builder0;
  builder0.collect();
  naja::DNL::destroy();
  univ->setTopDesign(top1_);
  BuildPrimaryOutputClauses builder1;
  builder1.collect();

  // normalize inputs and outputs
  auto inputs0sort = builder0.getInputs();
  auto inputs1sort = builder1.getInputs();
  auto outputs0sort = builder0.getOutputs();
  auto outputs1sort = builder1.getOutputs();
  printf("size of inputs0: %lu\n", inputs0sort.size());
  printf("size of inputs1: %lu\n", inputs1sort.size());
  printf("size of outputs0: %lu\n", outputs0sort.size());
  printf("size of outputs1: %lu\n", outputs1sort.size());
  normalizeInputs(inputs0sort, inputs1sort, builder0.getInputsMap(), builder1.getInputsMap());
  normalizeOutputs(outputs0sort, outputs1sort, builder0.getOutputsMap(), builder1.getOutputsMap());
  //return false;
  builder0.setInputs(inputs0sort);
  builder1.setInputs(inputs1sort);
  builder0.setOutputs(outputs0sort);
  builder1.setOutputs(outputs1sort);
  naja::DNL::destroy();
  univ->setTopDesign(top0_);
  builder0.build();
  const auto& PIs0 = builder0.getInputs();
  const auto& POs0 = builder0.getPOs();
  auto outputs0 = builder0.getOutputs();
  auto inputs2inputsIDs0 = builder0.getInputs2InputsIDs();
  auto outputs2outputsIDs0 = builder0.getOutputs2OutputsIDs();
  naja::DNL::destroy();
  univ->setTopDesign(top1_);
  builder1.build();
  const auto& PIs1 = builder1.getInputs();
  const auto& POs1 = builder1.getPOs();
  auto outputs1 = builder1.getOutputs();
  auto inputs2inputsIDs1 = builder1.getInputs2InputsIDs();
  auto outputs2outputsIDs1 = builder1.getOutputs2OutputsIDs();

  std::vector<naja::DNL::DNLID> outputs2DnlIds = builder1.getOutputs();

  if (topInit_ != nullptr) {
    univ->setTopDesign(topInit_);
  }

  if (POs0.empty() || POs1.empty()) {
    return false;
  }
  // if (outputs0 != outputs1) {
  //   for (size_t i = 0; i < outputs0.size(); ++i) {
  //     printf("Output %zu: %lu vs %lu\n", i, outputs0[i], outputs1[i]);
  //   }
  //   assert(outputs0 == outputs1);
  // }
  // print POs0 and POs1 in 2 separate for loops
  // for (size_t i = 0; i < POs0.size(); ++i) {
  //   printf("POs0[%zu]: %s\n", i, POs0[i]->toString().c_str());
  // }
  // for (size_t i = 0; i < POs1.size(); ++i) {
  //   printf("POs1[%zu]: %s\n", i, POs1[i]->toString().c_str());
  // }
  

  // Check if both sets inputs are the same
  // if (inputs2inputsIDs0.size() != inputs2inputsIDs1.size() ||
  //     outputs2outputsIDs0.size() != outputs2outputsIDs1.size()) {
  //   printf("Miter inputs/outputs must match in size: %zu vs %zu\n",
  //          inputs2inputsIDs0.size(), inputs2inputsIDs1.size());
  //   return false;
  // }

  // for (const auto& [id0, ids0] : inputs2inputsIDs0) {
  //   auto it = inputs2inputsIDs1.find(id0);
  //   if (it == inputs2inputsIDs1.end() || it->second.first != ids0.first ||
  //       it->second.second != ids0.second) {
  //     printf("Miter inputs must match in IDs: %lu\n", id0);
  //     return false;
  //   }
  // }

  // // Check if both sets outputs are the same
  // if (outputs2outputsIDs0.size() != outputs2outputsIDs1.size()) {
  //   printf("Miter outputs must match in size: %zu vs %zu\n",
  //          outputs2outputsIDs0.size(), outputs2outputsIDs1.size());
  //   return false;
  // }

  // for (const auto& [id0, ids0] : outputs2outputsIDs0) {
  //   auto it = outputs2outputsIDs1.find(id0);
  //   if (it == outputs2outputsIDs1.end() || it->second.first != ids0.first ||
  //       it->second.second != ids0.second) {
  //     printf("Miter outputs must match in IDs: %lu\n", id0);
  //     return false;
  //   }
  // }

  // build the Boolean‐miter expression
  auto miter = buildMiter(POs0, POs1);

  // Now SAT check via Glucose
  Glucose::SimpSolver solver;

  // mappings for Tseitin encoding
  std::unordered_map<const BoolExpr*, int> node2var;
  std::unordered_map<std::string, int> varName2idx;

  // Tseitin‐encode & get the literal for the root
  Glucose::Lit rootLit = tseitinEncode(solver, miter, node2var, varName2idx);

  // Assert root == true
  solver.addClause(rootLit);

  // solve with no assumptions
  printf("Started Glucose solving\n");
  bool sat = solver.solve();
  printf("Finished Glucose solving: %s\n", sat ? "SAT" : "UNSAT");

  if (sat) {
    printf("Miter failed: analyzing individual POs\n");
    for (size_t i = 0; i < POs0.size(); ++i) {
      tbb::concurrent_vector<BoolExpr*> singlePOs0S;
      singlePOs0S.push_back(POs0[i]);
      tbb::concurrent_vector<BoolExpr*> singlePOs1S;
      singlePOs1S.push_back(POs1[i]);
      auto singleMiter = buildMiter(singlePOs0S, singlePOs1S);
      
      std::unordered_map<const BoolExpr*, int> singleNode2var;
      std::unordered_map<std::string, int> singleVarName2idx;
      // Tseitin‐encode the single miter
      Glucose::SimpSolver singleSolver;
      Glucose::Lit singleRootLit = tseitinEncode(
          singleSolver, singleMiter, singleNode2var, singleVarName2idx);
      // print singleMiter
      //printf("Single miter for PO %zu: %s\n", i,
      //     singleMiter->toString().c_str());
      

      singleSolver.addClause(singleRootLit);
      if (singleSolver.solve()) {
        // print singleMitter
        // printf("------ PO failed for %zu mitter: %s\n", i,
        //      singleMiter->toString().c_str());


        // printf("------ PO failed for %zu: %s vs %s\n", i,
        //      singlePOs0S[0]->toString().c_str(),
        //      singlePOs1S[0]->toString().c_str());
        failedPOs_.push_back(outputs0[i]);
        // If UNSAT, print the single miter
        printf("Check failed for PO: %zu\n", i);
        std::vector<naja::NL::SNLDesign*> topModels;
        topModels.push_back(top0_);
        topModels.push_back(top1_);
        std::vector<std::vector<naja::DNL::DNLID>> PIs;
        PIs.push_back(PIs0);
        PIs.push_back(PIs1);
        naja::NL::SNLEquipotential::Terms terms0;
        naja::NL::SNLEquipotential::Terms terms1;
        naja::NL::SNLEquipotential::InstTermOccurrences insTerms0;
        naja::NL::SNLEquipotential::InstTermOccurrences insTerms1;
        for (size_t j = 0; j < topModels.size(); ++j) {
          DNL::destroy();
          NLUniverse::get()->setTopDesign(topModels[j]);
          SNLLogicCone cone(outputs0[i], PIs[j]);
          cone.run();
          std::string dotFileNameEquis(
              std::string(prefix_ + "_" + DNL::get()->getDNLTerminalFromID(outputs0[i]).getSnlBitTerm()->getName().getString() + std::to_string(outputs0[i]) + "_" +std::to_string(j) + std::string(".dot")));
          std::string svgFileNameEquis(
              std::string(prefix_ + "_" + DNL::get()->getDNLTerminalFromID(outputs0[i]).getSnlBitTerm()->getName().getString() + std::to_string(outputs0[i]) + "_" + std::to_string(j) + std::string(".svg")));
          SnlVisualiser snl2(topModels[j], cone.getEquipotentials());
          for (const auto& equi : cone.getEquipotentials()) {
            for (const auto& term : equi.getTerms()) {
              if (j == 0) {
                terms0.insert(term);
              } else {
                terms1.insert(term);
              }
            }
            for (const auto& termOcc : equi.getInstTermOccurrences()) {
              if (j == 0) {
                insTerms0.insert(termOcc);
              } else {
                insTerms1.insert(termOcc);
              }
            }
          }
          
          //snl2.process();
          //snl2.getNetlistGraph().dumpDotFile(dotFileNameEquis.c_str());
          //  executeCommand(std::string(std::string("dot -Tsvg ") +
          //                             dotFileNameEquis + std::string(" -o ") +
          //                             svgFileNameEquis)
          //                     .c_str());
          //printf("svg file name: %s\n", svgFileNameEquis.c_str());
        }
        // find intersection and diff of terms0 and terms1
        naja::NL::SNLEquipotential::Terms termsCommon;
        naja::NL::SNLEquipotential::Terms termsDiff;
        for (const auto& term0 : terms0) {
          bool found = false;
          for (const auto& term1 : terms1) {
            if (term0->getID() == term1->getID() && term0->getBit() == term1->getBit()) {
              found = true;
              break;
            }
          }
          if (found) {
            termsCommon.insert(term0);
          } else {
            termsDiff.insert(term0);
          }
        }
        for (const auto& term1 : terms1) {
          bool found = false;
          for (const auto& term0 : terms0) {
            if (term0->getID() == term1->getID() && term0->getBit() == term1->getBit()) {
              found = true;
              break;
            }
          }
          if (!found) {
            termsDiff.insert(term1);
          }
        }
        // print termsDiff
        for (const auto& term : termsDiff) {
          if (term->getDirection() == naja::NL::SNLBitTerm::Direction::Output) {
            continue;
          }
          printf("Diff term: %s\n",
                 term->getString().c_str());
        }
        // find intersection aprintf("size of intersend diff of insTerms0 and insTerms1
        naja::NL::SNLEquipotential::InstTermOccurrences insTermsCommon;
        naja::NL::SNLEquipotential::InstTermOccurrences insTermsDiff;
        for (const auto& term0 : insTerms0) {
          bool found = false;
          for (const auto& term1 : insTerms1) {
            if (term0.getPath().getPathIDs() == term1.getPath().getPathIDs() && 
                term0.getInstTerm()->getInstance()->getID() == term1.getInstTerm()->getInstance()->getID() &&
                term0.getInstTerm()->getBitTerm()->getID() == term1.getInstTerm()->getBitTerm()->getID() &&
                term0.getInstTerm()->getBitTerm()->getBit() == term1.getInstTerm()->getBitTerm()->getBit()) {
              found = true;
              break;
            }
          }
          if (found) {
            insTermsCommon.insert(term0);
          } else {
            insTermsDiff.insert(term0);
            if (term0.getInstTerm()->getDirection() == naja::NL::SNLInstTerm::Direction::Input || 
              !term0.getInstTerm()->getInstance()->getModel()->getInstances().empty()) {
              continue;
            }
            printf("Diff 0 inst term %s with direction %s\n",
                  term0.getString().c_str(), term0.getInstTerm()->getDirection().getString().c_str());
          }
        }
        for (const auto& term1 : insTerms1) {
          bool found = false;
          for (const auto& term0 : insTerms0) {
            if (term0.getPath().getPathIDs() == term1.getPath().getPathIDs() && 
                term0.getInstTerm()->getInstance()->getID() == term1.getInstTerm()->getInstance()->getID() &&
                term0.getInstTerm()->getBitTerm()->getID() == term1.getInstTerm()->getBitTerm()->getID() &&
                term0.getInstTerm()->getBitTerm()->getBit() == term1.getInstTerm()->getBitTerm()->getBit()) {
              found = true;
              break;
            }
          }
          if (!found) {
            insTermsDiff.insert(term1);
            if (term1.getInstTerm()->getDirection() == naja::NL::SNLInstTerm::Direction::Input || 
              !term1.getInstTerm()->getInstance()->getModel()->getInstances().empty()) {
              continue;
            }
            printf("Diff 1 inst term %s with direction %s\n",
                  term1.getString().c_str(), term1.getInstTerm()->getDirection().getString().c_str());
          }
        }
        // print insTermsDiff
        // for (const auto& term : insTermsDiff) {
        //   // static cast to SNLInstTerm
        //   //naja::NL::SNLInstTerm* instTerm = static_cast<naja::NL::SNLInstTerm*>(term.getInstTerm());
        //   if (term.getInstTerm()->getDirection() == naja::NL::SNLInstTerm::Direction::Input || 
        //     !term.getInstTerm()->getInstance()->getModel()->getInstances().empty()) {
        //     continue;
        //   }
        //   printf("Diff inst term %s\n",
        //          term.getString().c_str());
        // }
        printf("size of intersection of terms: %zu\n", termsCommon.size());
        printf("size of diff of terms: %zu\n", termsDiff.size());
        printf("size of intersection of inst terms: %zu\n", insTermsCommon.size());
        printf("size of diff of inst terms: %zu\n", insTermsDiff.size());
      }
    }
  }
  if (topInit_ != nullptr) {
    univ->setTopDesign(topInit_);
  }
  // if UNSAT → miter can never be true → outputs identical
  printf("Circuits are %s\n", sat ? "DIFFERENT" : "IDENTICAL");
  return !sat;
}

BoolExpr* MiterStrategy::buildMiter(
    const tbb::concurrent_vector<BoolExpr*>& A,
    const tbb::concurrent_vector<BoolExpr*>& B) const {
  // if (A.size() != B.size()) {
  //   printf("Miter inputs must match in length: %zu vs %zu\n", A.size(),
  //          B.size());
  //   assert(false && "Miter inputs must match in length");
  // }

  // Empty miter = always‐false (no outputs to compare)
  if (A.empty()) {
    assert(false);
    return BoolExpr::createFalse();
  }
  //printf("A[0]: %s\n", A[0]->toString().c_str());
  //printf("B[0]: %s\n", B[0]->toString().c_str());
  // Start with the first XOR
  auto miter = BoolExpr::Xor(A[0], B[0]);
  //printf("Initial miter: %s\n", miter->toString().c_str());
  // OR in the rest
  for (size_t i = 1; i < A.size(); ++i) {
    if (B.size() <= i) {
      //printf("Miter different number of outputs: %zu vs %zu\n", A.size(),
      //       B.size());
      break;
    }
    auto diff = BoolExpr::Xor(A[i], B[i]);
    miter = BoolExpr::Or(miter, diff);
  }
  return miter;
}
