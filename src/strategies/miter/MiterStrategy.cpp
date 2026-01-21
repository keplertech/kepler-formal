// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "MiterStrategy.h"
#include "BoolExpr.h"
#include "BuildPrimaryOutputClauses.h"
#include "NLUniverse.h"
#include "SNLDesignModeling.h"
#include "SNLLogicCloud.h"

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

// spdlog
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_sinks.h>  // ensure console sink is available
#include <spdlog/spdlog.h>

//#define DEBUG_CHECKS

using namespace naja;
using namespace naja::NL;
using namespace KEPLER_FORMAL;

SNLDesign* MiterStrategy::top0_ = nullptr;
SNLDesign* MiterStrategy::top1_ = nullptr;
std::string MiterStrategy::logFileName_ = "";
namespace {

static std::shared_ptr<spdlog::logger> logger;

void ensureLoggerInitialized() {
  if (logger) return;

  try {
    // 1) Choose a default file name in the current working directory
    int logIndex = 0;
    while (true) {
      std::string candidate = "miter_log_" + std::to_string(logIndex) + ".txt";
      std::ifstream infile(candidate);
      if (infile.good()) {
        ++logIndex;
      } else {
        break;
      }
    }
    std::string chosenLogFile = "miter_log_" + std::to_string(logIndex) + ".txt";

    // 2) If user provided an explicit path, try to use it (with safe checks)
    if (!MiterStrategy::logFileName_.empty()) {
      std::filesystem::path p(MiterStrategy::logFileName_);
      auto parent = p.parent_path();

      // If parent is empty, treat the provided name as a filename in CWD
      if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
        if (ec) {
          // LCOV_EXCL_START
          // Failed to create requested directory; log and fall back
          std::cerr << "Warning: failed to create log directory '" << parent.string()
                    << "': " << ec.message() << " (" << ec.value() << "). Using fallback path.\n";
          // LCOV_EXCL_STOP
        } else {
          chosenLogFile = p.string();
        }
      } else {
        // Provided path had no parent; use it as-is
        chosenLogFile = p.string();
      }
    }

    // 3) Try to create file sink; if it fails, fall back to temp dir or stdout
    try {
      auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(chosenLogFile, true);
      logger = std::make_shared<spdlog::logger>("miter_logger", file_sink);
    } catch (const spdlog::spdlog_ex& ex) {
      // LCOV_EXCL_START
      // Try a safe fallback: temp directory
      std::error_code ec;
      auto tmp = std::filesystem::temp_directory_path(ec);
      if (!ec) {
        std::filesystem::path fallback = tmp / ("miter_log_fallback_" + std::to_string(::getpid()) + ".txt");
        try {
          auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(fallback.string(), true);
          logger = std::make_shared<spdlog::logger>("miter_logger", file_sink);
        } catch (...) {
          // Final fallback to stdout sink
          auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
          logger = std::make_shared<spdlog::logger>("miter_logger_fallback", console_sink);
          logger->set_level(spdlog::level::debug);
          spdlog::register_logger(logger);
          logger->error("spdlog initialization failed for file sink and temp fallback: {}", ex.what());
        }
      } else {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
        logger = std::make_shared<spdlog::logger>("miter_logger_fallback", console_sink);
        logger->set_level(spdlog::level::debug);
        spdlog::register_logger(logger);
        logger->error("spdlog initialization failed and temp_directory_path() failed: {}", ex.what());
      }
      // LCOV_EXCL_STOP
    }

    // 4) Finalize logger if created
    if (logger) {
      logger->set_level(spdlog::level::info);
      logger->flush_on(spdlog::level::info);
      spdlog::register_logger(logger);
    }
  } catch (const std::exception& ex) {
    // LCOV_EXCL_START
    // Last-resort fallback to stdout logger to avoid crashing tests
    auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
    logger = std::make_shared<spdlog::logger>("miter_logger_fallback", console_sink);
    logger->set_level(spdlog::level::debug);
    spdlog::register_logger(logger);
    logger->error("Unexpected exception initializing logger: {}", ex.what());
    // LCOV_EXCL_STOP
  }
}


// void executeCommand(const std::string& command) {
//   ensureLoggerInitialized();
//   int result = system(command.c_str());
//   if (result != 0) {
//     logger->error("Command execution failed: {} (exit code {})", command,
//                   result);
//   } else {
//     logger->debug("Command executed successfully: {}", command);
//   }
// }

//
// A tiny Tseitin-translator from BoolExpr -> Glucose CNF.
//
// Returns a Glucose::Lit that stands for `e`, and adds
// all necessary clauses to S so that Lit ↔ (e) holds.
//
// node2var      caches each subformula’s fresh variable index.
// varName2idx   coalesces all inputs of the same name to one var.
//

Glucose::Lit tseitinEncode(
    Glucose::SimpSolver& S,
    BoolExpr* root,
    std::unordered_map<BoolExpr*, int>& node2var,
    std::unordered_map<std::string, int>& varName2idx) {
  ensureLoggerInitialized();
  logger->debug("Starting Tseitin encode for root expr");

  auto getOrCreateVar = [&](const std::string& key) -> int {
    auto it = varName2idx.find(key);
    if (it != varName2idx.end())
      return it->second;
    int v = S.newVar();
    varName2idx[key] = v;
    logger->trace("Created new var {} for key '{}'", v, key);
    return v;
  };

  auto constLit = [&](bool value) -> Glucose::Lit {
    const std::string key = value ? "$__CONST_TRUE__" : "$__CONST_FALSE__";
    int v = getOrCreateVar(key);
    Glucose::Lit lv = Glucose::mkLit(v);
    S.addClause(value ? lv : ~lv);
    logger->trace("Added constant clause for {} as var {}", value, v);
    return lv;
  };

  struct Frame {
    BoolExpr* expr;
    bool visited = false;
    Glucose::Lit leftLit, rightLit;
  };

  std::stack<Frame> stk;
  stk.push({root, false, {}, {}});
  std::unordered_map<BoolExpr*, Glucose::Lit> result;

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
      if (name == "0" || name == "false" || name == "False" || name == "FALSE")
        lit = constLit(false);
      else if (name == "1" || name == "true" || name == "True" ||
               name == "TRUE")
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
      if (e->getRight())
        stk.push({e->getRight(), false, {}, {}});
      if (e->getLeft())
        stk.push({e->getLeft(), false, {}, {}});
      continue;
    }

    // Children have been processed; retrieve their lits
    if (e->getLeft())
      fr.leftLit = result[e->getLeft()];
    if (e->getRight())
      fr.rightLit = result[e->getRight()];

    // Create fresh var for this gate
    int v = S.newVar();
    Glucose::Lit lit_v = Glucose::mkLit(v);
    node2var[e] = v;
    result[e] = lit_v;

    logger->trace("Encoding node op={} as var {}", static_cast<int>(e->getOp()),
                  v);

    // Emit Tseitin clauses
    switch (e->getOp()) {
      case Op::NOT:
        S.addClause(~lit_v, ~fr.leftLit);
        S.addClause(lit_v, fr.leftLit);
        break;
      case Op::AND:
        S.addClause(~lit_v, fr.leftLit);
        S.addClause(~lit_v, fr.rightLit);
        S.addClause(lit_v, ~fr.leftLit, ~fr.rightLit);
        break;
      case Op::OR:
        S.addClause(~fr.leftLit, lit_v);
        S.addClause(~fr.rightLit, lit_v);
        S.addClause(~lit_v, fr.leftLit, fr.rightLit);
        break;
      case Op::XOR:
        S.addClause(~lit_v, ~fr.leftLit, ~fr.rightLit);
        S.addClause(~lit_v, fr.leftLit, fr.rightLit);
        S.addClause(lit_v, ~fr.leftLit, fr.rightLit);
        S.addClause(lit_v, fr.leftLit, ~fr.rightLit);
        break;
      default:
        logger->warn("Unhandled operator in tseitinEncode: {}",
                     static_cast<int>(e->getOp()));
        break;
    }

    stk.pop();
  }

  logger->debug("Finished Tseitin encode");
  return result.at(root);
}

}  // namespace

 MiterStrategy::MiterStrategy(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1, const std::string& logFileName, const std::string& prefix)
      : prefix_(prefix) {
    top0_ = top0;
    top1_ = top1;
    logFileName_ = logFileName;
  }

size_t MiterStrategy::normalizeInputs(
    std::vector<naja::DNL::DNLID>& inputs0,
    std::vector<naja::DNL::DNLID>& inputs1,
    const std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>&
        inputs0Map,
    const std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>&
        inputs1Map) {
  ensureLoggerInitialized();
  logger->info("normalizeInputs: starting");

  // find the intersection of inputs0 and inputs1 based on the getFullPathIDs of
  // DNLTerminal and the diffs
  
  std::set<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> paths0;
  std::set<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> paths1;
  std::set<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> pathsCommon;

  for (const auto& [path0, input0] : inputs0Map) {
    paths0.insert(path0);
  }
  for (const auto& [path1, input1] : inputs1Map) {
    paths1.insert(path1);
  }
  size_t index = 0;
  for (const auto& [path0, input0] : inputs0Map) {
    if (paths1.find(path0) != paths1.end()) {
      pathsCommon.insert(path0 );
    }
  }
  std::vector<naja::DNL::DNLID> diff0;
  for (const auto& [path0, input0] : inputs0Map) {
    if (pathsCommon.find(path0) == pathsCommon.end()) {
      diff0.emplace_back(input0);
      const auto&pathInstance = path0;
      std::string pathString = "";
      for (const auto& name : pathInstance.first) {
        pathString += name.getString() + ".";
      }
      logger->info("diff0 input: {}", pathString);
    }
  }
  std::vector<naja::DNL::DNLID> diff1;
  for (const auto& [path1, input1] : inputs1Map) {
    if (pathsCommon.find(path1) == pathsCommon.end()) {
      diff1.emplace_back(input1);
      const auto&pathInstance = path1;
      std::string pathString = "";
      for (const auto& name : pathInstance.first) {
        pathString += name.getString() + ".";
      }
      logger->info("diff1 input: {}", pathString);
    }
  }
  inputs0.clear();
  for (const auto& path : pathsCommon) {
    inputs0.emplace_back(inputs0Map.at(path));
  }
  inputs0.insert(inputs0.end(), diff0.begin(), diff0.end());
  #ifdef DEBUG_CHECKS
  for (size_t i = 0; i < inputs0.size(); ++i) {
    logger->debug("normalized input0[{}]: {}:{}", i, 
      naja::DNL::get()->getDNLTerminalFromID(inputs0[i]).getDNLInstance().getPath().getString(), 
      naja::DNL::get()->getDNLTerminalFromID(inputs0[i]).getSnlBitTerm()->getName().getString());
  }
  #endif
  inputs1.clear();
  for (const auto& path : pathsCommon) {
    inputs1.emplace_back(inputs1Map.at(path));
  }
  inputs1.insert(inputs1.end(), diff1.begin(), diff1.end());
  #ifdef DEBUG_CHECKS
  for (size_t i = 0; i < inputs1.size(); ++i) {
    logger->debug("normalized input1[{}]: DNLID {}:{}", i, 
      naja::DNL::get()->getDNLTerminalFromID(inputs1[i]).getDNLInstance().getPath().getString(), 
      naja::DNL::get()->getDNLTerminalFromID(inputs1[i]).getSnlBitTerm()->getName().getString());
  }
  #endif
  logger->info("size of common inputs: {}", pathsCommon.size());
  logger->info("size of diff0 inputs: {}", diff0.size());
  logger->info("size of diff1 inputs: {}", diff1.size());
  return pathsCommon.size();
}

void MiterStrategy::normalizeOutputs(
    std::vector<naja::DNL::DNLID>& outputs0,
    std::vector<naja::DNL::DNLID>& outputs1,
    const std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>&
        outputs0Map,
    const std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>&
        outputs1Map) {
  ensureLoggerInitialized();
  logger->info("normalizeOutputs: starting");

  // find the intersection of outputs0 and outputs1 based on the getFullPathIDs
  // of DNLTerminal and the diffs
  
  std::set<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> paths0;
  std::set<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> paths1;
  std::set<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> pathsCommon;
  for (const auto& [path1, output1] : outputs1Map) {
    paths1.insert(path1);
  }
  size_t index = 0;
  for (const auto& [path0, output0] : outputs0Map) {
    if (paths1.find(path0) != paths1.end()) {
      pathsCommon.insert(path0 );
    }
  }
  std::vector<naja::DNL::DNLID> diff0;
  for (const auto& [path0, output0] : outputs0Map) {
    if (pathsCommon.find(path0) == pathsCommon.end()) {
      diff0.emplace_back(output0);
      std::string fullName;
      for (const auto& name : path0.first) {
        fullName += name.getString() + ".";
      }
      fullName += std::to_string(path0.second[0]) + ".";
      fullName += std::to_string(path0.second[1]);
      logger->info("Will ignore the analysis for: {} from netlist 0 as it does not exist in netlist 1", fullName);
    }
  }
  std::vector<naja::DNL::DNLID> diff1;
  for (const auto& [path1, output1] : outputs1Map) {
    if (pathsCommon.find(path1) == pathsCommon.end()) {
      std::string fullName;
      for (const auto& name : path1.first) {
        fullName += name.getString() + ".";
      }
      fullName += std::to_string(path1.second[0]) + ".";
      fullName += std::to_string(path1.second[1]);
      logger->info("Will ignore the analysis for: {} from netlist 1 as it does not exist in netlist 0", fullName);
      diff1.emplace_back(output1);
    }
  }
  outputs0.clear();
  for (const auto& path : pathsCommon) {
    outputs0.emplace_back(outputs0Map.at(path));
  }
  //outputs0.insert(outputs0.end(), diff0.begin(), diff0.end());
  outputs1.clear();
  for (const auto& path : pathsCommon) {
    outputs1.emplace_back(outputs1Map.at(path));
  }
  //outputs1.insert(outputs1.end(), diff1.begin(), diff1.end());
  logger->info("size of common outputs: {}", pathsCommon.size());
  logger->info("size of diff0 outputs: {}", diff0.size());
  logger->info("size of diff1 outputs: {}", diff1.size());
  #ifdef DEBUG_CHECKS
  if (outputs0.size() == outputs1.size()) {
    if (outputs0 != outputs1) {
      // build the paths vector for outputs0 and outputs1
      std::vector<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> paths0;
      std::vector<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>> paths1;
      for (const auto& output0 : outputs0) {
        std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>> path;
        for (const auto& [path0, output0m] : outputs0Map) {
          if (output0m == output0) {
            path = path0;
            break;
          }
        }
        paths0.emplace_back(path);
      }
      for (const auto& output1 : outputs1) {
        std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>> path;
        for (const auto& [path1, output1m] : outputs1Map) {
          if (output1m == output1) {
            path = path1;
            break;
          }
        }
        paths1.emplace_back(path);
      }
      if (paths0 != paths1) {
        logger->error("Miter outputs must match in order");
        // for (size_t i = 0; i < paths0.size(); ++i) {
        //   naja::NL::SNLPath p0 = naja::NL::SNLPath(
        //       top0_, paths0[i]);
        //   naja::NL::SNLPath p1 = naja::NL::SNLPath(
        //       top1_, paths1[i]);
        //   throw std::runtime_error("Output " + std::to_string(i) +
        //                            " mismatch: " + p0.getString() + " vs " +
        //                            p1.getString());
        // }
        assert(false && "Miter outputs must match in order");
      }
    }
  }
  #endif
}

void MiterStrategy::init() {
  ensureLoggerInitialized();
  logger->info("MiterStrategy::run starting");
  // build both sets of POs
  topInit_ = NLUniverse::get()->getTopDesign();
  NLUniverse* univ = NLUniverse::get();
  naja::DNL::destroy();
  univ->setTopDesign(top0_);
  logger->info("Collecting POs for design 0: {}\n", top0_->getName().getString().c_str());
  builder0_.collect();
  logger->info("Collected {} POs for design 0\n", builder0_.getOutputs().size());
  PIs0_ = builder0_.getInputs();
  naja::DNL::destroy();
  univ->setTopDesign(top1_);
  logger->info("Collecting POs for design 1: {}\n", top1_->getName().getString().c_str());
  builder1_.collect();
  logger->info("Collected {} POs for design 1\n", builder1_.getOutputs().size());
  PIs1_ = builder1_.getInputs();
}

bool MiterStrategy::run() {
  NLUniverse* univ = NLUniverse::get();
  // normalize inputs and outputs
  std::vector<naja::DNL::DNLID> inputs0sort;
  std::vector<naja::DNL::DNLID> inputs1sort;
  std::vector<naja::DNL::DNLID> outputs0sort;
  std::vector<naja::DNL::DNLID> outputs1sort;
  size_t commonSize = normalizeInputs(inputs0sort, inputs1sort, builder0_.getInputsMap(),
                  builder1_.getInputsMap());
  builder0_.setLastCommonID(commonSize > 0 ? commonSize + 2 : 1);
  builder1_.setLastCommonID(commonSize > 0 ? commonSize + 2 : 1);
  lastCommonVarID_ = commonSize -1 > 0 ? (commonSize - 1) + 2/* 0, 1*/ : 1 /* 0, 1*/;
  normalizeOutputs(outputs0sort, outputs1sort, builder0_.getOutputsMap(),
                   builder1_.getOutputsMap());
  
  logger->info("size of PIs in circuit 0: {}", inputs0sort.size());
  logger->info("size of PIs in circuit 1: {}", inputs1sort.size());
  logger->info("size of POs in circuit 0: {}", outputs0sort.size());
  logger->info("size of POs in circuit 1: {}", outputs1sort.size());
  // return false;
  naja::DNL::destroy();
  univ->setTopDesign(top0_);
  builder0_.setInputs(inputs0sort);
  builder0_.setOutputs(outputs0sort);
  naja::DNL::destroy();
  univ->setTopDesign(top1_);
  builder1_.setInputs(inputs1sort);
  builder1_.setOutputs(outputs1sort);
  naja::DNL::destroy();
  univ->setTopDesign(top0_);
  builder0_.build();
  const auto& PIs0 = builder0_.getInputs();
  const auto& POs0 = builder0_.getPOs();
  const auto& outputs0 = builder0_.getOutputs();
  const auto& inputs2inputsIDs0 = builder0_.getInputs2InputsIDs();
  const auto&outputs2outputsIDs0 = builder0_.getOutputs2OutputsIDs();
  naja::DNL::destroy();
  univ->setTopDesign(top1_);
  builder1_.build();
  const auto& PIs1 = builder1_.getInputs();
  const auto& POs1 = builder1_.getPOs();
  const auto& outputs1 = builder1_.getOutputs();
  const auto& inputs2inputsIDs1 = builder1_.getInputs2InputsIDs();
  const auto& outputs2outputsIDs1 = builder1_.getOutputs2OutputsIDs();

  // print path to var names
  const auto & inputs2DnlIds = builder0_.getInputs();
  // var names for inputs
  const auto & varNames = builder0_.getTermDNLID2VarID();
  for (size_t i = 0; i < inputs2DnlIds.size(); ++i) {
    const auto&path = builder0_.getInputs2InputsIDs().at(builder0_.getDNLIDforInput(i));
    logger->info("VARID {} DNLID {}", varNames[inputs2DnlIds[i]], inputs2DnlIds[i]);
    for (const auto& name : path.first) {
      logger->info("{}.", name.getString().c_str());
    }
    for (const auto& id : path.second) {
      logger->info("bit: {}.", id);
    }
    logger->info("\n");
  }
  // same for builder1
  const auto & inputs2DnlIds1 = builder1_.getInputs();
  const auto & varNames1 = builder1_.getTermDNLID2VarID();
  for (size_t i = 0; i < inputs2DnlIds1.size(); ++i) {
    const auto& path = builder1_.getInputs2InputsIDs().at(builder1_.getDNLIDforInput(i));
    logger->info("VARID {} DNLID {}", varNames1[inputs2DnlIds1[i]], inputs2DnlIds1[i]);
    for (const auto& name : path.first) {
      logger->info("{}.", name.getString().c_str());
    }
    for (const auto& id : path.second) {
      logger->info("bit: {}.", id);
    }
    logger->info("\n");
  }

  if (topInit_ != nullptr) {
    univ->setTopDesign(topInit_);
  }

  if (POs0.empty() || POs1.empty()) {
    logger->warn(
        "No primary outputs found on one of the designs; aborting run");
    return false;
  }

  // build the Boolean-miter expression
  logger->info("Building miter expression");
  auto miter = buildMiter(POs0, POs1);
  logger->info("Finished building miter expression");

  // Now SAT check via Glucose
  Glucose::SimpSolver solver;

  // mappings for Tseitin encoding
  std::unordered_map<BoolExpr*, int> node2var;
  std::unordered_map<std::string, int> varName2idx;

  // Tseitin-encode & get the literal for the root
  Glucose::Lit rootLit = tseitinEncode(solver, miter, node2var, varName2idx);

  // Assert root == true
  solver.addClause(rootLit);

  // solve with no assumptions
  logger->info("Started Glucose solving");
  bool sat = solver.solve();
  logger->info("Finished Glucose solving: {}", sat ? "SAT" : "UNSAT");

  if (sat) {
    logger->warn("Miter found a difference -> moving to analyze individual POs");
    for (size_t i = 0; i < POs0.size(); ++i) {
      if (POs0[i] == POs1[i]) { // We can do this comparison because of the caching in, if they are the same, they are the same pointer
        logger->info("PO index {} expressions are equal; skipping", i);
        continue;
      }
      if (builder0_.getOutputs2OutputsIDs().at(builder0_.getDNLIDforOutput(i)) !=
          builder1_.getOutputs2OutputsIDs().at(builder1_.getDNLIDforOutput(i))) {
        // LCOV_EXCL_START
        const auto&path0 = builder0_.getOutputs2OutputsIDs().at(builder0_.getDNLIDforOutput(i));
        const auto&path1 = builder1_.getOutputs2OutputsIDs().at(builder1_.getDNLIDforOutput(i));
        // print path0
        for (const auto& name : path0.first) {
          logger->info("{}.", name.getString().c_str());
        }
        for (const auto& id : path0.second) {
          logger->info("bit: {}.", id);
        }
        logger->info("\n");
        // print path1
        for (const auto& name : path1.first) {
          logger->info("{}.", name.getString().c_str());
        }
        for (const auto& id : path1.second) {
          logger->info("bit: {}.", id);
        }
        logger->info("\n");
        throw std::runtime_error("Miter PO index " + std::to_string(i) +
                                 " DNLIDs do not match");
        // LCOV_EXCL_STOP
      }
      tbb::concurrent_vector<BoolExpr*> singlePOs0S;
      singlePOs0S.emplace_back(POs0[i]);
      tbb::concurrent_vector<BoolExpr*> singlePOs1S;
      singlePOs1S.emplace_back(POs1[i]);
      auto singleMiter = buildMiter(singlePOs0S, singlePOs1S);

      std::unordered_map<BoolExpr*, int> singleNode2var;
      std::unordered_map<std::string, int> singleVarName2idx;
      // Tseitin-encode the single miter
      Glucose::SimpSolver singleSolver;
      Glucose::Lit singleRootLit = tseitinEncode(
          singleSolver, singleMiter, singleNode2var, singleVarName2idx);

      singleSolver.addClause(singleRootLit);
      if (singleSolver.solve()) {
        bool unSupportedVar = false;
        const auto&varSupportA = POs0[i]->getSupportVars();
        for (const auto&var : varSupportA) {
          if (lastCommonVarID_ < var) {
            logger->warn("Unsupported var for PO0: {}", var);  
            unSupportedVar = true;
          }
        }
        const auto&varSupportB = POs1[i]->getSupportVars();
        for (const auto&var : varSupportB) {
          if (lastCommonVarID_ < var) {
            logger->warn("Unsupported var for PO1: {}", var);  
            unSupportedVar = true;
          }
        }
        if (unSupportedVar) {
          logger->warn("buildMiter skipping output index {} due to unsupported variable", i);
          continue;
        }
        failedPOs_.emplace_back(i);
        logger->info("Found difference for PO: {}", i);
        logger->debug("Clause 0 {}", POs0[i]->toString());
        logger->debug("Clause 1 {}", POs1[i]->toString());
        // print path of index i
        const auto&path0 = builder0_.getOutputs2OutputsIDs().at(builder0_.getDNLIDforOutput(i));
        std::string pathString = "";
        for (const auto& name : path0.first) {
          pathString += name.getString() + ".";
        }
        for (const auto& id : path0.second) {
          pathString += std::to_string(id) + ".";
        }
        const auto&terminal0 = naja::DNL::get()->getDNLTerminalFromID(outputs0[i]);
        logger->info("Path of differing PO {}: {}", i, pathString);
        const auto&path1 = builder1_.getOutputs2OutputsIDs().at(builder1_.getDNLIDforOutput(i));
        std::string pathString1 = "";
        for (const auto& name : path1.first) {
          pathString1 += name.getString() + ".";
        }
        for (const auto& id : path1.second) {
          pathString1 += std::to_string(id) + ".";
        }
        logger->info("Path of differing PO {}: {}", i, pathString1);
        std::vector<naja::NL::SNLDesign*> topModels;
        topModels.emplace_back(top0_);
        topModels.emplace_back(top1_);
        std::vector<std::vector<naja::DNL::DNLID>> PIs;
        PIs.emplace_back(PIs0);
        PIs.emplace_back(PIs1);
        naja::NL::SNLEquipotential::Terms terms0;
        naja::NL::SNLEquipotential::Terms terms1;
        std::unordered_set<std::string> insTerms0;
        std::unordered_set<std::string> insTerms1;
        for (size_t j = 0; j < topModels.size(); ++j) {
          DNL::destroy();
          NLUniverse::get()->setTopDesign(topModels[j]);
          // if (j == 0) {
          //   //logger->info("$$$ 0 term {} of model {}",
          //   naja::DNL::get()->getDNLTerminalFromID(outputs0[i]).getSnlBitTerm()->getName().getString().c_str(),
          //   naja::DNL::get()->getDNLTerminalFromID(outputs0[i]).getSnlBitTerm()->getDesign()->getName().getString().c_str());
          // } else {
          //   //logger->info("### 0 term {} of model {}",
          //   naja::DNL::get()->getDNLTerminalFromID(outputs1[i]).getSnlBitTerm()->getName().getString().c_str(),
          //   naja::DNL::get()->getDNLTerminalFromID(outputs1[i]).getSnlBitTerm()->getDesign()->getName().getString().c_str());

          // }
          // if (dnls_.size() <= j) {
          //   dnls_.emplace_back(*naja::DNL::get());
          // }
          SNLLogicCone cone(j == 0 ? outputs0[i] : outputs1[i], PIs[j],
                            naja::DNL::get());
          cone.run();
          // std::string dotFileNameEquis(
          //     std::string(prefix_ + "_" +
          //     DNL::get()->getDNLTerminalFromID(outputs0[i]).getSnlBitTerm()->getName().getString()
          //     + std::to_string(outputs0[i]) + "_" +std::to_string(j) +
          //     std::string(".dot")));
          // std::string svgFileNameEquis(
          //     std::string(prefix_ + "_" +
          //     DNL::get()->getDNLTerminalFromID(outputs0[i]).getSnlBitTerm()->getName().getString()
          //     + std::to_string(outputs0[i]) + "_" + std::to_string(j) +
          //     std::string(".svg")));
          // SnlVisualiser snl2(topModels[j], cone.getEquipotentials());
          // for (const auto& equi : cone.getEquipotentials()) {
          //   for (const auto& term : equi.getTerms()) {
          //     if (j == 0) {
          //       terms0.insert(term);
          //       // logger->info("$$$ Term 0: {}", term->getString().c_str());
          //     } else {
          //       terms1.insert(term);
          //       // logger->info("### Term 1: {}", term->getString().c_str());
          //     }
          //   }
          //   for (const auto& termOcc : equi.getInstTermOccurrences()) {
          //     if (j == 0) {
          //       insTerms0.insert(termOcc);
          //       // logger->info("$$$ Inst Term 0: {}",
          //       // termOcc.getString().c_str());
          //     } else {
          //       insTerms1.insert(termOcc);
          //       // logger->info("### Inst Term 1: {}",
          //       // termOcc.getString().c_str());
          //     }
          //   }
          // }
          for (const auto& DNLID : cone.getCollectedTerms()) {
            const naja::DNL::DNLTerminalFull& termFull =
                naja::DNL::get()->getDNLTerminalFromID(DNLID);
            if (termFull.isTopPort()) {
              if (j == 0) {
                terms0.insert(termFull.getSnlBitTerm());
              } else {
                terms1.insert(termFull.getSnlBitTerm());
              }
            } else {
              std::string fullPath;
              SNLDesign* design = NLUniverse::get()->getTopDesign();
              fullPath += design->getName().getString() + ".";
              const auto& idPath = termFull.getFullPathIDs();
              for (size_t i = 0 ; i < idPath.size() - 2; ++i) {
                fullPath += design->getInstance(idPath[i])->getName().getString() + ".";
                design = design->getInstance(idPath[i])->getModel();         
              }
              fullPath += design->getTerm(idPath[idPath.size() - 2])->getName().getString() + ".";
              fullPath +=
                  std::to_string(termFull.getFullPathIDs().back());
              if (j == 0) {
                insTerms0.insert(fullPath);
              } else {
                insTerms1.insert(fullPath);
              }
            }
          }
          // snl2.process();
          // snl2.getNetlistGraph().dumpDotFile(dotFileNameEquis.c_str());
          // executeCommand(std::string(std::string("dot -Tsvg ") +
          //                            dotFileNameEquis + std::string(" -o ") +
          //                            svgFileNameEquis).c_str());
          // logger->info("svg file name: {}", svgFileNameEquis);
        }

        // find intersection and diff of terms0 and terms1
        naja::NL::SNLEquipotential::Terms termsCommon;
        naja::NL::SNLEquipotential::Terms termsDiff;
        for (const auto& term0 : terms0) {
          bool found = false;
          for (const auto& term1 : terms1) {
            if (term0->getName().getString() == term1->getName().getString() &&
                term0->getBit() == term1->getBit()) {
              found = true;
              break;
            }
          }
          if (found) {
            termsCommon.insert(term0);
          } else {
            termsDiff.insert(term0);
            if (term0->getDirection() ==
                naja::NL::SNLBitTerm::Direction::Output) {
              continue;
            }
            logger->info("Diff 0 term: {}", term0->getString());
          }
        }
        for (const auto& term1 : terms1) {
          bool found = false;
          for (const auto& term0 : terms0) {
            if (term0->getName().getString() == term1->getName().getString() &&
                term0->getBit() == term1->getBit()) {
              found = true;
              break;
            }
          }
          if (!found) {
            termsDiff.insert(term1);
            if (term1->getDirection() ==
                naja::NL::SNLBitTerm::Direction::Output) {
              continue;
            }
            logger->info("Diff 1 term: {}", term1->getString());
          }
        }
        // print termsDiff
        // for (const auto& term : termsDiff) {
        //   if (term->getDirection() ==
        //   naja::NL::SNLBitTerm::Direction::Output) {
        //     continue;
        //   }
        //   logger->info("Diff term: {}", term->getString());
        // }
        // find intersection and diff of insTerms0 and insTerms1
        std::set<std::string> insTermsCommon;
        std::set<std::string> insTermsDiff;
        for (const auto& term0 : insTerms0) {
          bool found = false;
          if (insTerms1.find(term0) != insTerms1.end()) {
            found = true;
          }
          if (found) {
            insTermsCommon.insert(term0);
          } else {
            insTermsDiff.insert(term0);
            logger->info("Diff 0 inst term {}",
                         term0);
          }
        }
        for (const auto& term1 : insTerms1) {
          bool found = false;
          if (insTerms0.find(term1) != insTerms0.end()) {
            found = true;
          }
          if (!found) {
            insTermsDiff.insert(term1);
            logger->info("Diff 1 inst term {}",
                         term1);
          }
        }

        logger->info("size of intersection of terms: {}", termsCommon.size());
        logger->info("size of diff of terms: {}", termsDiff.size());
        logger->info("size of intersection of inst terms: {}",
                      insTermsCommon.size());
        logger->info("size of diff of inst terms: {}", insTermsDiff.size());
      }
    }
  }
  if (topInit_ != nullptr) {
    univ->setTopDesign(topInit_);
  }
  // if UNSAT → miter can never be true → outputs identical
  logger->info("Circuits are {}", sat ? "DIFFERENT" : "IDENTICAL");
  return !sat;
}

BoolExpr* MiterStrategy::buildMiter(
    const tbb::concurrent_vector<BoolExpr*>& A,
    const tbb::concurrent_vector<BoolExpr*>& B) const {
  ensureLoggerInitialized();
  logger->debug("buildMiter: A.size={} B.size={}", A.size(), B.size());

  // Empty miter = always-false (no outputs to compare)
  if (A.empty()) {
    logger->error("buildMiter called with empty A");
    assert(false);
    return BoolExpr::createFalse();
  }

  // Start with the first XOR
  auto miter = BoolExpr::createFalse();

  // OR in the rest
  for (size_t i = 0; i < A.size(); ++i) {
    if (B.size() <= i) {
      logger->warn("Miter different number of outputs: {} vs {}", A.size(),
                   B.size());
      break;
    }
    // bool unSupportedVar = false;
    // const auto&varSupportA = A[i]->getSupportVars();
    // for (const auto&var : varSupportA) {
    //   if (lastCommonVarID_ < var) {
    //     logger->warn("Unsupported var: {}", var);  
    //     unSupportedVar = true;
    //   }
    // }
    // const auto&varSupportB = B[i]->getSupportVars();
    // for (const auto&var : varSupportB) {
    //   if (lastCommonVarID_ < var) {
    //     logger->warn("Unsupported var: {}", var);  
    //     unSupportedVar = true;
    //   }
    // }
    // if (unSupportedVar) {
    //   logger->warn("buildMiter skipping output index {} due to unsupported variable", i);
    //   continue;
    // }
    auto diff = BoolExpr::Xor(A[i], B[i]);
    miter = BoolExpr::Or(miter, diff);
  }
  // logger->trace("buildMiter produced expression: {}", miter->toString());
  return miter;
}
