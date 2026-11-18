// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "BoolExprCache.h"
#include "DNL.h"
#include "DesignBoundary.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstTerm.h"
#include "SNLInstance.h"
#include "SNLSVConstructor.h"
#include "SNLScalarNet.h"
#include "Tree2BoolExpr.h"
#include "model/SequentialDesignModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC {
namespace {

using namespace naja::NL;

// A 16 x 8 RAM behind a registered read port. `valid` gives every design a
// reset register, so the reset sequence has something to initialize.
std::string ramModule(const std::string& name, const std::string& body) {
  return "module " + name + R"((
  input  wire       clk,
  input  wire       rst,
  input  wire       we,
  input  wire [3:0] waddr,
  input  wire [7:0] wdata,
  input  wire [3:0] raddr,
  output reg  [7:0] rdata,
  output reg        valid
);
)" + body + R"(
  always @(posedge clk) begin
    if (rst) valid <= 1'b0;
    else valid <= we;
  end
endmodule
)";
}

const char* const kBaseBody = R"(
  reg [7:0] mem [0:15];
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata;
    rdata <= mem[raddr];
  end)";

// The same RAM written differently: reversed array range and renamed wires.
const char* const kRewrittenBody = R"(
  reg [7:0] mem [15:0];
  wire [3:0] wa = waddr[3:0];
  wire [3:0] ra = raddr + 4'd0;
  always @(posedge clk) begin
    if (we == 1'b1) mem[wa] <= wdata;
    rdata <= mem[ra];
  end)";

const char* const kWriteDataBugBody = R"(
  reg [7:0] mem [0:15];
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata ^ 8'h01;
    rdata <= mem[raddr];
  end)";

// The write address differs from the base only while nothing is written.
const char* const kIdleWriteAddressBody = R"(
  reg [7:0] mem [0:15];
  wire [3:0] wa = we ? waddr : 4'd0;
  always @(posedge clk) begin
    if (we) mem[wa] <= wdata;
    rdata <= mem[raddr];
  end)";

const char* const kRenamedMemoryBody = R"(
  reg [7:0] storage [0:15];
  always @(posedge clk) begin
    if (we) storage[waddr] <= wdata;
    rdata <= storage[raddr];
  end)";

// Word 0 is written with a constant while reset is asserted.
std::string resetWriteBody(const char* constant) {
  return std::string(R"(
  reg [7:0] mem [0:15];
  always @(posedge clk) begin
    if (rst) mem[4'd0] <= )") + constant + R"(;
    else if (we) mem[waddr] <= wdata;
    rdata <= mem[raddr];
  end)";
}

// Equivalent to the base after the reset sequence, but each write reaches the
// memory one cycle later and a bypass covers the gap.
const char* const kLaggedWriteBody = R"(
  reg [7:0] mem [0:15];
  reg       we_q;
  reg [3:0] waddr_q;
  reg [7:0] wdata_q;
  always @(posedge clk) begin
    if (rst) we_q <= 1'b0;
    else we_q <= we;
    waddr_q <= waddr;
    wdata_q <= wdata;
    if (we_q) mem[waddr_q] <= wdata_q;
    rdata <= (we_q && waddr_q == raddr) ? wdata_q : mem[raddr];
  end)";

// The read address is the read data itself.
const char* const kReadLoopBody = R"(
  reg [7:0] mem [0:15];
  wire [7:0] looped = mem[looped[3:0]];
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata;
    rdata <= looped;
  end)";

// One read addresses another: an order exists, so this is no loop.
const char* const kChainedReadBody = R"(
  reg [7:0] mem [0:15];
  wire [7:0] first = mem[raddr];
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata;
    rdata <= mem[first[3:0]];
  end)";

// A second write process whose enable is the given literal.
std::string secondWriteBody(const char* enable) {
  return std::string(R"(
  reg [7:0] mem [0:15];
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata;
    rdata <= mem[raddr];
  end
  always @(posedge clk) begin
    if ()") + enable + R"() mem[4'd0] <= 8'h55;
  end)";
}

// The write data depends on a register that nothing ever loads.
const char* const kUndefinedWriteDataBody = R"(
  reg [7:0] mem [0:15];
  reg [7:0] key;
  always @(posedge clk) begin
    key <= key;
    if (we) mem[waddr] <= wdata ^ key;
    rdata <= mem[raddr];
  end)";

// Two RAMs in submodules, the second one with a write-data bug on demand.
std::string twoRamModules(const std::string& name, const char* flip) {
  return "module " + name + R"(_ram #(parameter FLIP = 8'h00) (
  input  wire       clk,
  input  wire       we,
  input  wire [3:0] waddr,
  input  wire [7:0] wdata,
  input  wire [3:0] raddr,
  output reg  [7:0] rdata
);
  reg [7:0] mem [0:15];
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata ^ FLIP;
    rdata <= mem[raddr];
  end
endmodule

module )" + name + R"((
  input  wire       clk,
  input  wire       rst,
  input  wire       we,
  input  wire       sel,
  input  wire [3:0] waddr,
  input  wire [7:0] wdata,
  input  wire [3:0] raddr,
  output wire [7:0] rdata,
  output reg        valid
);
  wire [7:0] ra, rb;
  reg sel_q;
  )" + name + R"(_ram u_a (.clk(clk), .we(we & ~sel), .waddr(waddr),
      .wdata(wdata), .raddr(raddr), .rdata(ra));
  )" + name + R"(_ram #(.FLIP()" + flip + R"()) u_b (.clk(clk), .we(we & sel),
      .waddr(waddr), .wdata(wdata), .raddr(raddr), .rdata(rb));
  assign rdata = sel_q ? rb : ra;
  always @(posedge clk) begin
    sel_q <= sel;
    if (rst) valid <= 1'b0;
    else valid <= we;
  end
endmodule
)";
}

class MemoryAbstractionTests : public ::testing::Test {
 protected:
  void SetUp() override {
    auto* universe = NLUniverse::create();
    auto* db = NLDB::create(universe);
    designs_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    // The name of a parameterized test contains a slash.
    std::string name =
        ::testing::UnitTest::GetInstance()->current_test_info()->name();
    std::replace(name.begin(), name.end(), '/', '_');
    directory_ =
        std::filesystem::temp_directory_path() / ("kf_memory_abstraction_" + name);
    std::filesystem::remove_all(directory_);
    std::filesystem::create_directories(directory_);
  }

  void TearDown() override {
    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    Tree2BoolExpr::iso2boolExpr_.clear();
    BoolExprCache::destroy();
    std::filesystem::remove_all(directory_);
  }

  SNLDesign* load(const std::string& name, const std::string& body) {
    const auto path = directory_ / (name + ".sv");
    std::ofstream(path) << ramModule(name, body);
    SNLSVConstructor constructor(designs_);
    constructor.construct(path);
    auto* top = designs_->getSNLDesign(NLName(name));
    EXPECT_NE(top, nullptr);
    return top;
  }

  SNLDesign* loadSource(const std::string& name, const std::string& source) {
    const auto path = directory_ / (name + ".sv");
    std::ofstream(path) << source;
    SNLSVConstructor constructor(designs_);
    constructor.construct(path);
    auto* top = designs_->getSNLDesign(NLName(name));
    EXPECT_NE(top, nullptr);
    return top;
  }

  static SecResetSpec resetSpec() {
    SecResetSpec reset;
    reset.cycles = 2;
    reset.ports.push_back({"rst", true});
    return reset;
  }

  SequentialEquivalenceStrategy strategy(SNLDesign* left, SNLDesign* right,
                                         SecEncoding encoding,
                                         SecResetSpec reset = resetSpec()) {
    SequentialEquivalenceStrategy result(
        left, right, Config::SolverType::KISSAT, SecEngine::Pdr, encoding,
        std::move(reset));
    result.setMemoryAbstraction(true);
    return result;
  }

  NLLibrary* designs_ = nullptr;
  std::filesystem::path directory_;
};

class MemoryAbstractionEncodingTests
    : public MemoryAbstractionTests,
      public ::testing::WithParamInterface<SecEncoding> {};

std::vector<std::string> displayNames(const SequentialDesignModel& model,
                                      const std::vector<SignalKey>& keys) {
  std::vector<std::string> names;
  for (const auto& key : keys) {
    names.push_back(model.displayNameByKey.at(key));
  }
  return names;
}

TEST_F(MemoryAbstractionTests, ExtractionReplacesTheMemoryCellsByItsPins) {
  auto* top = load("base", kBaseBody);

  const auto bitLevel = SequentialDesignModel::extract(top);
  ASSERT_FALSE(bitLevel.hasUnsupportedFeatures());
  ASSERT_EQ(bitLevel.bitLevelMemories.size(), 1u);
  EXPECT_EQ(bitLevel.bitLevelMemories.front().displayName, "mem_mem");
  EXPECT_EQ(bitLevel.bitLevelMemories.front().depth, 16u);
  EXPECT_EQ(bitLevel.bitLevelMemories.front().width, 8u);
  EXPECT_TRUE(bitLevel.abstractedMemories.empty());
  // 16 x 8 cells, the read register and `valid`.
  EXPECT_EQ(bitLevel.stateBits.size(), 16u * 8u + 8u + 1u);

  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto abstracted = SequentialDesignModel::extract(top, {}, 0, options);
  ASSERT_FALSE(abstracted.hasUnsupportedFeatures());
  EXPECT_TRUE(abstracted.bitLevelMemories.empty());
  ASSERT_EQ(abstracted.abstractedMemories.size(), 1u);
  const auto& memory = abstracted.abstractedMemories.front();
  EXPECT_EQ(memory.displayName, "mem_mem");
  EXPECT_EQ(memory.depth, 16u);
  EXPECT_EQ(memory.width, 8u);
  EXPECT_NE(memory.signature.find("naja_mem__w8_d16"), std::string::npos)
      << memory.signature;
  EXPECT_EQ(abstracted.stateBits.size(), 8u + 1u);

  // Read data is an input of the model; every other pin is an observed output.
  EXPECT_EQ(memory.readDataKeys.size(), 8u);
  for (const auto& key : memory.readDataKeys) {
    EXPECT_NE(std::find(abstracted.environmentInputs.begin(),
                        abstracted.environmentInputs.end(), key),
              abstracted.environmentInputs.end());
  }
  const auto pins = displayNames(abstracted, memory.interfaceKeys);
  for (const auto* pin : {"mem_mem.WE[0]", "mem_mem.WADDR[3]",
                          "mem_mem.WDATA[7]", "mem_mem.RADDR[0]"}) {
    EXPECT_NE(std::find(pins.begin(), pins.end(), pin), pins.end()) << pin;
  }
  for (const auto& key : memory.interfaceKeys) {
    EXPECT_NE(std::find(abstracted.allObservedOutputs.begin(),
                        abstracted.allObservedOutputs.end(), key),
              abstracted.allObservedOutputs.end());
  }

  // Write address and data only matter while the port writes; the enable and
  // the read address always do.
  const auto guardsOf = [&](const std::string& pin) {
    for (const auto& key : memory.interfaceKeys) {
      if (abstracted.displayNameByKey.at(key) == pin) {
        const auto found = abstracted.observedOutputGuardKeysByKey.find(key);
        return found == abstracted.observedOutputGuardKeysByKey.end()
                   ? std::vector<std::string>{}
                   : displayNames(abstracted, found->second);
      }
    }
    ADD_FAILURE() << "no pin " << pin;
    return std::vector<std::string>{};
  };
  EXPECT_EQ(guardsOf("mem_mem.WADDR[0]"),
            std::vector<std::string>{"mem_mem.WE[0]"});
  EXPECT_FALSE(guardsOf("mem_mem.WDATA[0]").empty());
  EXPECT_EQ(guardsOf("mem_mem.WDATA[0]").front(), "mem_mem.WE[0]");
  EXPECT_TRUE(guardsOf("mem_mem.WE[0]").empty());
  EXPECT_TRUE(guardsOf("mem_mem.RADDR[0]").empty());

  // Only the listed memories are abstracted when a list is given.
  options.abstractOnly.emplace();
  const auto none = SequentialDesignModel::extract(top, {}, 0, options);
  EXPECT_TRUE(none.abstractedMemories.empty());
  EXPECT_EQ(none.bitLevelMemories.size(), 1u);
  options.abstractOnly->insert(memory.instanceKey);
  const auto listed = SequentialDesignModel::extract(top, {}, 0, options);
  EXPECT_EQ(listed.abstractedMemories.size(), 1u);
}

TEST_F(MemoryAbstractionTests, UserSelectedMemoryStaysAUserBoundary) {
  auto* top = load("base", kBaseBody);
  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto model =
      SequentialDesignModel::extract(top, {{"mem_mem", "mem_mem"}}, 0, options);
  ASSERT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_TRUE(model.abstractedMemories.empty());
  EXPECT_TRUE(model.bitLevelMemories.empty());
  // The pins of a user boundary carry its synthetic names.
  const bool hasBoundaryName = std::any_of(
      model.displayNameByKey.begin(), model.displayNameByKey.end(),
      [](const auto& entry) {
        return entry.second.rfind("__kepler_boundary_p0_", 0) == 0;
      });
  EXPECT_TRUE(hasBoundaryName);
}

TEST_P(MemoryAbstractionEncodingTests, ProvesRewrittenRamThroughItsInterface) {
  auto* left = load("base", kBaseBody);
  auto* right = load("rewritten", kRewrittenBody);
  const auto result = strategy(left, right, GetParam()).run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
      << result.reason;
  EXPECT_EQ(result.abstractedMemories, std::vector<std::string>{"mem_mem"});
  EXPECT_TRUE(result.memoryAbstractionFailure.empty());
  // The memory pins are proved but are no outputs of the designs.
  EXPECT_EQ(result.totalOutputs, 9u);
  EXPECT_EQ(result.coveredOutputs, 9u);
}

// Issue 308: a RAM against itself reads a word that was never written. Each
// design's cells used to be independent free values in the binary encoding.
TEST_P(MemoryAbstractionEncodingTests, RamAgainstItselfIsEquivalent) {
  auto* left = load("base", kBaseBody);
  auto* right = load("copy", kBaseBody);
  const auto result = strategy(left, right, GetParam()).run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
      << result.reason;
  EXPECT_EQ(result.abstractedMemories.size(), 1u);
}

TEST_P(MemoryAbstractionEncodingTests, WriteDataBugIsFoundOnTheBitLevelModel) {
  auto* left = load("base", kBaseBody);
  auto* right = load("bug", kWriteDataBugBody);
  const auto result = strategy(left, right, GetParam()).run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different)
      << result.reason;
  // The interface comparison failed and said so; the verdict is bit level.
  EXPECT_FALSE(result.memoryAbstractionFailure.empty());
  EXPECT_TRUE(result.abstractedMemories.empty());
  EXPECT_NE(result.reason.find("rdata["), std::string::npos) << result.reason;
}

TEST_P(MemoryAbstractionEncodingTests, IdleWriteAddressDoesNotBlockTheAbstraction) {
  auto* left = load("base", kBaseBody);
  auto* right = load("idle", kIdleWriteAddressBody);
  const auto result = strategy(left, right, GetParam()).run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
      << result.reason;
  EXPECT_EQ(result.abstractedMemories.size(), 1u)
      << result.memoryAbstractionFailure;
}

TEST_P(MemoryAbstractionEncodingTests, WritesDuringResetAreCompared) {
  auto* zero = load("zero", resetWriteBody("8'h00"));
  auto* zeroCopy = load("zero_copy", resetWriteBody("8'h00"));
  auto* ones = load("ones", resetWriteBody("8'hFF"));

  const auto same = strategy(zero, zeroCopy, GetParam()).run(8);
  EXPECT_EQ(same.status, SequentialEquivalenceStatus::Equivalent) << same.reason;
  EXPECT_EQ(same.abstractedMemories.size(), 1u);

  // After the reset sequence both designs drive the memory alike, and no
  // engine looks at the sequence itself.
  const auto different = strategy(zero, ones, GetParam()).run(8);
  EXPECT_EQ(different.status, SequentialEquivalenceStatus::Different)
      << different.reason;
  EXPECT_NE(different.memoryAbstractionFailure.find("reset sequence"),
            std::string::npos)
      << different.memoryAbstractionFailure;
}

TEST_P(MemoryAbstractionEncodingTests, ExtractedModelsReportAFailedComparison) {
  auto* zero = load("zero", resetWriteBody("8'h00"));
  auto* ones = load("ones", resetWriteBody("8'hFF"));
  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto model0 = SequentialDesignModel::extract(zero, {}, 0, options);
  const auto model1 = SequentialDesignModel::extract(ones, {}, 1, options);
  const auto result =
      strategy(zero, ones, GetParam()).runExtractedModels(model0, model1, 8);
  // Not a verdict on the designs: the caller repeats on bit-level models.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_NE(result.memoryAbstractionFailure.find("reset sequence"),
            std::string::npos)
      << result.memoryAbstractionFailure;
}

TEST_F(MemoryAbstractionTests, MemoriesWithoutATwinKeepTheirCells) {
  auto* left = load("base", kBaseBody);
  auto* right = load("renamed", kRenamedMemoryBody);
  const auto result =
      strategy(left, right, SecEncoding::DualRailSteady).run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
      << result.reason;
  EXPECT_TRUE(result.abstractedMemories.empty());

  // Models that abstract different memories are not a proof problem.
  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto model0 = SequentialDesignModel::extract(left, {}, 0, options);
  const auto model1 = SequentialDesignModel::extract(right, {}, 1, options);
  const auto unpaired = strategy(left, right, SecEncoding::DualRailSteady)
                            .runExtractedModels(model0, model1, 8);
  EXPECT_EQ(unpaired.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_NE(unpaired.memoryAbstractionFailure.find("same memories"),
            std::string::npos)
      << unpaired.memoryAbstractionFailure;
}

TEST_F(MemoryAbstractionTests, OptionOffModelsEveryMemoryBitByBit) {
  auto* left = load("base", kBaseBody);
  auto* right = load("rewritten", kRewrittenBody);
  auto bitLevel = strategy(left, right, SecEncoding::DualRailSteady);
  bitLevel.setMemoryAbstraction(false);
  const auto result = bitLevel.run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
      << result.reason;
  EXPECT_TRUE(result.abstractedMemories.empty());
  EXPECT_TRUE(result.memoryAbstractionFailure.empty());
}

TEST_F(MemoryAbstractionTests,
       BinaryCounterexampleOnUnwrittenMemoryIsNotADifference) {
  auto* left = load("base", kBaseBody);
  auto* right = load("copy", kBaseBody);
  auto* bug = load("bug", kWriteDataBugBody);

  // Bit level, each design's unwritten cells are free values: the engine
  // finds a trace, but no output differs once those cells are unknown.
  auto self = strategy(left, right, SecEncoding::Binary);
  self.setMemoryAbstraction(false);
  const auto result = self.run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive)
      << result.reason;
  EXPECT_TRUE(result.counterexampleDependsOnUninitializedState);
  EXPECT_NE(result.reason.find("uninitialized"), std::string::npos)
      << result.reason;

  // A written word that differs is a difference in either reading, and the
  // trace reported leads to a defined difference.
  auto real = strategy(left, bug, SecEncoding::Binary);
  real.setMemoryAbstraction(false);
  const auto difference = real.run(8);
  EXPECT_EQ(difference.status, SequentialEquivalenceStatus::Different)
      << difference.reason;
  EXPECT_FALSE(difference.counterexampleDependsOnUninitializedState);
  EXPECT_NE(difference.reason.find("Observed output mismatches"),
            std::string::npos)
      << difference.reason;
  EXPECT_NE(difference.reason.find("rdata["), std::string::npos)
      << difference.reason;
}

// Without a reset sequence the binary encoding assumes its property on the
// very first cycle when some state has no initial value. The memory pins are
// checked there instead: `valid` has no initial value, and an enable that
// depends on it is undefined on that cycle.
TEST_F(MemoryAbstractionTests, BinaryChecksThePinsOnTheFirstCycleWithoutReset) {
  const char* const gatedBody = R"(
  reg [7:0] mem [0:15];
  always @(posedge clk) begin
    if (we & valid) mem[waddr] <= wdata;
    rdata <= mem[raddr];
  end)";
  auto* left = load("gated", gatedBody);
  auto* right = load("gated_copy", gatedBody);
  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto model0 = SequentialDesignModel::extract(left, {}, 0, options);
  const auto model1 = SequentialDesignModel::extract(right, {}, 1, options);
  const auto binary = strategy(left, right, SecEncoding::Binary, SecResetSpec{})
                          .runExtractedModels(model0, model1, 4);
  EXPECT_EQ(binary.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_FALSE(binary.memoryAbstractionFailure.empty());

  // The dual-rail encoding assumes nothing: an X enable is no difference.
  const auto dualRail =
      strategy(left, right, SecEncoding::DualRailSteady, SecResetSpec{})
          .runExtractedModels(model0, model1, 4);
  EXPECT_EQ(dualRail.status, SequentialEquivalenceStatus::Equivalent)
      << dualRail.reason << dualRail.memoryAbstractionFailure;
  EXPECT_EQ(dualRail.abstractedMemories.size(), 1u);
}

TEST_F(MemoryAbstractionTests, ReadLoopIsRejectedLikeOnTheBitLevelModel) {
  auto* looped = load("looped", kReadLoopBody);
  auto* loopedCopy = load("looped_copy", kReadLoopBody);
  const char* const reason = "Combinational loop through a memory read port";

  const auto bitLevel = SequentialDesignModel::extract(looped);
  ASSERT_EQ(bitLevel.unsupportedReasons.size(), 1u);
  EXPECT_EQ(bitLevel.unsupportedReasons.front(), reason);
  // Read data is a free input once the memory is abstracted: without this
  // check nothing would relate it to its own address.
  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto abstracted = SequentialDesignModel::extract(looped, {}, 0, options);
  ASSERT_EQ(abstracted.abstractedMemories.size(), 1u);
  ASSERT_EQ(abstracted.unsupportedReasons.size(), 1u);
  EXPECT_EQ(abstracted.unsupportedReasons.front(), reason);

  const auto result =
      strategy(looped, loopedCopy, SecEncoding::DualRailSteady).run(4);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.reason, reason);

  // A read addressed by another read has an order.
  auto* chained = load("chained", kChainedReadBody);
  auto* chainedCopy = load("chained_copy", kChainedReadBody);
  const auto chain = SequentialDesignModel::extract(chained, {}, 0, options);
  EXPECT_FALSE(chain.hasUnsupportedFeatures());
  EXPECT_EQ(chain.abstractedMemories.size(), 1u);
  const auto proved =
      strategy(chained, chainedCopy, SecEncoding::DualRailSteady).run(8);
  EXPECT_EQ(proved.status, SequentialEquivalenceStatus::Equivalent)
      << proved.reason;
  EXPECT_EQ(proved.abstractedMemories.size(), 1u);
}

// The bit-level model reads an X write enable as a port that never writes.
// Compared as a pin, the X would agree with the constant-one enable.
TEST_F(MemoryAbstractionTests, XWriteEnableKeepsTheMemoryCells) {
  auto* undefined = load("undefined_enable", secondWriteBody("1'bx"));
  auto* always = load("always_enable", secondWriteBody("1'b1"));

  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto model = SequentialDesignModel::extract(undefined, {}, 0, options);
  ASSERT_FALSE(model.hasUnsupportedFeatures());
  ASSERT_EQ(model.abstractedMemories.size(), 1u);
  EXPECT_NE(model.abstractedMemories.front().obstacle.find("X literal"),
            std::string::npos)
      << model.abstractedMemories.front().obstacle;

  for (const bool abstraction : {true, false}) {
    auto compared = strategy(undefined, always, SecEncoding::DualRailSteady);
    compared.setMemoryAbstraction(abstraction);
    const auto result = compared.run(8);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different)
        << abstraction << " " << result.reason;
    EXPECT_TRUE(result.abstractedMemories.empty());
  }
}

// A falling-edge register reads after the rising-edge registers have moved:
// through a registered address it sees another word than the rising edge
// does, which one shared read-data value per cycle cannot express.
TEST_P(MemoryAbstractionEncodingTests, FallingEdgeReadKeepsTheMemoryCells) {
  const auto twoReads = [](const char* secondEdge) {
    return std::string(R"(
  reg [7:0] mem [0:15];
  reg [3:0] ra_q;
  reg [7:0] first, second;
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata;
    ra_q <= raddr;
    first <= mem[ra_q];
  end
  always @()") + secondEdge + R"( clk) second <= mem[ra_q];
  always @(posedge clk) rdata <= first ^ second;)";
  };
  auto* falling = load("falling", twoReads("negedge"));
  auto* rising = load("rising", twoReads("posedge"));
  auto* risingCopy = load("rising_copy", twoReads("posedge"));

  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto model = SequentialDesignModel::extract(falling, {}, 0, options);
  ASSERT_EQ(model.abstractedMemories.size(), 1u);
  EXPECT_NE(model.abstractedMemories.front().obstacle.find("falling clock edge"),
            std::string::npos)
      << model.abstractedMemories.front().obstacle;
  const auto risingModel =
      SequentialDesignModel::extract(rising, {}, 0, options);
  ASSERT_EQ(risingModel.abstractedMemories.size(), 1u);
  EXPECT_TRUE(risingModel.abstractedMemories.front().obstacle.empty());

  // The two designs differ exactly there, and the abstraction must not hide it.
  std::vector<std::string> progress;
  auto compared = strategy(falling, rising, GetParam());
  compared.setProgressLogger([&](SecProgressLevel, const std::string& message) {
    progress.push_back(message);
  });
  const auto different = compared.run(8);
  EXPECT_EQ(different.status, SequentialEquivalenceStatus::Different)
      << different.reason;
  EXPECT_TRUE(different.abstractedMemories.empty());
  EXPECT_TRUE(std::any_of(
      progress.begin(), progress.end(), [](const std::string& message) {
        return message.find("memory `mem_mem` of design 1 keeps its cells: its "
                            "read data is sampled on the falling clock edge") !=
               std::string::npos;
      }));

  const auto same = strategy(rising, risingCopy, GetParam()).run(8);
  EXPECT_EQ(same.status, SequentialEquivalenceStatus::Equivalent) << same.reason;
  EXPECT_EQ(same.abstractedMemories.size(), 1u);
}

// A register that reads a memory addressed from another clock domain is
// outside the SEC timing model: the bit-level model skips it, and comparing
// the memory at its pins must not bring it back.
TEST_F(MemoryAbstractionTests, ReadAddressedFromAnotherClockDomainIsSkipped) {
  const auto source = [](const std::string& name) {
    return "module " + name + R"((
  input  wire       clk,
  input  wire       clk2,
  input  wire       rst,
  input  wire       we,
  input  wire [3:0] waddr,
  input  wire [7:0] wdata,
  input  wire [3:0] raddr,
  output reg  [7:0] rdata,
  output reg        valid
);
  reg [7:0] mem [0:15];
  reg [3:0] ra_q;
  always @(posedge clk) begin
    if (we) mem[waddr] <= wdata;
    ra_q <= raddr;
    if (rst) valid <= 1'b0; else valid <= we;
  end
  always @(posedge clk2) rdata <= mem[ra_q];
endmodule
)";
  };
  auto* left = loadSource("crossing", source("crossing"));
  auto* right = loadSource("crossing_copy", source("crossing_copy"));

  for (const bool abstraction : {false, true}) {
    auto compared = strategy(left, right, SecEncoding::DualRailSteady);
    compared.setMemoryAbstraction(abstraction);
    const auto result = compared.run(4);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved)
        << abstraction << " " << result.reason;
    // Only `valid` stays inside one clock domain.
    EXPECT_EQ(result.totalOutputs, 9u);
    EXPECT_EQ(result.coveredOutputs, 1u);
    EXPECT_EQ(result.multiClockDomainSkippedOutputs.size(), 8u);
    EXPECT_EQ(result.abstractedMemories.size(), abstraction ? 1u : 0u);
  }
}

TEST_P(MemoryAbstractionEncodingTests, MemoriesInSubmodulesPairByPath) {
  auto* left = loadSource("pair_a", twoRamModules("pair_a", "8'h00"));
  auto* right = loadSource("pair_b", twoRamModules("pair_b", "8'h00"));
  auto* bug = loadSource("pair_c", twoRamModules("pair_c", "8'h10"));

  auto result = strategy(left, right, GetParam()).run(8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
      << result.reason;
  std::sort(result.abstractedMemories.begin(), result.abstractedMemories.end());
  EXPECT_EQ(result.abstractedMemories,
            (std::vector<std::string>{"u_a/mem_mem", "u_b/mem_mem"}));

  const auto different = strategy(left, bug, GetParam()).run(8);
  EXPECT_EQ(different.status, SequentialEquivalenceStatus::Different)
      << different.reason;
}

// Write data that depends on a register nothing loads is X on both sides:
// no binary-defined difference. The binary encoding has no X and would have
// to assume the two registers equal, so there the memory keeps its cells.
TEST_F(MemoryAbstractionTests, BinaryNeedsDefinedPinsAfterReset) {
  auto* left = load("undefined_data", kUndefinedWriteDataBody);
  auto* right = load("undefined_data_copy", kUndefinedWriteDataBody);

  const auto dualRail =
      strategy(left, right, SecEncoding::DualRailSteady).run(8);
  EXPECT_EQ(dualRail.status, SequentialEquivalenceStatus::Equivalent)
      << dualRail.reason;
  EXPECT_EQ(dualRail.abstractedMemories.size(), 1u);

  const auto binary = strategy(left, right, SecEncoding::Binary).run(8);
  EXPECT_NE(binary.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_NE(binary.status, SequentialEquivalenceStatus::Different)
      << binary.reason;
  EXPECT_NE(binary.memoryAbstractionFailure.find(
                "first cycle after the reset sequence"),
            std::string::npos)
      << binary.memoryAbstractionFailure;
}

// The designs differ only in what the read register samples during the last
// reset cycle, so the outputs differ on the first cycle after reset and never
// again. The binary encoding assumes the outputs agree on that cycle; with
// the memory compared at its pins nothing else would differ, and the
// difference would be lost.
TEST_P(MemoryAbstractionEncodingTests, DifferenceOnTheFirstCycleAfterResetIsKept) {
  const auto bypassed = [](const char* enable) {
    return std::string(R"(
  reg [7:0] mem [0:15];
  wire en = )") + enable + R"(;
  always @(posedge clk) begin
    if (rst) mem[4'd0] <= wdata;
    else if (en) mem[waddr] <= wdata;
    rdata <= (en && waddr == raddr) ? wdata : mem[raddr];
  end)";
  };
  auto* gated = load("gated", bypassed("we & ~rst"));
  auto* plain = load("plain", bypassed("we"));

  for (const bool abstraction : {true, false}) {
    auto compared = strategy(gated, plain, GetParam());
    compared.setMemoryAbstraction(abstraction);
    const auto result = compared.run(8);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different)
        << abstraction << " " << result.reason;
    EXPECT_EQ(result.bound, 0u) << abstraction;
    EXPECT_TRUE(result.abstractedMemories.empty());
    if (abstraction && GetParam() == SecEncoding::Binary) {
      EXPECT_NE(result.memoryAbstractionFailure.find("assumes away"),
                std::string::npos)
          << result.memoryAbstractionFailure;
    }
  }
}

// The lagged design drops a write issued while reset is asserted, which the
// base design performs: the two differ only if reset is asserted again after
// the reset sequence. SEC holds the reset ports released from then on, in
// every PDR query, including for a reset input that joins a cached solver
// after the solver was created. The designs are therefore equivalent.
TEST_F(MemoryAbstractionTests, PdrHoldsResetReleasedAfterTheResetSequence) {
  auto* left = load("base", kBaseBody);
  auto* right = load("lagged", kLaggedWriteBody);
  SecResetSpec reset = resetSpec();
  reset.cycles = 4;
  for (const bool abstraction : {true, false}) {
    auto compared = strategy(left, right, SecEncoding::DualRailSteady, reset);
    compared.setMemoryAbstraction(abstraction);
    const auto result = compared.run(8);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
        << abstraction << " " << result.reason;
  }
}

TEST_F(MemoryAbstractionTests, UnknownResetPortIsFoundWithoutAModel) {
  auto* top = load("base", kBaseBody);
  SecResetSpec reset;
  EXPECT_FALSE(
      SequentialEquivalenceStrategy::findUnknownResetPort(top, reset).has_value());
  reset.cycles = 2;
  reset.ports.push_back({"rst", true});
  EXPECT_FALSE(
      SequentialEquivalenceStrategy::findUnknownResetPort(top, reset).has_value());
  // One bit of a bus input, and a bus named without a bit: the latter is
  // rejected later with its own message.
  reset.ports = {{"waddr[3]", true}, {"wdata", false}};
  EXPECT_FALSE(
      SequentialEquivalenceStrategy::findUnknownResetPort(top, reset).has_value());
  // An output, an out-of-range bit and a name that does not exist.
  for (const auto* name : {"valid", "waddr[4]", "arst_i"}) {
    reset.ports = {{"rst", true}, {name, true}};
    const auto error =
        SequentialEquivalenceStrategy::findUnknownResetPort(top, reset);
    ASSERT_TRUE(error.has_value()) << name;
    EXPECT_EQ(*error, std::string("Reset bootstrap port `") + name +
                          "` was not found among aligned top-level inputs");
  }

  // run() reports it before extracting anything.
  auto* copy = load("copy", kBaseBody);
  std::vector<std::string> progress;
  auto rejected = strategy(top, copy, SecEncoding::DualRailSteady, reset);
  rejected.setProgressLogger(
      [&](SecProgressLevel, const std::string& message) {
        progress.push_back(message);
      });
  const auto result = rejected.run(4);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("arst_i"), std::string::npos) << result.reason;
  EXPECT_TRUE(progress.empty());
}

TEST_F(MemoryAbstractionTests, ProgressReportsModelSizesAndMemoryHandling) {
  auto* left = load("base", kBaseBody);
  auto* right = load("bug", kWriteDataBugBody);
  std::vector<std::pair<SecProgressLevel, std::string>> progress;
  auto logged = strategy(left, right, SecEncoding::DualRailSteady);
  logged.setProgressLogger(
      [&](SecProgressLevel level, const std::string& message) {
        progress.emplace_back(level, message);
      });
  EXPECT_EQ(logged.run(8).status, SequentialEquivalenceStatus::Different);

  const auto count = [&](SecProgressLevel level, const std::string& text) {
    return std::count_if(
        progress.begin(), progress.end(), [&](const auto& entry) {
          return entry.first == level &&
                 entry.second.find(text) != std::string::npos;
        });
  };
  EXPECT_EQ(count(SecProgressLevel::Info,
                  "comparing 1 paired memory at its interface: `mem_mem` "
                  "(16 words x 8 bits)"),
            1);
  EXPECT_EQ(count(SecProgressLevel::Warning,
                  "could not carry the proof"),
            1);
  // Two models with the memory abstracted, then two bit-level ones.
  EXPECT_EQ(count(SecProgressLevel::Info, "model has 9 state bits"), 2);
  EXPECT_EQ(count(SecProgressLevel::Info,
                  "137 state bits, 19 inputs and 9 observed outputs; 1 memory "
                  "modeled bit by bit (128 state bits)"),
            2);
  EXPECT_GE(count(SecProgressLevel::Info, "running the pdr engine"), 1);
}

TEST_F(MemoryAbstractionTests, LargeBitLevelMemoryIsReportedWithAWarning) {
  // A model that stops the run right after extraction: only its report runs.
  auto model = std::make_shared<SequentialDesignModel>();
  model->bitLevelMemories.push_back({"u_ram/mem_mem", 32, 16384});
  model->bitLevelMemories.push_back({"u_ram/small_mem", 8, 16});
  model->unsupportedReasons.push_back("stop after extraction");

  std::vector<std::pair<SecProgressLevel, std::string>> progress;
  SequentialEquivalenceStrategy reported(
      nullptr, nullptr, Config::SolverType::KISSAT, SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  reported.setProgressLogger(
      [&](SecProgressLevel level, const std::string& message) {
        progress.emplace_back(level, message);
      });
  const auto result = reported.runWithModelExtractor(
      [&](size_t, const SequentialExtractOptions&) { return model; }, 4);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.reason, "stop after extraction");

  std::vector<std::string> warnings;
  for (const auto& [level, message] : progress) {
    if (level == SecProgressLevel::Warning) {
      warnings.push_back(message);
    }
  }
  ASSERT_EQ(warnings.size(), 1u);
  EXPECT_NE(warnings.front().find("`u_ram/mem_mem`"), std::string::npos);
  EXPECT_NE(warnings.front().find("16384 words x 32 bits"), std::string::npos);
  EXPECT_NE(warnings.front().find("524288 state bits"), std::string::npos);
  EXPECT_GT(size_t{16384} * 32, kLargeBitLevelMemoryBits);
  EXPECT_LT(size_t{16} * 8, kLargeBitLevelMemoryBits);
}

// A write enable on a net nothing drives: the bit-level model reads it as a
// port that never writes, so the memory keeps its cells and the log says why.
TEST_F(MemoryAbstractionTests, UndrivenWriteEnableIsLoggedAsTheReason) {
  auto* left = load("base", kBaseBody);
  auto* right = load("copy", kBaseBody);
  for (auto* top : {left, right}) {
    auto* memory = top->getInstance(NLName("mem_mem"));
    ASSERT_NE(memory, nullptr);
    const auto interface = SNLDesignModeling::getMemoryInterface(memory);
    ASSERT_FALSE(interface.writePorts.empty());
    ASSERT_FALSE(interface.writePorts.front().enables.empty());
    memory->getInstTerm(interface.writePorts.front().enables.front())
        ->setNet(SNLScalarNet::create(top, NLName("undriven_enable")));
  }

  SequentialExtractOptions options;
  options.abstractMemories = true;
  const auto asked = SequentialDesignModel::extract(left, {}, 0, options);
  EXPECT_TRUE(asked.abstractedMemories.empty());
  ASSERT_EQ(asked.bitLevelMemories.size(), 1u);
  EXPECT_NE(asked.bitLevelMemories.front().abstractionObstacle.find(
                "write enable is unconnected, undriven"),
            std::string::npos)
      << asked.bitLevelMemories.front().abstractionObstacle;
  // Nothing to explain when the extraction was not asked to abstract.
  const auto notAsked = SequentialDesignModel::extract(left);
  ASSERT_EQ(notAsked.bitLevelMemories.size(), 1u);
  EXPECT_TRUE(notAsked.bitLevelMemories.front().abstractionObstacle.empty());

  std::vector<std::string> progress;
  auto logged = strategy(left, right, SecEncoding::DualRailSteady);
  logged.setProgressLogger(
      [&](SecProgressLevel, const std::string& message) {
        progress.push_back(message);
      });
  const auto result = logged.run(4);
  EXPECT_TRUE(result.abstractedMemories.empty());
  for (const auto* design : {"design 1", "design 2"}) {
    const std::string expected =
        std::string("SEC memory abstraction: memory `mem_mem` of ") + design +
        " keeps its cells: a write enable is unconnected";
    EXPECT_EQ(std::count_if(
                  progress.begin(), progress.end(),
                  [&](const std::string& message) {
                    return message.find(expected) != std::string::npos;
                  }),
              1)
        << design;
  }
}

// An unsupported second design ends the run: the first design, extracted
// with its memory abstracted, is not extracted a second time for nothing.
TEST_F(MemoryAbstractionTests, UnsupportedSecondDesignNeedsNoSecondExtraction) {
  auto* top = load("base", kBaseBody);
  auto unsupported = std::make_shared<SequentialDesignModel>();
  unsupported->unsupportedReasons.push_back("design 2 is unsupported");

  size_t extractions[2] = {0, 0};
  SequentialEquivalenceStrategy counted(
      nullptr, nullptr, Config::SolverType::KISSAT, SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  counted.setMemoryAbstraction(true);
  const auto result = counted.runWithModelExtractor(
      [&](size_t side, const SequentialExtractOptions& options)
          -> std::shared_ptr<const SequentialDesignModel> {
        ++extractions[side];
        if (side == 1) {
          return unsupported;
        }
        return std::make_shared<const SequentialDesignModel>(
            SequentialDesignModel::extract(top, {}, 0, options));
      },
      4);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.reason, "design 2 is unsupported");
  EXPECT_EQ(extractions[0], 1u);
  EXPECT_EQ(extractions[1], 1u);
}

// The clock of a memory is read from the netlist: the top-level inputs that
// reach its clock pin through combinational cells.
TEST_F(MemoryAbstractionTests, MemoryClockIsFoundOnTheNetlist) {
  const auto clockModule = [](const std::string& name, const char* declarations,
                              const char* writeClock) {
    return "module " + name + R"((
  input  wire       clk,
  input  wire       en,
  input  wire       rst,
  input  wire       we,
  input  wire [3:0] waddr,
  input  wire [7:0] wdata,
  input  wire [3:0] raddr,
  output reg  [7:0] rdata,
  output reg        valid
);
  reg [7:0] mem [0:15];
)" + declarations + "  always @(posedge " + writeClock + R"() begin
    if (we) mem[waddr] <= wdata;
  end
  always @(posedge clk) begin
    rdata <= mem[raddr];
    if (rst) valid <= 1'b0;
    else valid <= we;
  end
endmodule
)";
  };
  const auto sourcesOf = [](const SequentialDesignModel& model, bool abstracted) {
    return abstracted ? model.abstractedMemories.at(0).clockSources
                      : model.bitLevelMemories.at(0).clockSources;
  };
  const auto names = [](const SequentialDesignModel& model,
                        const MemoryClockSources& sources) {
    auto result = displayNames(model, sources.topInputs);
    std::sort(result.begin(), result.end());
    return result;
  };
  auto* direct = loadSource("direct", clockModule("direct", "", "clk"));
  auto* gated = loadSource(
      "gated", clockModule("gated", "  wire gclk = clk & en;\n", "gclk"));
  auto* divided = loadSource(
      "divided",
      clockModule("divided",
                  "  reg div;\n  always @(posedge clk) div <= ~div;\n",
                  "div"));

  for (const bool abstracted : {false, true}) {
    SequentialExtractOptions options;
    options.abstractMemories = abstracted;

    const auto directModel = SequentialDesignModel::extract(direct, {}, 0, options);
    const auto directSources = sourcesOf(directModel, abstracted);
    EXPECT_FALSE(directSources.internalSource);
    EXPECT_EQ(names(directModel, directSources),
              std::vector<std::string>({"clk[0]"}));

    // Through a gate: both inputs of the gate; `clk` is the clock carrier.
    const auto gatedModel = SequentialDesignModel::extract(gated, {}, 0, options);
    const auto gatedSources = sourcesOf(gatedModel, abstracted);
    EXPECT_FALSE(gatedSources.internalSource);
    EXPECT_EQ(names(gatedModel, gatedSources),
              std::vector<std::string>({"clk[0]", "en[0]"}));

    // From a register: no top-level input stands for this clock.
    const auto dividedModel =
        SequentialDesignModel::extract(divided, {}, 0, options);
    const auto dividedSources = sourcesOf(dividedModel, abstracted);
    EXPECT_TRUE(dividedSources.internalSource);
    EXPECT_TRUE(dividedSources.topInputs.empty());
  }

  // A memory on a gated version of the clock of its reader, or on a clock no
  // input stands for, is compared as before.
  for (auto* top : {gated, divided}) {
    const auto result =
        strategy(top, top, SecEncoding::DualRailSteady).run(4);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
        << top->getName().getString() << " " << result.reason;
    EXPECT_TRUE(result.multiClockDomainSkippedOutputs.empty());
  }
}

TEST_F(MemoryAbstractionTests, AbstractedLeafFollowsTheBoundaryLeafRules) {
  auto* top = load("base", kBaseBody);
  NLUniverse::get()->setTopDesign(top);
  naja::DNL::destroy();
  auto* dnl = naja::DNL::get();
  naja::DNL::DNLID memoryLeaf = naja::DNL::DNLID_MAX;
  for (const auto leafID : dnl->getLeaves()) {
    const auto& instance = dnl->getDNLInstanceFromID(leafID);
    if (instance.getSNLInstance()->getName().getString() == "mem_mem") {
      memoryLeaf = leafID;
    }
  }
  ASSERT_NE(memoryLeaf, naja::DNL::DNLID_MAX);
  EXPECT_EQ(LeafBoundary::getAbstractedLeafObstacle(*dnl, memoryLeaf), "");

  const LeafBoundary boundary(*dnl, {}, 0, {memoryLeaf});
  EXPECT_TRUE(boundary.containsInstance(memoryLeaf));
  EXPECT_TRUE(boundary.isAbstractedLeaf(memoryLeaf));
  EXPECT_FALSE(boundary.getInputs().empty());
  EXPECT_FALSE(boundary.getOutputs().empty());
  for (const auto& port : boundary.getPorts()) {
    EXPECT_TRUE(port.pathKeyed);
    EXPECT_TRUE(port.topTermName.empty());
  }

  // A leaf a user pair selects is not abstracted a second time.
  const LeafBoundary selected(*dnl, {{"mem_mem", "mem_mem"}}, 0, {memoryLeaf});
  EXPECT_TRUE(selected.containsInstance(memoryLeaf));
  EXPECT_FALSE(selected.isAbstractedLeaf(memoryLeaf));
  for (const auto& port : selected.getPorts()) {
    EXPECT_FALSE(port.pathKeyed);
  }

  // A read-data pin that shares its net with another driver breaks the rule
  // that a boundary output is the only driver of its net.
  auto* memory = top->getInstance(NLName("mem_mem"));
  ASSERT_NE(memory, nullptr);
  SNLInstTerm* readData = nullptr;
  SNLInstTerm* otherReadData = nullptr;
  for (auto* instTerm : memory->getInstTerms()) {
    if (instTerm->getBitTerm()->getDirection() != SNLTerm::Direction::Output ||
        instTerm->getNet() == nullptr) {
      continue;
    }
    (readData == nullptr ? readData : otherReadData) = instTerm;
  }
  ASSERT_NE(readData, nullptr);
  ASSERT_NE(otherReadData, nullptr);
  otherReadData->setNet(readData->getNet());
  naja::DNL::destroy();
  dnl = naja::DNL::get();
  for (const auto leafID : dnl->getLeaves()) {
    if (dnl->getDNLInstanceFromID(leafID).getSNLInstance() == memory) {
      EXPECT_NE(LeafBoundary::getAbstractedLeafObstacle(*dnl, leafID).find(
                    "only driver"),
                std::string::npos);
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    Encodings, MemoryAbstractionEncodingTests,
    ::testing::Values(SecEncoding::DualRailSteady, SecEncoding::Binary),
    [](const ::testing::TestParamInfo<SecEncoding>& info) {
      return info.param == SecEncoding::Binary ? "Binary" : "DualRail";
    });

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
