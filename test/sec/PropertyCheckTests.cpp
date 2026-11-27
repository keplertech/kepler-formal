// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "BoolExpr.h"
#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB.h"
#include "NLLibrary.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLSVConstructor.h"
#include "SNLUtils.h"
#include "Tree2BoolExpr.h"
#include "property/PropertyCheckStrategy.h"
#include "property/PropertyExpression.h"
#include "property/SvaPropertyTranslator.h"

using namespace KEPLER_FORMAL;
using namespace KEPLER_FORMAL::SEC;
using namespace naja::NL;

namespace {

constexpr size_t kLeafVarBase = 100;
constexpr size_t kMonitorVarBase = 1000;

std::vector<BoolExpr*> makeLeaves(size_t count) {
  std::vector<BoolExpr*> leaves;
  for (size_t i = 0; i < count; ++i) {
    leaves.push_back(BoolExpr::Var(kLeafVarBase + i));
  }
  return leaves;
}

CompiledProperty compile(const std::string& expression, size_t leafCount,
                         bool cover = false) {
  size_t nextVar = kMonitorVarBase;
  PropertyCompiler compiler(makeLeaves(leafCount),
                            [&nextVar]() { return nextVar++; });
  const auto term = parsePropertyTerm(expression);
  return cover ? compiler.compileCover(term) : compiler.compileViolation(term);
}

// Runs a compiled monitor on leaf waveforms such as {"1100", "0110"}: one
// string per leaf, one character per step. Returns the signal of every step.
std::string simulate(const CompiledProperty& compiled,
                     const std::vector<std::string>& waves) {
  std::unordered_map<size_t, bool> state;
  for (const auto& bit : compiled.monitor) {
    state[bit.var] = false;
  }
  std::string signal;
  const size_t steps = waves.empty() ? 0 : waves.front().size();
  for (size_t step = 0; step < steps; ++step) {
    auto env = state;
    for (size_t leaf = 0; leaf < waves.size(); ++leaf) {
      env[kLeafVarBase + leaf] = waves[leaf][step] == '1';
    }
    signal += compiled.signal->evaluate(env) ? '1' : '0';
    std::unordered_map<size_t, bool> next;
    for (const auto& bit : compiled.monitor) {
      next[bit.var] = bit.next->evaluate(env);
    }
    state = std::move(next);
  }
  return signal;
}

std::string violations(const std::string& expression,
                       const std::vector<std::string>& waves) {
  const auto compiled = compile(expression, waves.size());
  EXPECT_TRUE(compiled.supported()) << compiled.unsupportedReason;
  return compiled.supported() ? simulate(compiled, waves) : std::string();
}

std::string unsupportedReason(const std::string& expression, size_t leafCount,
                              bool cover = false) {
  return compile(expression, leafCount, cover).unsupportedReason;
}

struct Outcome {
  PropertyStatus status = PropertyStatus::Unsupported;
  size_t bound = 0;
  std::string detail;
};

using Expected = std::map<std::string, PropertyStatus>;

constexpr PropertyStatus kProved = PropertyStatus::Proved;
constexpr PropertyStatus kFailed = PropertyStatus::Failed;
constexpr PropertyStatus kCovered = PropertyStatus::Covered;
constexpr PropertyStatus kUnreachable = PropertyStatus::Unreachable;
constexpr PropertyStatus kAssumed = PropertyStatus::Assumed;
constexpr PropertyStatus kUnsupported = PropertyStatus::Unsupported;

SecResetSpec resetBy(const std::string& port, bool activeValue, size_t cycles) {
  SecResetSpec spec;
  spec.cycles = cycles;
  spec.ports.push_back({port, activeValue});
  return spec;
}

const char* kCounter = R"sv(
// Saturating counter with properties of each supported kind.
module top(input logic clk, input logic rst, input logic en, output logic [2:0] cnt);
  always_ff @(posedge clk) begin
    if (rst) cnt <= 3'd0;
    else if (en && cnt != 3'd5) cnt <= cnt + 3'd1;
  end

  // True invariant.
  a_range: assert property (@(posedge clk) disable iff (rst) cnt <= 3'd5);
  // False invariant: the counter reaches 5.
  a_range_bad: assert property (@(posedge clk) disable iff (rst) cnt <= 3'd4);
  // True same-cycle implication.
  a_impl: assert property (@(posedge clk) disable iff (rst) (cnt == 3'd5) |-> (cnt != 3'd0));
  // True next-cycle implication.
  a_next: assert property (@(posedge clk) disable iff (rst) en |=> cnt != 3'd0);
  // False next-cycle implication.
  a_next_bad: assert property (@(posedge clk) disable iff (rst) en |=> cnt == 3'd1);
  // Reachable cover.
  c_limit: cover property (@(posedge clk) disable iff (rst) cnt == 3'd5);
  // Unreachable cover.
  c_never: cover property (@(posedge clk) disable iff (rst) cnt == 3'd7);
endmodule
)sv";

const char* kHandshake = R"sv(
// A request starts a transaction that is busy for two cycles, then pulses done.
module top(input logic clk, input logic rst_n, input logic req,
           output logic busy, output logic done);
  logic [1:0] state;  // 0 idle, 1 and 2 busy, 3 done
  always_ff @(posedge clk) begin
    if (!rst_n) state <= 2'd0;
    else case (state)
      2'd0: if (req) state <= 2'd1;
      2'd1: state <= 2'd2;
      2'd2: state <= 2'd3;
      2'd3: state <= 2'd0;
    endcase
  end
  assign busy = (state == 2'd1) || (state == 2'd2);
  assign done = (state == 2'd3);
  wire start = (state == 2'd0) && req;

  default clocking cb @(posedge clk); endclocking
  default disable iff (!rst_n);

  a_seq:        assert property (start |=> busy ##1 busy ##1 done);          // true
  a_rep:        assert property (start |=> busy[*2] ##1 done);               // true
  a_window:     assert property (start |-> ##[1:3] done);                    // true
  a_window_bad: assert property (start |-> ##[1:2] done);                    // false
  a_until:      assert property ((state == 2'd1) |-> busy until done);       // true
  a_until_bad:  assert property ((state == 2'd1) |-> (state == 2'd1) until done); // false
  a_pulse:      assert property (done |=> !done);                            // true
  a_never:      assert property (not (busy && done));                        // true
  a_ante_seq:   assert property (busy ##1 busy |=> done);                    // true
  a_ante_bad:   assert property (busy ##1 busy |=> busy);                    // false
  a_nexttime:   assert property (start |-> nexttime busy);                   // true
  a_goto:       assert property (start ##1 done[->1] |=> (state == 2'd0));   // true
  c_txn:        cover property (start ##1 busy ##1 busy ##1 done);           // covered
  c_double:     cover property (done ##1 done);                              // unreachable
  a_live:       assert property (req |-> s_eventually done);                 // true (liveness)
  a_rose:       assert property ($rose(busy) |-> ##2 done);                  // true (sampled value)
endmodule
)sv";

const char* kLiveness = R"sv(
// A divider that never stops, and a job that advances only while its
// environment does not hold it.
module top(input logic clk, input logic rst, input logic go, input logic hold,
           output logic tick, output logic done);
  logic [1:0] div;
  logic [1:0] left;
  logic busy;
  always_ff @(posedge clk) begin
    if (rst) begin
      div <= 2'd0;
      left <= 2'd0;
      busy <= 1'b0;
    end else begin
      div <= div + 2'd1;
      if (!busy) begin
        if (go) begin
          busy <= 1'b1;
          left <= 2'd2;
        end
      end else if (!hold) begin
        if (left == 2'd0) busy <= 1'b0;
        else left <= left - 2'd1;
      end
    end
  end
  assign tick = (div == 2'd3);
  assign done = busy && !hold && (left == 2'd0);

  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);

  l_tick:       assert property (s_eventually tick);                 // true: the divider never stops
  l_tick_after: assert property (go |=> s_eventually tick);          // true
  l_seq:        assert property (go ##1 !go |-> s_eventually tick);  // true, sequence on the left
  l_go_bad:     assert property (s_eventually go);                   // false: go may stay low
  l_both_bad:   assert property (go |-> s_eventually (tick && go));  // false: go may drop for good
  m_release:    assume property (hold |-> s_eventually !hold);       // a hold is always released
  l_done:       assert property (go |-> s_eventually done);          // true only under m_release
  l_idle:       assert property (busy |-> s_eventually !busy);       // true only under m_release
  l_disable:    assert property (disable iff (rst || hold) busy |-> s_eventually !busy);  // true: a hold drops the attempt
  a_left:       assert property (left <= 2'd2);                      // safety, true
  a_done_bad:   assert property (!done);                             // safety, false
endmodule
)sv";

const char* kLivenessForms = R"sv(
module top(input logic clk, input logic rst, input logic a, output logic q);
  always_ff @(posedge clk) q <= rst ? 1'b0 : a;
  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);
  m_next:  assume property (a |=> s_eventually !a);        // liveness assumption behind a delay
  l_range: assert property (s_eventually [1:3] q);         // range
  l_seq:   assert property (a |-> s_eventually (q ##1 q)); // not a Boolean
  l_until: assert property (a |-> q s_until a);            // another liveness operator
  l_q:     assert property (a |=> s_eventually q);         // true
endmodule
)sv";

const char* kAssumptions = R"sv(
// A register that loads its input; the environment promises small, spaced data.
module top(input logic clk, input logic rst, input logic [3:0] d, input logic valid,
           output logic [3:0] acc);
  logic [3:0] acc_q;
  always_ff @(posedge clk) begin
    if (rst) begin
      acc <= 4'd0;
      acc_q <= 4'd0;
    end else begin
      acc_q <= acc;
      if (valid) acc <= d;
    end
  end
  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);

  m_small:  assume property (valid |-> d <= 4'd9);          // data is at most 9
  m_gap:    assume property (valid |=> !valid);             // never valid twice in a row
  a_acc:    assert property (acc <= 4'd9);                  // true only under m_small
  a_stable: assert property (valid |=> ##1 (acc == acc_q)); // true only under m_gap
  a_bad:    assert property (acc != 4'd9);                  // false, and the trace must keep d <= 9
  c_nine:   cover property (acc == 4'd9);                   // covered
  c_big:    cover property (acc == 4'd12);                  // unreachable under m_small
endmodule
)sv";

const char* kHierarchy = R"sv(
module child(input logic clk_i, input logic rst_ni, input logic a, output logic q);
  always_ff @(posedge clk_i) begin
    if (!rst_ni) q <= 1'b0;
    else q <= a;
  end
  // Shell used by a common open-source assertion macro.
  QFollowsA_A: assert property (@(posedge clk_i) disable iff ((!rst_ni) !== '0) (a |=> q));
  // Sampled on the falling edge: not supported yet.
  Neg_A: assert property (@(negedge clk_i) disable iff (!rst_ni) (q || !q));
endmodule

module top(input logic clk_i, input logic rst_ni, input logic x, input logic y,
           output logic qx, output logic qy);
  child u_x(.clk_i, .rst_ni, .a(x), .q(qx));
  child u_y(.clk_i, .rst_ni, .a(y), .q(qy));
  logic both;
  always_comb begin
    both = qx & qy;
    Imm_A: assert (both == (qx && qy));   // immediate assertion: the frontend ignores it
  end
  Both_A: assert property (@(posedge clk_i) disable iff (!rst_ni) (x && y) |=> both);
  // No label: named after its kind and line.
  assert property (@(posedge clk_i) disable iff (!rst_ni) both |-> qx);
endmodule
)sv";

const char* kFirstCycle = R"sv(
// flag is 1 only in the first cycle after reset.
module top(input logic clk, input logic rst, input logic a, output logic flag, output logic seen);
  always_ff @(posedge clk) begin
    if (rst) begin
      flag <= 1'b1;
      seen <= 1'b0;
    end else begin
      flag <= 1'b0;
      seen <= seen | a;
    end
  end
  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);
  a_first_bad:  assert property (!flag);              // false, in the very first cycle
  a_first_ok:   assert property (flag |=> !flag);     // true
  a_input_bad:  assert property (!a);                 // false in any cycle, input only
  a_comb_ok:    assert property (a || !a);            // trivially true
  c_first:      cover property (flag);                // covered in the first cycle
  c_second:     cover property (!flag ##1 flag);      // unreachable
  a_second_bad: assert property (flag |=> flag);      // false at the second cycle
endmodule
)sv";

const char* kFirstCycleInput = R"sv(
// flag is 1 only in the first cycle after reset; these properties can only
// fail there, and only for one value of the input.
module top(input logic clk, input logic rst, input logic a, output logic flag);
  always_ff @(posedge clk) begin
    if (rst) flag <= 1'b1;
    else flag <= 1'b0;
  end
  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);
  a_first_input_bad: assert property (flag |-> !a);   // false: first cycle with a = 1
  c_first_input:     cover property (flag && a);      // covered in the first cycle
  a_first_input_ok:  assert property (flag |-> (a || !a)); // true
endmodule
)sv";

const char* kFirstCycleHistory = R"sv(
// a_q has no reset: after the reset window it holds the input sampled during
// reset, which is free. flag is 1 only in the first cycle after reset.
module top(input logic clk, input logic rst, input logic a, output logic flag, output logic a_q);
  always_ff @(posedge clk) begin
    a_q <= a;
    if (rst) flag <= 1'b1;
    else flag <= 1'b0;
  end
  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);
  a_q_bad:    assert property (flag |-> !a_q);          // false: a was 1 during reset
  a_q_bad2:   assert property (flag |-> a_q);           // false: a was 0 during reset
  c_q_one:    cover property (flag && a_q);             // covered
  c_q_zero:   cover property (flag && !a_q);            // covered
  a_past_bad: assert property ($past(a) |-> !flag);     // false, same reason through $past
  a_past_ok:  assert property ($past(a) == a_q);        // true: both sample a
  a_past2_ok: assert property (!flag |-> ($past(a, 2) == $past(a_q)));  // true from the second cycle
endmodule
)sv";

const char* kSemantics = R"sv(
// Pairs of true and false properties around the corners of each operator.
module top(input logic clk, input logic rst, input logic [1:0] v, input logic go,
           output logic [2:0] t);
  // t counts 0..4 after go, then returns to 0.
  always_ff @(posedge clk) begin
    if (rst) t <= 3'd0;
    else if (t == 3'd0) t <= go ? 3'd1 : 3'd0;
    else if (t == 3'd4) t <= 3'd0;
    else t <= t + 3'd1;
  end
  wire start = (t == 3'd0) && go;
  wire run = (t != 3'd0);

  default clocking cb @(posedge clk); endclocking
  default disable iff (rst);

  // A multi-bit operand is true when it is not zero.
  a_truthy:      assert property ((v != 2'b00) |-> v);                  // true
  a_truthy_bad:  assert property (v |-> v[0]);                          // false: v = 2
  // Consecutive repetition in a consequent.
  a_rep4:        assert property (start |=> run[*4] ##1 !run);          // true
  a_rep5_bad:    assert property (start |=> run[*5]);                   // false
  a_rep3_bad:    assert property (start |=> run[*3] ##1 !run);          // false
  // Repetition and ranges on the left of an implication.
  a_ante_rep:    assert property (run[*4] |=> !run);                    // true
  a_ante_rep_bad: assert property (run[*3] |=> !run);                   // false
  a_ante_range:  assert property (start ##[4:5] 1'b1 |-> !run || t == 3'd4); // true
  a_ante_range_bad: assert property (start ##[3:4] 1'b1 |-> !run);      // false: t == 3 or 4
  a_ante_unb:    assert property (start ##1 run[*1:$] ##1 !run |-> t == 3'd0); // true
  // Windows.
  a_win:         assert property (start |-> ##[4:6] !run);              // true: ends at 5
  a_win_exact:   assert property (start |-> ##[5:5] !run);              // true
  a_win_bad:     assert property (start |-> ##[1:4] !run);              // false
  // nexttime and always.
  a_next3:       assert property (start |-> nexttime [3] (t == 3'd3));  // true
  a_next3_bad:   assert property (start |-> nexttime [3] (t == 3'd2));  // false
  a_always_bad:  assert property (start |-> always run);                // false
  a_always:      assert property ((t == 3'd7) |-> always (t == 3'd7));  // true, never triggered
  // until.
  a_until_with:  assert property ((t == 3'd1) |-> run until_with (t == 3'd4)); // true
  a_until_with_bad: assert property ((t == 3'd1) |-> run until_with !run);     // false
  // Conjunction and if/else of properties.
  a_and:         assert property (start |=> (run and (t == 3'd1)));     // true
  a_and_bad:     assert property (start |=> (run and (t == 3'd2)));     // false
  a_if:          assert property (if (run) (t != 3'd0) else (t == 3'd0)); // true
  a_if_bad:      assert property (if (run) (t == 3'd1) else (t == 3'd0)); // false
  // goto.
  a_goto2:       assert property (start ##1 (t[0])[->2] |-> t == 3'd3); // true: odd values 1 and 3
  a_goto2_bad:   assert property (start ##1 (t[0])[->2] |-> t == 3'd1); // false
  // never.
  a_never_seq:   assert property (not (start ##5 run));                 // true: idle again after 5
  a_never_bad:   assert property (not (start ##2 run));                 // false
  // Covers.
  c_seq:         cover property (start ##4 run ##1 !run);               // covered
  c_seq_bad:     cover property (start ##5 run);                        // unreachable
  c_impl:        cover property (start |-> run);                        // unsupported: implication
  // $past with a depth, $stable, $changed, $fell.
  a_past3:       assert property ((t == 3'd4) |-> $past(t, 3) == 3'd1); // true
  a_past3_bad:   assert property ((t == 3'd4) |-> $past(t, 3) == 3'd2); // false
  a_stable:      assert property ((t == 3'd0 && !go) |=> $stable(t));   // true
  a_changed:     assert property (run |=> $changed(t));                 // true
  a_changed_bad: assert property ($changed(t) |-> run);                 // false: 4 -> 0
  a_fell:        assert property ($fell(run) |-> t == 3'd0);            // true
endmodule
)sv";

const char* kExpressions = R"sv(
// Operand expressions over free inputs: each property is a statement about
// how SystemVerilog evaluates an expression.
module top(input logic clk, input logic en,
           input logic [3:0] a, input logic [3:0] b,
           input logic signed [3:0] s, input logic [7:0] w,
           input logic [1:0] sh, output logic q);
  typedef struct packed { logic [1:0] hi; logic [1:0] lo; } pair_t;
  pair_t p;
  assign p = a;
  logic [0:3] be;       // be[0] is the most significant bit
  assign be = a;
  logic [1:0][1:0] m;   // m[1] is the upper half
  assign m = a;
  always_ff @(posedge clk) q <= en;

  default clocking cb @(posedge clk); endclocking

  // Arithmetic is modulo the operand width.
  x_add:        assert property ((a + b) == (b + a));                       // true
  x_add_bad:    assert property ((a + b) >= a);                             // false: wraps
  x_sub:        assert property (((a - b) + b) == a);                       // true
  x_neg:        assert property ((-a) == (~a + 4'd1));                      // true
  // Comparison follows the sign of the operands.
  x_unsigned:   assert property (!(a < 4'd0));                              // true
  x_signed:     assert property ((s < 0) == s[3]);                          // true
  x_signed_bad: assert property (s >= 0);                                   // false
  // A signed operand is zero-extended once the expression is unsigned.
  x_zext:       assert property ((w + s) == (w + {4'b0000, s}));            // true
  x_zext_bad:   assert property ((w + s) == (w + $signed({{4{s[3]}}, s}))); // false: s negative
  x_sext:       assert property (($signed(w) + s) == ($signed(w) + $signed({{4{s[3]}}, s}))); // true
  // Shifts, with a constant and with a variable amount.
  x_shl:        assert property ((a << 1) == {a[2:0], 1'b0});               // true
  x_shr_var:    assert property (((a >> sh) << sh) <= a);                   // true
  x_shr_bad:    assert property (((a >> sh) << sh) == a);                   // false
  x_ashr:       assert property (((s >>> 1) < 0) == (s < 0));               // true
  x_lshr_bad:   assert property (((s >> 1) < 0) == (s < 0));                // false: s negative
  // Selects.
  x_part:       assert property (a[3:2] == {a[3], a[2]});                   // true
  x_indexed:    assert property (a[1 +: 2] == a[2:1]);                      // true
  x_bit_var:    assert property (a[sh] == ((a >> sh) & 4'd1));              // true
  x_field:      assert property (p.hi == a[3:2] && p.lo == a[1:0]);         // true
  x_field_bad:  assert property (p.hi == a[1:0]);                           // false
  x_be_bit:     assert property (be[0] == a[3] && be[3] == a[0]);           // true
  x_be_part:    assert property (be[0:1] == a[3:2]);                        // true
  x_be_var:     assert property (be[sh] == a[2'd3 - sh]);                   // true
  x_md:         assert property (m[1] == a[3:2] && m[0][1] == a[1]);        // true
  x_md_var:     assert property (m[sh[0]] == (sh[0] ? a[3:2] : a[1:0]));    // true
  x_md_bad:     assert property (m[sh[0]] == a[1:0]);                       // false: sh[0] = 1
  // Conditional, inside, reductions, replication.
  x_cond:       assert property (((en ? a : b) == a) || !en);               // true
  x_inside:     assert property ((a inside {4'd1, [4'd4:4'd6]}) ==
                                 (a == 4'd1 || (a >= 4'd4 && a <= 4'd6)));  // true
  // Compared as unsigned, so the signed operand is zero-extended.
  x_inside_mix: assert property ((s inside {w}) == ({4'b0000, s} == w));    // true
  x_inside_mix_bad: assert property ((s inside {w}) ==
                                     ($signed({{4{s[3]}}, s}) == $signed(w))); // false: s negative
  x_and_red:    assert property ((&a) == (a == 4'hF));                      // true
  x_xor_red:    assert property ((^a) == (a[0] ^ a[1] ^ a[2] ^ a[3]));      // true
  x_repl:       assert property ({2{a[1:0]}} == {a[1:0], a[1:0]});          // true
  // System functions.
  x_onehot:     assert property ($onehot(a) == (a == 4'd1 || a == 4'd2 || a == 4'd4 || a == 4'd8)); // true
  x_onehot0:    assert property ($onehot0(a) == ($onehot(a) || a == 4'd0)); // true
  x_countones:  assert property ($countones(a) == (a[0] + a[1] + a[2] + a[3])); // true
  x_count_bad:  assert property ($countones(a) < 4);                        // false: a = 15
  x_known:      assert property (!$isunknown(a));                           // true: two-valued
  x_sampled:    assert property ($sampled(a) == a);                         // true
  // Not translated.
  x_mul:        assert property ((a * b) == (b * a));                       // unsupported
  x_index_wide: assert property (a[w[2:0]] || !a[w[2:0]]);                  // unsupported: index up to 7
endmodule
)sv";

}  // namespace

class PropertyExpressionTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
};

TEST_F(PropertyExpressionTests, ParsesAtomsAndLists) {
  const auto term = parsePropertyTerm(" (implies (leaf 0)\n (concat (1 $ (leaf 1))))");
  ASSERT_FALSE(term.isAtom);
  EXPECT_EQ("implies", term.op());
  ASSERT_EQ(3u, term.items.size());
  EXPECT_EQ("leaf", term.items[1].op());
  const auto& element = term.items[2].items[1];
  // A concat element starts with its delay range, not with an operator.
  ASSERT_EQ(3u, element.items.size());
  EXPECT_EQ("1", element.items[0].atom);
  EXPECT_EQ("$", element.items[1].atom);
  EXPECT_EQ("leaf", element.items[2].op());
  EXPECT_TRUE(parsePropertyTerm("x").isAtom);
}

TEST_F(PropertyExpressionTests, RejectsMalformedText) {
  for (const char* text : {"", "(", "(leaf 0", ")", "()", "(leaf 0) x", "(a (b)"}) {
    EXPECT_THROW(parsePropertyTerm(text), std::invalid_argument) << text;
  }
}

TEST_F(PropertyExpressionTests, BooleanPropertiesNeedNoMonitor) {
  for (const char* expression :
       {"(leaf 0)", "(not (leaf 0))", "(implies (leaf 0) (leaf 1))",
        "(and (leaf 0) (or (leaf 1) (leaf 0)))", "(if (leaf 0) (leaf 1) (leaf 0))",
        "(iff (leaf 0) (leaf 1))", "(concat (0 0 (leaf 0)))"}) {
    const auto compiled = compile(expression, 2);
    EXPECT_TRUE(compiled.supported()) << expression;
    EXPECT_TRUE(compiled.monitor.empty()) << expression;
  }
  // a |-> b is violated exactly where a holds and b does not.
  EXPECT_EQ("0100", violations("(implies (leaf 0) (leaf 1))", {"1100", "1011"}));
  EXPECT_EQ("0110", violations("(iff (leaf 0) (leaf 1))", {"1100", "1010"}));
  EXPECT_EQ("1001", violations("(if (leaf 0) (leaf 1) (not (leaf 1)))", {"1100", "0101"}));
}

TEST_F(PropertyExpressionTests, ImplicationDelaysAndExactSequences) {
  // a |=> b
  EXPECT_EQ("0100", violations("(implies_next (leaf 0) (leaf 1))", {"1010", "0001"}));
  // a |-> ##2 b
  EXPECT_EQ("00100", violations("(implies (leaf 0) (concat (2 2 (leaf 1))))",
                                {"10100", "00001"}));
  // a |=> b ##1 c: the attempt dies at its first failing step.
  EXPECT_EQ("010000", violations(
                          "(implies_next (leaf 0) (concat (0 0 (leaf 1)) (1 1 (leaf 2))))",
                          {"100000", "000000", "000000"}));
  EXPECT_EQ("001000", violations(
                          "(implies_next (leaf 0) (concat (0 0 (leaf 1)) (1 1 (leaf 2))))",
                          {"100000", "010000", "000000"}));
  EXPECT_EQ("000000", violations(
                          "(implies_next (leaf 0) (concat (0 0 (leaf 1)) (1 1 (leaf 2))))",
                          {"100000", "010000", "001000"}));
  // Overlapping attempts are tracked independently.
  EXPECT_EQ("00010", violations("(implies (leaf 0) (concat (2 2 (leaf 1))))",
                                {"11000", "00100"}));
  // b[*2] ##1 c in a consequent.
  const char* repetition =
      "(implies_next (leaf 0) (concat (0 0 (repeat consecutive 2 2 (leaf 1))) (1 1 (leaf 2))))";
  EXPECT_EQ("00000", violations(repetition, {"10000", "01100", "00010"}));
  EXPECT_EQ("00100", violations(repetition, {"10000", "01000", "00010"}));
  EXPECT_EQ("00010", violations(repetition, {"10000", "01100", "00000"}));
  // nexttime [2]
  EXPECT_EQ("00100", violations("(implies (leaf 0) (nexttime 2 2 (leaf 1)))",
                                {"10000", "11011"}));
}

TEST_F(PropertyExpressionTests, BoundedWindowFailsAtItsDeadline) {
  const char* window = "(implies (leaf 0) (concat (1 3 (leaf 1))))";
  EXPECT_EQ("000100", violations(window, {"100000", "000000"}));
  EXPECT_EQ("000000", violations(window, {"100000", "000100"}));
  EXPECT_EQ("000000", violations(window, {"100000", "010000"}));
  // The operand at the start step itself is outside a ##[1:3] window.
  EXPECT_EQ("000100", violations(window, {"100000", "100010"}));
  // Two attempts served by one match, and a later one that is not served.
  EXPECT_EQ("000000", violations(window, {"110000", "000100"}));
  EXPECT_EQ("000001", violations(window, {"101000", "010000"}));
  // A weak unbounded wait can never be violated.
  EXPECT_EQ("000000", violations("(implies (leaf 0) (concat (1 $ (leaf 1))))",
                                 {"100000", "000000"}));
}

TEST_F(PropertyExpressionTests, AntecedentSequencesMatchEveryAlternative) {
  // a ##1 b |-> c
  EXPECT_EQ("0100", violations("(implies (concat (0 0 (leaf 0)) (1 1 (leaf 1))) (leaf 2))",
                               {"1000", "0100", "0000"}));
  // a ##[1:2] b |-> c: both delays trigger.
  EXPECT_EQ("0110", violations("(implies (concat (0 0 (leaf 0)) (1 2 (leaf 1))) (leaf 2))",
                               {"1000", "0110", "0000"}));
  // a ##[1:$] b |-> c: every later b triggers.
  EXPECT_EQ("00101", violations("(implies (concat (0 0 (leaf 0)) (1 $ (leaf 1))) (leaf 2))",
                                {"10000", "10101", "00000"}));
  // a[*2] |=> b and a[*2:$] |=> b
  EXPECT_EQ("00110", violations("(implies_next (repeat consecutive 2 2 (leaf 0)) (leaf 1))",
                                {"11100", "00000"}));
  EXPECT_EQ("00010", violations("(implies_next (repeat consecutive 3 $ (leaf 0)) (leaf 1))",
                                {"11100", "00000"}));
  // a[->2] |-> b: only the second occurrence triggers.
  EXPECT_EQ("00010", violations("(implies (repeat goto 2 2 (leaf 0)) (leaf 1))",
                                {"01010", "00000"}));
  // not (a ##1 b): violated where the sequence matches.
  EXPECT_EQ("0100", violations("(not (concat (0 0 (leaf 0)) (1 1 (leaf 1))))",
                               {"1010", "0100"}));
}

TEST_F(PropertyExpressionTests, UntilAndAlways) {
  const char* until = "(implies (leaf 0) (until (leaf 1) (leaf 2)))";
  EXPECT_EQ("00100", violations(until, {"10000", "11000", "00000"}));
  EXPECT_EQ("00000", violations(until, {"10000", "11000", "00100"}));
  // Released in the start step: nothing is required of the left side.
  EXPECT_EQ("00000", violations(until, {"10000", "00000", "10000"}));
  const char* untilWith = "(implies (leaf 0) (until_with (leaf 1) (leaf 2)))";
  EXPECT_EQ("00100", violations(untilWith, {"10000", "11000", "00100"}));
  EXPECT_EQ("00000", violations(untilWith, {"10000", "11100", "00100"}));
  // a |-> always b: b must hold from the start step on.
  EXPECT_EQ("00011", violations("(implies (leaf 0) (always (leaf 1)))",
                                {"01000", "01100"}));
}

TEST_F(PropertyExpressionTests, CoverReportsSequenceMatches) {
  auto compiled = compile("(concat (0 0 (leaf 0)) (1 1 (leaf 1)))", 2, true);
  ASSERT_TRUE(compiled.supported());
  EXPECT_EQ("0100", simulate(compiled, {"1010", "0100"}));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(implies (leaf 0) (leaf 1))", 2, true).find("sequences only"));
}

TEST_F(PropertyExpressionTests, StrongEventuallyBecomesAnObligation) {
  // a |-> s_eventually b: no finite trace violates it.
  auto compiled = compile("(implies (leaf 0) (s_eventually (leaf 1)))", 2);
  ASSERT_TRUE(compiled.supported()) << compiled.unsupportedReason;
  EXPECT_EQ("0000", simulate(compiled, {"1100", "0000"}));
  EXPECT_TRUE(compiled.monitor.empty());
  ASSERT_EQ(1u, compiled.obligations.size());
  EXPECT_EQ(BoolExpr::Var(kLeafVarBase), compiled.obligations[0].trigger);
  EXPECT_EQ(BoolExpr::Var(kLeafVarBase + 1), compiled.obligations[0].target);
  // a |=> s_eventually b: the obligation starts one step later.
  compiled = compile("(implies_next (leaf 0) (s_eventually (leaf 1)))", 2);
  ASSERT_EQ(1u, compiled.obligations.size());
  ASSERT_EQ(1u, compiled.monitor.size());
  EXPECT_EQ(BoolExpr::Var(compiled.monitor[0].var), compiled.obligations[0].trigger);
  // A safety part next to it is still reported step by step.
  compiled = compile("(implies (leaf 0) (and (leaf 1) (s_eventually (leaf 1))))", 2);
  EXPECT_EQ("0100", simulate(compiled, {"1100", "1011"}));
  EXPECT_EQ(1u, compiled.obligations.size());
  EXPECT_TRUE(compile("(implies (leaf 0) (leaf 1))", 2).obligations.empty());

  EXPECT_NE(std::string::npos,
            unsupportedReason("(s_eventually 1 3 (leaf 0))", 1).find("range"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(s_eventually (concat (0 0 (leaf 0)) (1 1 (leaf 1))))", 2)
                .find("Boolean operand"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(s_eventually (leaf 0))", 1, true).find("sequences only"));
}

TEST_F(PropertyExpressionTests, ReportsWhatItCannotCompile) {
  EXPECT_NE(std::string::npos,
            unsupportedReason("(implies (leaf 0) (s_until (leaf 1) (leaf 0)))", 2)
                .find("liveness"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(implies (leaf 0) (strong (concat (1 1 (leaf 1)))))", 2)
                .find("liveness"));
  EXPECT_NE(std::string::npos,
            unsupportedReason(
                "(implies (leaf 0) (concat (1 2 (leaf 1)) (1 1 (leaf 0))))", 2)
                .find("alternative tracking"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(implies (intersect (leaf 0) (leaf 1)) (leaf 0))", 2)
                .find("intersect"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(implies (repeat consecutive 0 2 (leaf 0)) (leaf 1))", 2)
                .find("empty match"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(or (implies_next (leaf 0) (leaf 1)) (leaf 0))", 2)
                .find("Boolean operands"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(not (implies_next (leaf 0) (leaf 1)))", 2).find("negation"));
  EXPECT_NE(std::string::npos, unsupportedReason("(leaf 7)", 2).find("leaf index"));
  EXPECT_NE(std::string::npos,
            unsupportedReason("(implies (leaf 0) (concat (5000 5000 (leaf 1))))", 2)
                .find("4096"));
  EXPECT_FALSE(unsupportedReason("(frobnicate (leaf 0))", 1).empty());
}

class PropertyCheckTests : public ::testing::Test {
 protected:
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) {
      universe->destroy();
    }
    Tree2BoolExpr::iso2boolExpr_.clear();
    BoolExprCache::destroy();
  }

  // Loads the SystemVerilog text with its assertions kept and checks them.
  std::map<std::string, Outcome> check(const std::string& text,
                                       const SecResetSpec& reset,
                                       SecEngine engine = SecEngine::Pdr,
                                       std::string* error = nullptr) {
    TearDown();
    const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
    const auto path = std::filesystem::temp_directory_path() /
                      ("kepler_fpv_" + std::string(info->name()) + "_" +
                       std::to_string(::getpid()) + "_" +
                       std::to_string(counter_++) + ".sv");
    {
      std::ofstream file(path);
      file << text;
    }
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* library = NLLibrary::create(db, NLName("DESIGN"));
    SNLSVConstructor constructor(library);
    SNLSVConstructor::ConstructOptions options;
    options.preserveAssertions = true;
    options.keepASTLink = true;
    options.diagnosticsReportPath.reset();
    try {
      constructor.construct(path, options);
    } catch (...) {
      std::filesystem::remove(path);
      throw;
    }
    std::filesystem::remove(path);
    auto* top = SNLUtils::findTop(library);
    EXPECT_NE(top, nullptr);
    db->setTopDesign(top);

    PropertyCheckStrategy strategy(top, makeSvaPropertyTranslator(),
                                   Config::SolverType::KISSAT, engine, reset);
    const auto result = strategy.run(32);
    if (error != nullptr) {
      *error = result.error;
    } else {
      EXPECT_TRUE(result.error.empty()) << result.error;
    }
    std::map<std::string, Outcome> outcomes;
    for (const auto& property : result.properties) {
      outcomes[property.name] = {property.status, property.bound, property.detail};
    }
    return outcomes;
  }

  static void expectStatuses(const std::map<std::string, Outcome>& outcomes,
                             const Expected& expected) {
    EXPECT_EQ(expected.size(), outcomes.size());
    for (const auto& [name, status] : expected) {
      const auto found = outcomes.find(name);
      if (found == outcomes.end()) {
        ADD_FAILURE() << "missing property " << name;
        continue;
      }
      EXPECT_STREQ(propertyStatusName(status),
                   propertyStatusName(found->second.status))
          << name << ": " << found->second.detail;
    }
  }

  size_t counter_ = 0;
};

TEST_F(PropertyCheckTests, CounterInvariantsImplicationsAndCovers) {
  const auto outcomes = check(kCounter, resetBy("rst", true, 1));
  expectStatuses(outcomes, {{"a_range", kProved},
                            {"a_range_bad", kFailed},
                            {"a_impl", kProved},
                            {"a_next", kProved},
                            {"a_next_bad", kFailed},
                            {"c_limit", kCovered},
                            {"c_never", kUnreachable}});
  // Steps are counted from the first cycle after reset.
  EXPECT_EQ(5u, outcomes.at("a_range_bad").bound);
  EXPECT_EQ(2u, outcomes.at("a_next_bad").bound);
  EXPECT_EQ(5u, outcomes.at("c_limit").bound);
  const auto& trace = outcomes.at("a_range_bad").detail;
  EXPECT_NE(std::string::npos, trace.find("  reset cycle 0:")) << trace;
  EXPECT_NE(std::string::npos, trace.find("  cycle 4: en[0]=1")) << trace;
  // The reset port is driven by the reset sequence and the clock is implicit.
  EXPECT_EQ(std::string::npos, trace.find("rst")) << trace;
  EXPECT_EQ(std::string::npos, trace.find("clk")) << trace;
}

TEST_F(PropertyCheckTests, SequencesWindowsAndSampledValues) {
  expectStatuses(check(kHandshake, resetBy("rst_n", false, 2)),
                 {{"a_seq", kProved},
                  {"a_rep", kProved},
                  {"a_window", kProved},
                  {"a_window_bad", kFailed},
                  {"a_until", kProved},
                  {"a_until_bad", kFailed},
                  {"a_pulse", kProved},
                  {"a_never", kProved},
                  {"a_ante_seq", kProved},
                  {"a_ante_bad", kFailed},
                  {"a_nexttime", kProved},
                  {"a_goto", kProved},
                  {"c_txn", kCovered},
                  {"c_double", kUnreachable},
                  {"a_live", kProved},
                  {"a_rose", kProved}});
}

TEST_F(PropertyCheckTests, OperatorCorners) {
  expectStatuses(check(kSemantics, resetBy("rst", true, 4)),
                 {{"a_truthy", kProved},
                  {"a_truthy_bad", kFailed},
                  {"a_rep4", kProved},
                  {"a_rep5_bad", kFailed},
                  {"a_rep3_bad", kFailed},
                  {"a_ante_rep", kProved},
                  {"a_ante_rep_bad", kFailed},
                  {"a_ante_range", kProved},
                  {"a_ante_range_bad", kFailed},
                  {"a_ante_unb", kProved},
                  {"a_win", kProved},
                  {"a_win_exact", kProved},
                  {"a_win_bad", kFailed},
                  {"a_next3", kProved},
                  {"a_next3_bad", kFailed},
                  {"a_always_bad", kFailed},
                  {"a_always", kProved},
                  {"a_until_with", kProved},
                  {"a_until_with_bad", kFailed},
                  {"a_and", kProved},
                  {"a_and_bad", kFailed},
                  {"a_if", kProved},
                  {"a_if_bad", kFailed},
                  {"a_goto2", kProved},
                  {"a_goto2_bad", kFailed},
                  {"a_never_seq", kProved},
                  {"a_never_bad", kFailed},
                  {"c_seq", kCovered},
                  {"c_seq_bad", kUnreachable},
                  {"c_impl", kUnsupported},
                  {"a_past3", kProved},
                  {"a_past3_bad", kFailed},
                  {"a_stable", kProved},
                  {"a_changed", kProved},
                  {"a_changed_bad", kFailed},
                  {"a_fell", kProved}});
}

TEST_F(PropertyCheckTests, AssumptionsConstrainEveryTrace) {
  const auto outcomes = check(kAssumptions, resetBy("rst", true, 1));
  expectStatuses(outcomes, {{"m_small", kAssumed},
                            {"m_gap", kAssumed},
                            {"a_acc", kProved},
                            {"a_stable", kProved},
                            {"a_bad", kFailed},
                            {"c_nine", kCovered},
                            {"c_big", kUnreachable}});
  // The failing trace loads 9, the largest value the assumption allows.
  const auto& trace = outcomes.at("a_bad").detail;
  EXPECT_NE(std::string::npos, trace.find("d[0]=1, d[1]=0, d[2]=0, d[3]=1")) << trace;

  // Without the assumptions the same assertions fail.
  std::string unconstrained = kAssumptions;
  for (const char* keyword : {"m_small:  assume", "m_gap:    assume"}) {
    const auto position = unconstrained.find(keyword);
    ASSERT_NE(std::string::npos, position);
    unconstrained.replace(position, std::string(keyword).size(),
                          std::string(keyword).substr(0, 10) + "cover ");
  }
  const auto free = check(unconstrained, resetBy("rst", true, 1));
  EXPECT_EQ(kFailed, free.at("a_acc").status);
  EXPECT_EQ(kFailed, free.at("a_stable").status);
  EXPECT_EQ(kCovered, free.at("c_big").status);
}

TEST_F(PropertyCheckTests, LivenessObligationsAndAssumptions) {
  const auto outcomes = check(kLiveness, resetBy("rst", true, 1));
  expectStatuses(outcomes, {{"l_tick", kProved},
                            {"l_tick_after", kProved},
                            {"l_seq", kProved},
                            {"l_go_bad", kFailed},
                            {"l_both_bad", kFailed},
                            {"m_release", kAssumed},
                            {"l_done", kProved},
                            {"l_idle", kProved},
                            {"l_disable", kProved},
                            {"a_left", kProved},
                            {"a_done_bad", kFailed}});
  // A liveness failure is a trace that ends in a loop where go stays low.
  const auto& loop = outcomes.at("l_go_bad").detail;
  EXPECT_NE(std::string::npos, loop.find("Loop: cycle ")) << loop;
  EXPECT_NE(std::string::npos, loop.find("go[0]=0")) << loop;
  // The variables added for the check are not part of the trace.
  EXPECT_EQ(std::string::npos, loop.find("__kepler")) << loop;
  EXPECT_GE(outcomes.at("l_go_bad").bound, 1u);
  // A liveness assumption says nothing about a finite trace.
  EXPECT_NE(std::string::npos,
            outcomes.at("a_done_bad").detail.find("liveness assumption"));
  EXPECT_EQ(std::string::npos, outcomes.at("a_done_bad").detail.find("Loop:"));

  // Without the assumption the job may be held forever.
  std::string unconstrained = kLiveness;
  const std::string keyword = "m_release:    assume";
  const auto position = unconstrained.find(keyword);
  ASSERT_NE(std::string::npos, position);
  unconstrained.replace(position, keyword.size(), "m_release:    cover ");
  const auto free = check(unconstrained, resetBy("rst", true, 1));
  EXPECT_EQ(kFailed, free.at("l_done").status);
  EXPECT_EQ(kFailed, free.at("l_idle").status);
  EXPECT_NE(std::string::npos, free.at("l_done").detail.find("hold[0]=1"));
  EXPECT_EQ(kProved, free.at("l_disable").status);
  EXPECT_EQ(kProved, free.at("l_tick").status);
  EXPECT_EQ(kUnsupported, free.at("m_release").status);
}

TEST_F(PropertyCheckTests, InterpolationEngineFindsLivenessFailures) {
  const auto outcomes = check(kLiveness, resetBy("rst", true, 1), SecEngine::Imc);
  EXPECT_EQ(kFailed, outcomes.at("l_go_bad").status);
  EXPECT_EQ(kFailed, outcomes.at("l_both_bad").status);
  EXPECT_NE(std::string::npos, outcomes.at("l_go_bad").detail.find("Loop: cycle "));
  // A true liveness property is proved or left open, never failed.
  for (const char* name :
       {"l_tick", "l_tick_after", "l_seq", "l_done", "l_idle", "l_disable"}) {
    const auto actual = outcomes.at(name).status;
    EXPECT_TRUE(actual == kProved || actual == PropertyStatus::Inconclusive)
        << name << ": " << propertyStatusName(actual);
  }
}

TEST_F(PropertyCheckTests, LivenessFormsOutsideTheSubset) {
  const auto outcomes = check(kLivenessForms, resetBy("rst", true, 1));
  expectStatuses(outcomes, {{"m_next", kUnsupported},
                            {"l_range", kUnsupported},
                            {"l_seq", kUnsupported},
                            {"l_until", kUnsupported},
                            {"l_q", kProved}});
  EXPECT_NE(std::string::npos, outcomes.at("m_next").detail.find("liveness assumption"));
  EXPECT_NE(std::string::npos, outcomes.at("l_range").detail.find("range"));
  EXPECT_NE(std::string::npos, outcomes.at("l_seq").detail.find("Boolean operand"));
  EXPECT_NE(std::string::npos, outcomes.at("l_until").detail.find("s_until"));
}

TEST_F(PropertyCheckTests, HierarchyAndUnsupportedForms) {
  // An assertion without a label is named after its kind and line.
  const std::string text = kHierarchy;
  const auto position = text.find("  assert property");
  ASSERT_NE(std::string::npos, position);
  const auto line = 1 + std::count(text.begin(), text.begin() + position, '\n');
  const auto outcomes = check(kHierarchy, resetBy("rst_ni", false, 1));
  expectStatuses(outcomes, {{"u_x/QFollowsA_A", kProved},
                            {"u_x/Neg_A", kUnsupported},
                            {"u_y/QFollowsA_A", kProved},
                            {"u_y/Neg_A", kUnsupported},
                            {"Both_A", kProved},
                            {"assert_line" + std::to_string(line), kProved}});
  EXPECT_NE(std::string::npos, outcomes.at("u_x/Neg_A").detail.find("falling clock edge"));
}

// Operands are compiled by SystemVerilog's width and sign rules.
TEST_F(PropertyCheckTests, OperandExpressions) {
  const auto outcomes = check(kExpressions, SecResetSpec{});
  expectStatuses(outcomes, {{"x_add", kProved},
                            {"x_add_bad", kFailed},
                            {"x_sub", kProved},
                            {"x_neg", kProved},
                            {"x_unsigned", kProved},
                            {"x_signed", kProved},
                            {"x_signed_bad", kFailed},
                            {"x_zext", kProved},
                            {"x_zext_bad", kFailed},
                            {"x_sext", kProved},
                            {"x_shl", kProved},
                            {"x_shr_var", kProved},
                            {"x_shr_bad", kFailed},
                            {"x_ashr", kProved},
                            {"x_lshr_bad", kFailed},
                            {"x_part", kProved},
                            {"x_indexed", kProved},
                            {"x_bit_var", kProved},
                            {"x_field", kProved},
                            {"x_field_bad", kFailed},
                            {"x_be_bit", kProved},
                            {"x_be_part", kProved},
                            {"x_be_var", kProved},
                            {"x_md", kProved},
                            {"x_md_var", kProved},
                            {"x_md_bad", kFailed},
                            {"x_cond", kProved},
                            {"x_inside", kProved},
                            {"x_inside_mix", kProved},
                            {"x_inside_mix_bad", kFailed},
                            {"x_and_red", kProved},
                            {"x_xor_red", kProved},
                            {"x_repl", kProved},
                            {"x_onehot", kProved},
                            {"x_onehot0", kProved},
                            {"x_countones", kProved},
                            {"x_count_bad", kFailed},
                            {"x_known", kProved},
                            {"x_sampled", kProved},
                            {"x_mul", kUnsupported},
                            {"x_index_wide", kUnsupported}});
  EXPECT_NE(std::string::npos, outcomes.at("x_mul").detail.find("multiplication"));
  EXPECT_NE(std::string::npos,
            outcomes.at("x_index_wide").detail.find("outside its vector"));
}

// The first cycle after reset is checked like any other: nothing is assumed
// about it, whatever the registers hold then.
TEST_F(PropertyCheckTests, FirstCycleAfterResetIsChecked) {
  for (const auto engine : {SecEngine::Pdr, SecEngine::KInduction}) {
    const auto first = check(kFirstCycle, resetBy("rst", true, 1), engine);
    expectStatuses(first, {{"a_first_bad", kFailed},
                           {"a_first_ok", kProved},
                           {"a_input_bad", kFailed},
                           {"a_comb_ok", kProved},
                           {"c_first", kCovered},
                           {"c_second", kUnreachable},
                           {"a_second_bad", kFailed}});
    EXPECT_EQ(0u, first.at("a_first_bad").bound);
    EXPECT_EQ(0u, first.at("a_input_bad").bound);
    EXPECT_EQ(1u, first.at("a_second_bad").bound);

    expectStatuses(check(kFirstCycleInput, resetBy("rst", true, 2), engine),
                   {{"a_first_input_bad", kFailed},
                    {"c_first_input", kCovered},
                    {"a_first_input_ok", kProved}});

    // Registers without reset keep what they sampled during the reset window.
    const auto history = check(kFirstCycleHistory, resetBy("rst", true, 2), engine);
    expectStatuses(history, {{"a_q_bad", kFailed},
                             {"a_q_bad2", kFailed},
                             {"c_q_one", kCovered},
                             {"c_q_zero", kCovered},
                             {"a_past_bad", kFailed},
                             {"a_past_ok", kProved},
                             {"a_past2_ok", kProved}});
    const auto& trace = history.at("a_past_bad").detail;
    EXPECT_NE(std::string::npos, trace.find("initial ")) << trace;
  }
}

TEST_F(PropertyCheckTests, WithoutResetEveryRegisterStartsFree) {
  const auto outcomes = check(kCounter, SecResetSpec{});
  // From an arbitrary start the counter may already be above its limit.
  EXPECT_EQ(kFailed, outcomes.at("a_range").status);
  EXPECT_EQ(0u, outcomes.at("a_range").bound);
  EXPECT_EQ(kCovered, outcomes.at("c_never").status);
  EXPECT_EQ(kProved, outcomes.at("a_impl").status);
}

TEST_F(PropertyCheckTests, InterpolationEngineNeverContradicts) {
  const auto outcomes = check(kCounter, resetBy("rst", true, 1), SecEngine::Imc);
  const Expected expected = {{"a_range", kProved},
                             {"a_range_bad", kFailed},
                             {"a_impl", kProved},
                             {"a_next", kProved},
                             {"a_next_bad", kFailed},
                             {"c_limit", kCovered},
                             {"c_never", kUnreachable}};
  for (const auto& [name, status] : expected) {
    const auto actual = outcomes.at(name).status;
    EXPECT_TRUE(actual == status || actual == PropertyStatus::Inconclusive)
        << name << ": " << propertyStatusName(actual);
  }
}

TEST_F(PropertyCheckTests, ReportsDesignLevelProblems) {
  std::string error;
  auto outcomes = check(kCounter, resetBy("no_such_port", true, 1), SecEngine::Pdr, &error);
  EXPECT_TRUE(outcomes.empty());
  EXPECT_NE(std::string::npos, error.find("no_such_port")) << error;

  outcomes = check(kCounter, resetBy("rst", true, 100000), SecEngine::Pdr, &error);
  EXPECT_NE(std::string::npos, error.find("reset window")) << error;

  // A design without properties gives an empty, error-free result.
  outcomes = check("module top(input logic a, output logic b); assign b = a; endmodule",
                   SecResetSpec{}, SecEngine::Pdr, &error);
  EXPECT_TRUE(outcomes.empty());
  EXPECT_TRUE(error.empty()) << error;
}
