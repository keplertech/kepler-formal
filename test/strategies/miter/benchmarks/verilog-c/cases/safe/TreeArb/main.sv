// A simple synchronous tree arbiter.
//
// The arbiter is made up of identical cells.  Each cell---except the
// root cell---acts toward its parent like a client process, asserting
// req0 to issue a request on behalf of its children, and awaiting
// ack0 before acknowledging its children's requests.  Suppose that
// starting from a quiescent condition, req1 of C0_0 is asserted
// and req2 is not (P0 is requesting and P1 is not).  Then, the
// arbiter cell relays the request to its parent C1_0 by asserting
// req0 in the next clock cycle.  When C1_0 asserts its ack1
// output (which is the ack0 input of C0_0) C0_0 asserts its
// own ack1 output (with one clock cycle delay) and P0 is allowed
// to proceed.  When P0 returns to the idle state, C0_0 lowers
// its req0 output to let C1_0 know that the arbitration cycle is
// over.  (Once again, this takes one clock cycle.)

// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main (clk, choice);
    input clk;
    input [7:0] choice;

    wire  [0:0] a3;
    wire  [1:0] a2;
    wire  [3:0] a1;
    wire  [7:0] a0;
    wire  [0:0] r3;
    wire  [1:0] r2;
    wire  [3:0] r1;
    wire  [7:0] r0;

    assign a3[0] = r3[0];       // The root's request is always acknowledged

    cell C2_0 (clk, r2[0], r2[1], a2[0], a2[1], r3[0], a3[0]);

    cell C1_0 (clk, r1[0], r1[1], a1[0], a1[1], r2[0], a2[0]);
    cell C1_1 (clk, r1[2], r1[3], a1[2], a1[3], r2[1], a2[1]);

    cell C0_0 (clk, r0[0], r0[1], a0[0], a0[1], r1[0], a1[0]);
    cell C0_1 (clk, r0[2], r0[3], a0[2], a0[3], r1[1], a1[1]);
    cell C0_2 (clk, r0[4], r0[5], a0[4], a0[5], r1[2], a1[2]);
    cell C0_3 (clk, r0[6], r0[7], a0[6], a0[7], r1[3], a1[3]);

    proc P0 (clk, a0[0], r0[0], choice[0]);
    proc P1 (clk, a0[1], r0[1], choice[1]);
    proc P2 (clk, a0[2], r0[2], choice[2]);
    proc P3 (clk, a0[3], r0[3], choice[3]);
    proc P4 (clk, a0[4], r0[4], choice[4]);
    proc P5 (clk, a0[5], r0[5], choice[5]);
    proc P6 (clk, a0[6], r0[6], choice[6]);
    proc P7 (clk, a0[7], r0[7], choice[7]);

 // #PASS: Mutual exclusion.
 assert property (!((P0.state==2 || P1.state==2 || P2.state==2 || P3.state==2) &&
 (P4.state==2 || P5.state==2 || P6.state==2 || P7.state==2) ||
 (P0.state==2 || P1.state==2 || P4.state==2 || P5.state==2) &&
 (P2.state==2 || P3.state==2 || P6.state==2 || P7.state==2) ||
 (P0.state==2 || P2.state==2 || P4.state==2 || P6.state==2) &&
 (P1.state==2 || P3.state==2 || P5.state==2 || P7.state==2)));


 // #PASS: Absence of starvation.  Needs fairness constraints.
 assert property (P0.state==1 |-> ##[0:$] P0.state==2);
 assert property (P1.state==1 |-> ##[0:$] P1.state==2);
 assert property (P2.state==1 |-> ##[0:$] P2.state==2);
 assert property (P3.state==1 |-> ##[0:$] P3.state==2);
 assert property (P4.state==1 |-> ##[0:$] P4.state==2);
 assert property (P5.state==1 |-> ##[0:$] P5.state==2);
 assert property (P6.state==1 |-> ##[0:$] P6.state==2);
 assert property (P7.state==1 |-> ##[0:$] P7.state==2);

endmodule // arbiter


// A client process loops through three states: idle, requesting, and
// locking.  The transitions from idle to requesting, and from locking
// to idle are nondeterministic.

//typedef enum {idle, requesting, locking} ProcPhase;

module proc (clk, ack, req, choice);
    input  clk;
    input  ack;
    input  choice;
    output req;

    parameter idle = 0;
    parameter requesting = 1;
    parameter locking = 2;

    reg state;

    initial state = idle;

    always @ (posedge clk) begin
        case (state)
          idle:       if (choice) state = requesting;
          requesting: if (ack)    state = locking;
          locking:    if (choice) state = idle;
        endcase
    end

    assign req = state == requesting || state == locking;

endmodule // proc


// The arbiter cell is a Moore machine with two inputs from children and
// two outputs to chidren; one input from parent, and one output to parent.
// If both children of a cell issue a request in the same clock cycle,
// the cell chooses one of the children, giving priority to the one
// that did not receive the last grant.  Before any grant is
// given---for instance in the initial state---the child connected to
// req1 and ack1 is given priority in case of simultaneous requests.

//typedef enum {I1, I2, R1, R2, A1, A2} CellState;

module cell (clk, req1, req2, ack1, ack2, req0, ack0);
    input  clk;
    input  ack0, req1, req2;
    output req0, ack1, ack2;

    parameter I1 = 0;
    parameter I2 = 1;
    parameter R1 = 2;
    parameter R2 = 3;
    parameter A1 = 4;
    parameter A2 = 5;

    reg state;
    initial state = I1;

    always @ (posedge clk) begin
        case (state)
          I1: if (req1)  state = R1; else if (req2) state = R2;
          R1: if (ack0)  state = A1;
          A1: if (!req1) state = I2;
          I2: if (req2)  state = R2; else if (req1) state = R1;
          R2: if (ack0)  state = A2;
          A2: if (!req2) state = I1;
        endcase
    end

    assign req0 = state == R1 || state == R2 || state == A1 || state == A2;
    assign ack1 = state == A1;
    assign ack2 = state == A2;

endmodule // cell
