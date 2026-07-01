// Verilog model for the round-robin arbiter described in
// the CHARME99 paper by Katz, Grumberg, and Geist.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main(clock,ir0,ir1,ack0,ack1);
    input  clock;
    input  ir0, ir1;
    output ack0, ack1;

    reg    req0, req1, ack0, ack1, robin;

    initial begin
        ack0 = 0; ack1 = 0; robin = 0;
        req0 = ir0; req1 = ir1; // nondeterministic initial requests
    end

    always @ (posedge clock) begin
        if (~req0)
          ack0 = 0;             // no request -> no ack
        else if (~req1)
          ack0 = 1;             // a single request
        else if (~ack0 & ~ack1)
          ack0 = ~robin;        // simultaneous request assertions
        else
          ack0 = ~ack0;         // both requesting: toggle ack
    end

    always @ (posedge clock) begin
        if (~req1)
          ack1 = 0;             // no request -> no ack
        else if (~req0)
          ack1 = 1;             // a single request
        else if (~ack0 & ~ack1)
          ack1 = robin;         // simultaneous request assertions
        else
          ack1 = ~ack1;         // both requesting: toggle ack
    end

    always @ (posedge clock) begin
        if (req0 & req1 & ~ack0 & ~ack1)
          robin = ~robin;       // simultaneous request assertions
    end

    // Latched inputs.
    always @ (posedge clock) begin
        req0 = ir0;
        req1 = ir1;
    end

  // Mutual exclusion.
  assert property1: (ack0==0 || ack1==0);
  assert property (req0==0 && req1==0 |-> ##1 (ack0==0 && ack1==0));
  assert property (req0==1 && req1==0 |-> ##1 ack0==1);
  assert property (req0==0 && req1==1 |-> ##1 ack1==1);
  assert property (req0==1 && ack1==1 |-> ##1 ack0==1);
  assert property (req1==1 && ack0==1 |-> ##1 ack1==1);

endmodule // rrobin
