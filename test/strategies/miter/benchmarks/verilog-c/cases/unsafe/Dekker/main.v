// Model of Dekker's algorithm for mutual exclusion of two processes.
//
// Author: Fabio Somenzi <Fabio@Colorado.EDU>
// Modified by: Rajdeep Mukherjee <rajdeep.mukherjee@cs.ox.ac.uk>

module main(clock,select,pause);
    input     clock;
    input     select;
    input     pause;
    parameter L0 = 0;
    parameter L1 = 1;
    parameter L2 = 2;
    parameter L3 = 3;
    parameter L4 = 4;
    parameter L5 = 5;
    parameter L6 = 6;


    reg       c[0:1];
    reg  turn, self;
    reg [2:0]  pc[0:1];

    initial begin
        pc[0] = L0; pc[1] = L0;
        c[0] = 1; c[1] = 1;
        //turn = $ND(0,1);
   //self = select; // changed to compile with yosys
    end

    always @ (posedge clock) begin
        self = select;
        case (pc[self])
          L0: if (!pause) pc[self] = L1;        // noncritical section
          L1: begin c[self] = 0; pc[self] = L2; end
          L2: if (c[~self] == 1) pc[self] = L5; else pc[self] = L3;
          L3:
            if (turn == self)
              pc[self] = L2;
            else begin
                c[self] = 1;
                pc[self] = L4;
            end
          L4: if (turn == self) begin c[self] = 0; pc[self] = L2; end
          L5: if (!pause) pc[self] = L6;        // critical section
          L6: begin c[self] = 1; turn = ~self; pc[self] = L0; end
        endcase
    end
   // Mutual exclusion property, both the processes can not be in the critical
   // section at the same time
   assert property1: ((pc[0] == L5 && pc[1] == L5));
   // fair scheduling
   //assert property2: ((self[0] == 1 || self[0] == 0));
endmodule // dekker
