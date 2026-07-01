// Model for the RCU (read-copy update) mutual exclusion mechanism.
// Translated from the Promela model of Paul McKenney.
// This model assumes one update process.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>

//typedef enum {L0, L1, L2, L3, L4, L5, L6, L7} locR;
//typedef enum {L0, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10} locU;

// Modified by: Rajdeep Mukherjee <rajdeep.mukherjee@cs.ox.ac.uk>
module main(clock, select);
                                // 0 <= select < NRDR --> reader
                                // else               --> update

    parameter           PASSES = 10;
    parameter           NRDR = 4;       // number of reader processes
    parameter           NRDR_ELEM = NRDR + NRDR;
    parameter           SELMSB = 2;     // MSB for the select input
                                        // must satisfy 2^(SELMSB+1) > NRDR


    parameter L0 = 0;
    parameter L1 = 1;
    parameter L2 = 2;
    parameter L3 = 3;
    parameter L4 = 4;
    parameter L5 = 5;
    parameter L6 = 6;
    parameter L7 = 7;
    parameter L8 = 8;
    parameter L9 = 9;
    parameter L10 = 10;

    input               clock;
    input [SELMSB:0]    select; // nondeterministic scheduler
    reg                 flip;
    reg                 ctr[0:NRDR_ELEM-1];
    reg [7:0]           passctr;
    // latched version of select, to which we can refer in properties
    reg [SELMSB:0]      self;
    // local variables for the reader processes
    reg [2:0]   pc [0:NRDR-1];
    reg                 lclFlip [0:NRDR-1];
    reg                 both [0:NRDR-1];
    // local variables for the update process
    reg [3:0]   pcu;
    reg [7:0]           lclPassctr;
    reg [SELMSB:0]      cpunum;

    integer             i;

    initial begin
        flip = 0;
        passctr = 0;
        cpunum = 0;
        for (i = 0; i < NRDR_ELEM; i = i + 1)
          ctr[i] = 0;
        for (i = 0; i < NRDR; i = i + 1) begin
            lclFlip[i] = 0;
            both[i] = 0;
            pc[i] = L0;
        end
        pcu = L0;
        self = 0;
        lclPassctr = 0;
    end

    always @ (posedge clock) begin
        self = select;
        if (self >= NRDR) begin
            // upd process
            case (pcu)
              L0: if (passctr < PASSES) begin
                  lclPassctr = passctr;
                  pcu = L1;
              end
              L1: begin
                  if (~lclPassctr[0])
                    lclPassctr = 255;
                  pcu = L2;
              end
              L2: begin
                  cpunum = 0;
                  pcu = L3;
              end
              L3: if (cpunum < NRDR)
                pcu = L4;
              else
                pcu = L6;
              L4: if (ctr[{cpunum,~flip}] == 0)
                pcu = L5;
              L5: begin
                  cpunum = cpunum + 1;
                  pcu = L3;
              end
              L6: begin
                  flip = ~flip;
                  pcu = L7;
              end
              L7: begin
                  cpunum = 0;
                  pcu = L8;
              end
              L8: if (cpunum < NRDR)
                pcu = L9;
              else
                pcu = L0;
              L9: if (ctr[{cpunum,~flip}] == 0)
                pcu = L10;
              L10: begin
                  cpunum = cpunum + 1;
                  pcu = L8;
              end
            endcase
        end else begin
            // rdr process
            case (pc[self])
              L0: if (passctr < PASSES) begin
                  lclFlip[self] = flip;
                  pc[self] = L1;
              end
              L1: begin
                  ctr[{self,lclFlip[self]}] = ~ctr[{self,lclFlip[self]}];
                  pc[self] = L2;
              end
              L2: if (lclFlip[self] == flip) begin
                  both[self] = 0;
                  pc[self] = L4;
              end else begin
                  ctr[{self,~lclFlip[self]}] = ~ctr[{self,~lclFlip[self]}];
                  pc[self] = L3;
              end
              L3: begin
                  both[self] = 1;
                  pc[self] = L4;
              end
              L4: begin
                  passctr = passctr + 1;
                  pc[self] = L5;
              end
              L5: begin
                  passctr = passctr + 1;
                  pc[self] = L6;
              end
              L6: begin
                  ctr[{self,lclFlip[self]}] = ~ctr[{self,lclFlip[self]}];
                  pc[self] = L7;
              end
              L7: begin
                  if (both[self])
                   ctr[{self,~lclFlip[self]}] = ~ctr[{self,~lclFlip[self]}];
                  pc[self] = L0;
              end
            endcase
        end
    end

// Safety Invariant
assert property1: (!(pcu==L8 && cpunum[2:0]==4)) || (!(lclPassctr[7:0] == passctr[7:0]));
assert property2: (lclPassctr[7:4]==0 || lclPassctr[7:4] == 15);
//assert property3: (passctr[7:4]==0);
endmodule // rcu
