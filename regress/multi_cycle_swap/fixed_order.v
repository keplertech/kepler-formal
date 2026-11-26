// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

`define AW 4
`define W 4
// Two memory accesses per 4-phase step: save a word in the upper half, read a
// word of the lower half. SWAP changes their order only. The output is the
// same in every cycle; the upper half and `r` differ during phase 2.
// CSWAP swaps only in the steps whose mode bit m is set.
module top (
    input  wire           clk,
    input  wire           rst,
    input  wire [`AW-1:0] a_in,
    input  wire [`W-1:0]  d_in,
    output wire [`W-1:0]  o
);
    reg [1:0]     ph;
    reg [`AW-1:0] adr;                  // upper word to save into
    reg [`AW-1:0] adr2;                 // lower word to read
    reg [`W-1:0]  wd;
    reg           m;                    // mode of the step
    reg [`W-1:0]  r;
    reg [`W-1:0]  u;
    reg [`W-1:0]  mem [0:(2<<`AW)-1];   // never reset

    assign o = (ph == 2'd3) ? r : (ph == 2'd0) ? u : {`W{1'b0}};

`ifdef SWAP
    wire rd_ph = (ph == 2'd1);
    wire wr_ph = (ph == 2'd2);
`elsif CSWAP
    // the order is swapped only for the steps whose mode bit is set
    wire rd_ph = m ? (ph == 2'd1) : (ph == 2'd2);
    wire wr_ph = m ? (ph == 2'd2) : (ph == 2'd1);
`else
    wire rd_ph = (ph == 2'd2);
    wire wr_ph = (ph == 2'd1);
`endif
`ifdef BUG
    wire [`W-1:0] save = (adr == 0) ? ~wd : wd;
`else
    wire [`W-1:0] save = wd;
`endif

    always @(posedge clk) begin
        if (rst) begin
            ph <= 2'd0;
        end else begin
            ph <= ph + 2'd1;
            if (ph == 2'd0) begin
                adr  <= a_in;
                adr2 <= d_in[`AW-1:0];
                wd   <= d_in;
                m    <= d_in[`W-1];
                mem[{1'b0, a_in}] <= d_in;
            end else if (wr_ph) begin
                mem[{1'b1, adr}] <= save;
            end else if (rd_ph) begin
                r <= mem[{1'b0, adr2}];
            end else begin
                u <= mem[{1'b1, a_in}];
            end
        end
    end
endmodule
