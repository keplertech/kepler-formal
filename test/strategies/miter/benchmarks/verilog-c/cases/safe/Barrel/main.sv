// This is a system composed of a register file and a shifter
// with the same number of registers of the same width.
// At each transition, the shifter rotates its contents up by one
// position.  The register file holds its contents perpetually.  The
// verification problem consists in proving the following safety
// property (an invariant):
//
//   If one entry of the register file holds the same contents of an
//   entry of the shift register, then their neighbors hold the same
//   contents too.  In formula, letting n be the number of registers
//   in the shifter and the register file:
//       r(i) = b(j) || ( r(i+1 mod n) = b(j+1 mod n),
//   for all i,j in {0,...,n-1}.
//
// Note that this invariant imposes the following constraint. Let m be
// the minimum value of |j-k| such that j != k, r(i) = b(j), and
// r(i) = b(k).  Let m = n if no j and k satisfy the requirement.  Then
// the contents of r and b are periodic with period gcd(m,n).
//
// The width of each register is the logarithm base 2 of n
// (the number of registers). This guarantees enough values to fill
// all registers of either shifter or register file with distinct values.
//
// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main(clock);
    input clock;

    reg [2-1:0] b0, r0;
    reg [2-1:0] b1, r1;
    reg [2-1:0] b2, r2;
    reg [2-1:0] b3, r3;

    integer index;

    function valid;
        input [2-1:0] b0, r0;
        input [2-1:0] b1, r1;
        input [2-1:0] b2, r2;
        input [2-1:0] b3, r3;
        begin: _valid
            valid =
                    (b0 != r0 || b1 == r1) &&
                    (b0 != r1 || b1 == r2) &&
                    (b0 != r2 || b1 == r3) &&
                    (b0 != r3 || b1 == r0) &&
                    (b1 != r0 || b2 == r1) &&
                    (b1 != r1 || b2 == r2) &&
                    (b1 != r2 || b2 == r3) &&
                    (b1 != r3 || b2 == r0) &&
                    (b2 != r0 || b3 == r1) &&
                    (b2 != r1 || b3 == r2) &&
                    (b2 != r2 || b3 == r3) &&
                    (b2 != r3 || b3 == r0) &&
                    (b3 != r0 || b0 == r1) &&
                    (b3 != r1 || b0 == r2) &&
                    (b3 != r2 || b0 == r3) &&
                    (b3 != r3 || b0 == r0);
        end
    endfunction // valid

  initial begin
        // Start by assuming all states initial.
        for (index = 0; index < 2; index = index + 1) begin
            b0[index] = 0; r0[index] = 0;
            b1[index] = 0; r1[index] = 0;
            b2[index] = 0; r2[index] = 0;
            b3[index] = 0; r3[index] = 0;
        end
        // Then remap invalid initial states to valid one.
        if (!valid(
                   b0, r0,
                   b1, r1,
                   b2, r2,
                   b3, r3
                  )) begin
            b0 = 0; r0 = 0;
            b1 = 0; r1 = 0;
            b2 = 0; r2 = 0;
            b3 = 0; r3 = 0;
        end
    end

always @ (posedge clock) begin

        b0 = b1;
        b1 = b2;
        b2 = b3;
        b3 = b0;
end

assert property
(((b0[1:0] != r0[1:0]) || (b1[1:0] == r1[1:0])) &&
((b0[1:0] != r1[1:0]) || ( b1[1:0] == r2[1:0])) &&
((b0[1:0] != r2[1:0]) || ( b1[1:0] == r3[1:0])) &&
((b0[1:0] != r3[1:0]) || ( b1[1:0] == r0[1:0])) &&
((b1[1:0] != r0[1:0]) || ( b2[1:0] == r1[1:0])) &&
((b1[1:0] != r1[1:0]) || ( b2[1:0] == r2[1:0])) &&
((b1[1:0] != r2[1:0]) || ( b2[1:0] == r3[1:0])) &&
((b1[1:0] != r3[1:0]) || ( b2[1:0] == r0[1:0])) &&
((b2[1:0] != r0[1:0]) || ( b3[1:0] == r1[1:0])) &&
((b2[1:0] != r1[1:0]) || ( b3[1:0] == r2[1:0])) &&
((b2[1:0] != r2[1:0]) || ( b3[1:0] == r3[1:0])) &&
((b2[1:0] != r3[1:0]) || ( b3[1:0] == r0[1:0])) &&
((b3[1:0] != r0[1:0]) || ( b0[1:0] == r1[1:0])) &&
((b3[1:0] != r1[1:0]) || ( b0[1:0] == r2[1:0])) &&
((b3[1:0] != r2[1:0]) || ( b0[1:0] == r3[1:0])) &&
((b3[1:0] != r3[1:0]) || ( b0[1:0] == r0[1:0])));
endmodule // barrel
