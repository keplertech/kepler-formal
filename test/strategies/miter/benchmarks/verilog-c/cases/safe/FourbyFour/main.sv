// Slider puzzle reduced to 2 rows and four columns. +----+----+----+----+
//                                                   |  0 |  1 |  2 |  3 |
// The entries of the matrix are numbered thus:      +----+----+----+----+
//                                                   |  4 |  5 |  6 |  7 |
// Author: Fabio Somenzi <Fabio@Colorado.EDU>        +----+----+----+----+

module main(clock,from,to);
    input       clock;
    input [2:0] from;
    input [2:0] to;

    reg [2:0]   b[0:7];
    reg [2:0]   freg, treg;
    wire        valid, parity, permutation, oddInversions;

    // The initial states include illegal configurations, that is,
    // configurations that are not permutations of (0,...,7).  We use the
    // "permutation" predicate in the properties to restrict considetation
    // to legal initial states.
    initial begin
/*      b[0] = $ND(0,1,2,3,4,5,6,7);
        b[1] = $ND(0,1,2,3,4,5,6,7);
        b[2] = $ND(0,1,2,3,4,5,6,7);
        b[3] = $ND(0,1,2,3,4,5,6,7);
        b[4] = $ND(0,1,2,3,4,5,6,7);
        b[5] = $ND(0,1,2,3,4,5,6,7);
        b[6] = $ND(0,1,2,3,4,5,6,7);
        b[7] = $ND(0,1,2,3,4,5,6,7); */
        treg = 0;
        freg = 0;
    end

    // We latch the inputs so that we can refer to them in properties.
    // In particular, we could not refer to "valid" in properties without
    // these latches.
    always @ (posedge clock) begin
        freg = from;
        treg = to;
    end

    always @ (posedge clock) begin
        if (valid) begin
            b[treg] = b[freg];
            b[freg] = 0;
        end
    end

    // This predicate is true of the valid moves, that is, of moves that
    // swap the empty cell with one of its neighbors.
    assign valid = (b[treg] == 3'b000) &&
                   (((treg[1:0] == freg[1:0]) && !(treg[2] == freg[2])) ||
                    ((treg[2] == freg[2]) &&
                     (((treg[1:0]==2'b00) && (freg[1:0]==2'b01)) ||
                      ((treg[1:0]==2'b01) && (freg[0]==0)) ||
                      ((treg[1:0]==2'b10) && (freg[0]==1)) ||
                      ((treg[1:0]==2'b11) && (freg[1:0]==2'b10))
                      )
                     )
                    );

    // This is an invariant that is true of all board contents that are a
    // permutation of (0,...,7).  It's not really useful in this model,
    // because it is implied by "permutation," but it is not entirely
    // obvious, and hence makes for an interesting property to check.
    // The idea is to sum the Manhattan distances between each number and its
    // destination.  For the identity permutation this sum is 0, and every
    // transposition changes this sum by 2, 0, or -2.  Hence, for all
    // permutations, this sum is even.  To find the sum modulo 2, we just need
    // to know whether each distance is even or odd.  This in turn requires
    // only checking bits 2 and 0 of the cell index and its content.
    assign parity =
           (((b[0] & 5) == 1) | ((b[0] & 5) == 4)) ^
           (((b[1] & 5) == 0) | ((b[1] & 5) == 5)) ^
           (((b[2] & 5) == 1) | ((b[2] & 5) == 4)) ^
           (((b[3] & 5) == 0) | ((b[3] & 5) == 5)) ^
           (((b[4] & 5) == 1) | ((b[4] & 5) == 4)) ^
           (((b[5] & 5) == 0) | ((b[5] & 5) == 5)) ^
           (((b[6] & 5) == 1) | ((b[6] & 5) == 4)) ^
           (((b[7] & 5) == 0) | ((b[7] & 5) == 5));

    // This predicate is true of all board contents that represent a
    // permutation of (0,...,7).  As it is formulated, it seems to require
    // only that the sixteen numbers be different, but since these
    // numbers are on three bits, thy must be {0,...,7}.
    assign permutation = b[0]!=b[1] && b[0]!=b[2] && b[0]!=b[3] &&
           b[0]!=b[4] && b[0]!=b[5] && b[0]!=b[6] && b[0]!=b[7] &&
           b[1]!=b[2] && b[1]!=b[3] && b[1]!=b[4] && b[1]!=b[5] &&
           b[1]!=b[6] && b[1]!=b[7] && b[2]!=b[3] && b[2]!=b[4] &&
           b[2]!=b[5] && b[2]!=b[6] && b[2]!=b[7] && b[3]!=b[4] &&
           b[3]!=b[5] && b[3]!=b[6] && b[3]!=b[7] && b[4]!=b[5] &&
           b[4]!=b[6] && b[4]!=b[7] && b[5]!=b[6] && b[5]!=b[7] &&
           b[6]!=b[7];

    // This predicate is true of all board contents that:
    // 1. Have an odd number of inversions and the empty cell in
    //    an even-numbered row.
    // 2. Have an even number of inversions and the empty cell in
    //    an odd-numbered row.
    // The number of inversions is the number of pairs (i,j) such that:
    // 1. i > 0 and j > 0.
    // 2. i < j and b[i] > b[j].
    // The oddInversion is 0 for the identity permutation, and is invariant
    // under all legal moves.  In general, one can go from from one
    // permutation to another iff the two permutations have the same
    // value of oddInversions.
    assign oddInversions = (b[4]==0 || b[5]==0 || b[6]==0 || b[7]==0) ^
           (b[0]>b[1] && b[1]!=0) ^ (b[0]>b[2] && b[2]!=0) ^
           (b[0]>b[3] && b[3]!=0) ^ (b[0]>b[4] && b[4]!=0) ^
           (b[0]>b[5] && b[5]!=0) ^ (b[0]>b[6] && b[6]!=0) ^
           (b[0]>b[7] && b[7]!=0) ^ (b[1]>b[2] && b[2]!=0) ^
           (b[1]>b[3] && b[3]!=0) ^ (b[1]>b[4] && b[4]!=0) ^
           (b[1]>b[5] && b[5]!=0) ^ (b[1]>b[6] && b[6]!=0) ^
           (b[1]>b[7] && b[7]!=0) ^ (b[2]>b[3] && b[3]!=0) ^
           (b[2]>b[4] && b[4]!=0) ^ (b[2]>b[5] && b[5]!=0) ^
           (b[2]>b[6] && b[6]!=0) ^ (b[2]>b[7] && b[7]!=0) ^
           (b[3]>b[4] && b[4]!=0) ^ (b[3]>b[5] && b[5]!=0) ^
           (b[3]>b[6] && b[6]!=0) ^ (b[3]>b[7] && b[7]!=0) ^
           (b[4]>b[5] && b[5]!=0) ^ (b[4]>b[6] && b[6]!=0) ^
           (b[4]>b[7] && b[7]!=0) ^ (b[5]>b[6] && b[6]!=0) ^
           (b[5]>b[7] && b[7]!=0) ^ (b[6]>b[7] && b[7]!=0);


//#PASS:
assert property (permutation==1 |-> parity==0);
//#PASS:
assert property (permutation==1 |-> ##[0:100] permutation==1);
//#PASS:
assert property (permutation==1 && oddInversions==0 |-> ##[0:100] oddInversions==0);
//#PASS:
assert property (parity==0 |-> ##[0:100] parity==0);

endmodule // twoByFour
