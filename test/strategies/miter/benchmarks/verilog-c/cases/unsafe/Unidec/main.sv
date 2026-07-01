// Nondeterministic checker for unique decipherability.
//
// This model illustrates how a search problem may be turned into a
// reachability problem for a nondeterministic machine.
//
// Here we want to know if a given code is uniquely decipherable.  For that,
// there must be no concatenation of code words that can be divided in code
// words in more than one way.
//
// Consider the following code taken from "Information theory" by Robert
// V. Ash, Dover,p. 29.  It consists of seven code words:
//
//  a, c, ad, abb, bad, deb, bbcde.
//
// The following message is not uniquely decipherable: abbcdebad.  Indeed, it
// can be parsed as either a,bbcde,bad or abb,c,deb,ad.  Let us review the
// ideas behind the search for such an ambiguous message.
//
// Since a and abb are both code words, if a message in our code starts
// with abb, we do not know whether the first word is a or abb.
// Suppose the first word is a, then the code is ambiguous if there is a
// message bbS, such that S is also a valid message.
//
//   a|bb
//   abb|
//
// A message starting with bb must indeed start with bbcde.  Hence, S must
// start with cde.
//
//   a|bbcde|
//   abb|cde
//
// A message that starts with cde must start with c.  Hence, for such a
// message to exist, there must exist a message starting with de.
//
//   a|bbcde|
//   abb|c|de
//
// Such a message must start with deb.  The situation so far is this:
//
//   a|bbcde|b
//   abb|c|deb|
//
// Therefore, the first interpretation must continue with b.  Suppose it
// continues with bad.  Then, the second interpretation can be completed
// with ad.
//
//   a|bbcde|bad|
//   abb|c|deb|ad|
//
// At each round there is an open suffix in one of the two interpretations,
// that we try to close by removing a prefix that is a code word
// (cde -> c|de in our example) or by appending more letters.  If the prefix
// we remove is the entire suffix (as in our example) we are done.
// Appending more letters, on the other hand, creates an open suffix in the
// other interpretation.
//
// If a suffix repeats, we have to backtrack.  For instance, in the example
// above, if we had chosen bbcde instead of bad, we would have found it
// impossible to continue, and we should have backtracked.
// This is just the intuition behind the algorithm.  More details can be
// found in the cited book.  The book does not talk about backtracking.
// Rather, it employs breadth-first search.
// The nondeterministic machine obviously does not backtrack, but rather
// searches all paths in parallel, and therefore also does BFS, though with
// minor differences from what done in the book.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main(clk, sel1, sel2, found);
    input       clk;
    input [2:0] sel1;
    input [1:0] sel2;
    output      found;

    wire [15:0] other;
    reg [15:0]  word;
    reg         found;
    reg         init;

    // This function returns a code word.  Each character is three bits
    // (a: 000, b: 001, c: 010, d: 011, and e: 100).  Each word is up to
    // five characters plus a "stop" bit.  Characters are stored in a word
    // in reverse order so that a right shift produces a suffix.
  function [15:0] code;
        input [2:0] sel;
        begin: _code
            case (sel)
              0: code = 16'b0000000000001000; // a
              1: code = 16'b0000000000001010; // c
              2: code = 16'b0000000001011000; // ad
              3: code = 16'b0000001001001000; // abb
              4: code = 16'b0000001011000001; // bad
              5: code = 16'b0000001001100011; // deb
              6: code = 16'b1100011010001001; // bbcde
              7: code = 16'b0000000000001000; // a
            endcase
        end
    endfunction // code

    // This function extracts a proper prefix of lengh sel+1 from word.
    // If word does not have more than sel+1 characters, it returns
    // an invalid word.
  function [15:0] prefix;
        input [15:0] word;
        input [1:0] sel;
        begin: _prefix
            case (sel)
              0: if (word[15:4] == 0)  prefix = 16'b0111111111111111;
              else prefix = {13'b1, word[2:0]};
              1: if (word[15:7] == 0)  prefix = 16'b0111111111111111;
              else prefix = {10'b1, word[5:0]};
              2: if (word[15:10] == 0) prefix = 16'b0111111111111111;
              else prefix = { 7'b1, word[8:0]};
              3: if (word[15:13] == 0) prefix = 16'b0111111111111111;
              else prefix = { 4'b1, word[11:0]};
            endcase
        end
    endfunction // prefix

    // This function returns a suffix of word dropping the first sel+1
    // characters.
  function [15:0] suffix;
        input [15:0] word;
        input [1:0] sel;
        begin: _suffix
            case (sel)
              0: suffix = { 3'b0, word[15:3]};
              1: suffix = { 6'b0, word[15:6]};
              2: suffix = { 9'b0, word[15:9]};
              3: suffix = {12'b0, word[15:12]};
            endcase
        end
    endfunction // suffix


 initial begin
        word = code(sel1);
        found = 0;
        init = 1;
    end

    assign other = code(sel1);

    always @ (posedge clk) begin
        found = !init && word == other;
        init = 0;
        if (other == prefix(word,sel2)) begin
            // There is a code word that is a prefix of the current word.
            // Make the suffix of the current word the next word.
            word = suffix(word,sel2);
        end else if (prefix(other,sel2) == word) begin
            // The current word is a prefix of another code word.
            // Make the suffix of the other word the next word.
            word = suffix(other,sel2);
        end else begin
            // Neither applies.  Go to trap state.
            word = 0;
        end
    end

  //#PASS: The trap state lives up to its name.
  //assert property ((word[15:0]==0) |-> ##1 (word[15:0]==0)); // FAILING

  //#PASS: Eventually the trap is inevitable.
  //assert property (true |-> ##[0:$] word[15:0]==0); // FAILING

  // FAIL
  assert property (found==0);
endmodule // unidec
