// Verilog translation of the original b08 circuit from the ITC99
// benchmark set.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main(CLOCK, START, I, O);
    input        CLOCK;
    input        START;
    input [7:0]  I;
    output [3:0] O;

    reg [3:0]    O;
    reg    STATO;

    parameter start_st = 0;
    parameter init = 1;
    parameter loop_st = 2;
    parameter the_end = 3;

    wire [7:0]   ROM_1, ROM_2;
    reg [7:0]    IN_R;
    reg [2:0]    MAR;
    reg [3:0]    OUT_R;
    wire [3:0]   ROM_OR;
    wire [19:0]  ROM_OUT;

    function [19:0] ROM;
        input [2:0] a;
        begin: _ROM
            if (a == 0)
              ROM = 20'b01111111100101111010;
            else if (a == 1)
              ROM = 20'b00111001110101100010;
            else if (a == 2)
              ROM = 20'b10101000111111111111;
            else if (a == 3)
              ROM = 20'b11111111011010111010;
            else if (a == 4)
              ROM = 20'b11111111111101101110;
            else if (a == 5)
              ROM = 20'b11111111101110101000;
            else if (a == 6)
              ROM = 20'b11001010011101011011;
            else
              ROM = 20'b00101111111111110100;
        end // block: _ROM
    endfunction // ROM

    initial begin
        STATO = start_st;
        MAR = 0;
        IN_R = 0;
        OUT_R = 0;
        O = 0;
    end

    assign ROM_OUT = ROM(MAR);
    assign ROM_1 = ROM_OUT[19:12];
    assign ROM_2 = ROM_OUT[11:4];
    assign ROM_OR = ROM_OUT[3:0];

    always @ (posedge CLOCK) begin
        case (STATO)
          start_st: begin
              if (START) STATO = init;
          end
          init: begin
              IN_R  = I;
              OUT_R = 0;
              MAR   = 0;
              STATO = loop_st;
          end
          loop_st: begin
              if (((ROM_2 & ~IN_R) | (ROM_1 & IN_R) | (ROM_2 & ROM_1)) ==
                  8'b11111111) begin
                  OUT_R = OUT_R | ROM_OR;
              end
              STATO = the_end;
          end
          the_end: begin
              if (MAR != 7) begin
                  MAR = MAR + 1;
                  STATO = loop_st;
              end else if (!START) begin
                  O = OUT_R;
                  STATO = start_st;
              end
          end
        endcase
    end
  assert property (ROM_OR[3:0]!=0);
endmodule // b08
