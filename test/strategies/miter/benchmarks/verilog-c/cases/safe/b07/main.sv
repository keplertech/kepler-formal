// Verilog translation of the original b07 circuit from the ITC99
// benchmark set.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>


module main(PUNTI_RETTA, START, clock);
    output [7:0] PUNTI_RETTA;
    input        START;
    input        clock;

    reg [7:0]    PUNTI_RETTA;
    reg    stato;
    reg [7:0]    cont, x;
    reg [6:0]    y, t;
    reg [3:0]    mar;
    wire [7:0]   mem_mar;
    parameter S_RESET = 0;
    parameter S_START = 1;
    parameter S_LOAD_X = 2;
    parameter S_UPDATE_MAR = 3;
    parameter S_LOAD_Y = 4;
    parameter S_CALC_RETTA = 5;
    parameter S_INCREMENTA = 6;

    function [7:0] mem;
        input [3:0] a;
        begin: _mem
            if      (a ==  0) mem = 1;
            else if (a ==  1) mem = 255;
            else if (a ==  2) mem = 0;
            else if (a ==  3) mem = 0;
            else if (a ==  4) mem = 0;
            else if (a ==  5) mem = 2;
            else if (a ==  6) mem = 0;
            else if (a ==  7) mem = 0;
            else if (a ==  8) mem = 0;
            else if (a ==  9) mem = 2;
            else if (a == 10) mem = 255;
            else if (a == 11) mem = 5;
            else if (a == 12) mem = 0;
            else if (a == 13) mem = 2;
            else if (a == 14) mem = 0;
            else              mem = 2;
        end
    endfunction // mem

    assign mem_mar = mem(mar);

    initial begin
        stato = S_RESET;
        PUNTI_RETTA = 0;
        cont = 0;
        mar = 0;
        x = 0;
        y = 0;
        t = 0;
    end

    always @ (posedge clock) begin
        case (stato)
          S_RESET: begin
              stato = S_START;
          end
          S_START: begin
              if (START) begin
                  cont = 0;
                  mar = 0;
                  stato = S_LOAD_X;
              end else begin
                  stato = S_START;
                  PUNTI_RETTA = 0;
              end
          end
          S_LOAD_X: begin
              x = mem_mar;
              stato = S_UPDATE_MAR;
          end
          S_UPDATE_MAR: begin
              mar = mar + 1;
              t = {x[5:0], 1'b0};
              stato = S_LOAD_Y;
          end
          S_LOAD_Y: begin
              y = mem_mar[6:0];
              x = {1'b0, x[6:0]} + {1'b0, t};
              stato = S_CALC_RETTA;
          end
          S_CALC_RETTA: begin
              x = {1'b0, x[6:0]} + {1'b0, y};
              stato = S_INCREMENTA;
          end
          S_INCREMENTA: begin
              if (mar != 15) begin
                  if (x == 2)
                    cont = cont + 1;
                  mar = mar + 1;
                  stato = S_LOAD_X;
              end else begin
                  if (START == 0) begin
                      if (x == 2)
                        PUNTI_RETTA = cont + 1;
                      else
                        PUNTI_RETTA = cont;
                      stato = S_START;
                  end else begin
                      stato = S_INCREMENTA;
                  end
              end
          end
        endcase
    end
 assert property (x[7:0]!=148);

endmodule // b07
