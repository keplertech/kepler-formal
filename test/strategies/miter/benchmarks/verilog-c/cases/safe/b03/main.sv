// Verilog translation of the original b03 circuit from the ITC99
// benchmark set.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main(clock, REQUEST1, REQUEST2, REQUEST3, REQUEST4, GRANT_O);
    input        clock;
    input        REQUEST1, REQUEST2, REQUEST3, REQUEST4;
    output [3:0] GRANT_O;

    parameter    U1 = 3'b100;
    parameter    U2 = 3'b010;
    parameter    U3 = 3'b001;
    parameter    U4 = 3'b111;
    parameter INIT = 0;
    parameter ANALISI_REQ = 1;
    parameter ASSEGNA = 2;


    reg [3:0]    GRANT_O;
    reg [1:0]    stato;
    reg [2:0]    coda0, coda1, coda2, coda3;
    reg          ru1, ru2, ru3, ru4, fu1, fu2, fu3, fu4;
    reg [3:0]    grant;

    initial begin
        stato = INIT;
        coda0 = 0;
        coda1 = 0;
        coda2 = 0;
        coda3 = 0;
        ru1 = 0;
        fu1 = 0;
        ru2 = 0;
        fu2 = 0;
        ru3 = 0;
        fu3 = 0;
        ru4 = 0;
        fu4 = 0;
        grant = 0;
        GRANT_O = 0;
    end

    always @ (posedge clock) begin
        case (stato)
          INIT: begin
              ru1 = REQUEST1;
              ru2 = REQUEST2;
              ru3 = REQUEST3;
              ru4 = REQUEST4;
              stato = ANALISI_REQ;
          end
          ANALISI_REQ: begin
              GRANT_O = grant;
              if (ru1) begin
                  if (!fu1) begin
                      coda3 = coda2;
                      coda2 = coda1;
                      coda1 = coda0;
                      coda0 = U1;
                  end
              end else if (ru2) begin
                  if (!fu2) begin
                      coda3 = coda2;
                      coda2 = coda1;
                      coda1 = coda0;
                      coda0 = U2;
                  end
              end else if (ru3) begin
                  if (!fu3) begin
                      coda3 = coda2;
                      coda2 = coda1;
                      coda1 = coda0;
                      coda0 = U3;
                  end
              end else if (ru4) begin
                  if (!fu4) begin
                      coda3 = coda2;
                      coda2 = coda1;
                      coda1 = coda0;
                      coda0 = U4;
                  end
              end

              fu1 = ru1;
              fu2 = ru2;
              fu3 = ru3;
              fu4 = ru4;

              stato = ASSEGNA;
          end
          ASSEGNA: begin
              if (fu1 | fu2 | fu3 | fu4) begin
                  case (coda0)
                    U1: grant = 4'b1000;
                    U2: grant = 4'b0100;
                    U3: grant = 4'b0010;
                    U4: grant = 4'b0001;
                    default: grant = 4'b0000;
                  endcase
                  coda0 = coda1;
                  coda1 = coda2;
                  coda2 = coda3;
                  coda3 = 0;
              end
              ru1 = REQUEST1;
              ru2 = REQUEST2;
              ru3 = REQUEST3;
              ru4 = REQUEST4;
              stato = ANALISI_REQ;
          end
        endcase
    end
  assert property (GRANT_O==0 || GRANT_O==8 || GRANT_O==4 || GRANT_O==2 || GRANT_O==1);
endmodule // b03
