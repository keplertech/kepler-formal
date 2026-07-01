// Verilog translation of the original b06 circuit from the ITC99
// benchmark set.

// Author: Fabio Somenzi <Fabio@Colorado.EDU>


module main(CC_MUX, EQL, USCITE, clock, ENABLE_COUNT, ACKOUT, CONT_EQL);
    output [2:1] CC_MUX;
    input        EQL;
    output [2:1] USCITE;
    input        clock;
    output       ENABLE_COUNT;
    output       ACKOUT;
    input        CONT_EQL;

    parameter    cc_nop   = 2'b01;
    parameter    cc_enin  = 2'b01;
    parameter    cc_intr  = 2'b10;
    parameter    cc_ackin = 2'b11;
    parameter    out_norm = 2'b01;

    parameter s_init=0;
    parameter s_wait=1;
    parameter s_enin=2;
    parameter s_enin_w=3;
    parameter s_intr=4;
    parameter s_intr_1=5;
    parameter s_intr_w=6;

    reg [2:1]    CC_MUX, USCITE;
    reg          ENABLE_COUNT, ACKOUT;
    reg    state;

    initial begin
        state = s_init;
        CC_MUX = 0;
        ENABLE_COUNT = 0;
        ACKOUT = 0;
        USCITE = 0;
    end

    always @ (posedge clock) begin
        if (CONT_EQL) begin
            ACKOUT = 0;
            ENABLE_COUNT = 0;
        end else begin
            ACKOUT = 1;
            ENABLE_COUNT = 1;
        end

        case (state)

          s_init: begin
                CC_MUX = cc_enin;
                USCITE = out_norm;
                state = s_wait;
          end
          s_wait: begin
              if (EQL) begin
                  USCITE = 0;
                  CC_MUX = cc_ackin;
                  state = s_enin;
              end else begin
                  USCITE = out_norm;
                  CC_MUX = cc_intr;
                  state = s_intr_1;
              end
          end
          s_intr_1: begin
              if (EQL) begin
                  USCITE = 0;
                  CC_MUX = cc_ackin;
                  state = s_intr;
              end else begin
                  USCITE = out_norm;
                  CC_MUX = cc_enin;
                  state = s_wait;
              end
          end
          s_enin: begin
              if (EQL) begin
                  USCITE = 0;
                  CC_MUX = cc_ackin;
                  state = s_enin;
              end else begin
                  USCITE = 1;
                  ACKOUT = 1;
                  ENABLE_COUNT = 1;
                  CC_MUX = cc_enin;
                  state = s_enin_w;
              end
          end
          s_enin_w: begin
              if (EQL) begin
                  USCITE = 1;
                  CC_MUX = cc_enin;
                  state = s_enin_w;
              end else begin
                  USCITE = out_norm;
                  CC_MUX = cc_enin;
                  state = s_wait;
              end
          end
          s_intr: begin
              if (EQL) begin
                  USCITE = 0;
                  CC_MUX = cc_ackin;
                  state = s_intr;
              end else begin
                  USCITE = 3;
                  CC_MUX = cc_intr;
                  state = s_intr_w;
              end
          end
          s_intr_w: begin
              if (EQL) begin
                  USCITE = 3;
                  CC_MUX = cc_intr;
                  state = s_intr_w;
              end else begin
                  USCITE = out_norm;
                  CC_MUX = cc_enin;
                  state = s_wait;
              end
          end
        endcase
    end

  assert property (ENABLE_COUNT==ACKOUT);
  assert property (USCITE[2:1]!=2);
endmodule // b06
