// SEC wrapper for one Gemini trial's rtl/c4004_core.v (instantiated unchanged).
// The same wrapper is used for every trial, so both SEC sides have the top
//   cpu_top(clk, rst, test, d_in[3:0] -> bus_oe, bus_hiz, bus_d[3:0], sync,
//           cm_rom, cm_ram[3:0]).
// All ten trials have the same c4004_core ports (results/inventory.txt), and
// all reset to phase A1 with the PC taken from buf_l/m/h = 0, so no startup
// alignment is needed (unlike keplertech_demo/sec/wrap_ref.v).
// SRAM128 is sec/sram128_sec.v: no initial contents (X under dual_rail_steady).
module cpu_top (
    input  wire       clk,
    input  wire       rst,
    input  wire       test,
    input  wire [3:0] d_in,
    output wire       bus_oe,
    output wire       bus_hiz,
    output wire [3:0] bus_d,
    output wire       sync,
    output wire       cm_rom,
    output wire [3:0] cm_ram
);
    wire [3:0] d_core;
    wire       test_core;
    wire [2:0] env_ph;
    wire       env_ph_valid;
    wire [3:0] d_out;
    wire       d_oe, d_hiz, sync_raw;

    env_filter u_env (
        .clk(clk), .rst(rst), .d_in(d_in), .test(test), .sync(sync_raw),
        .d_core(d_core), .test_core(test_core), .ph(env_ph), .ph_valid(env_ph_valid)
    );

    c4004_core u_core (
        .clk(clk), .rst(rst), .test(test_core), .d_in(d_core),
        .d_out(d_out), .d_oe(d_oe), .d_hiz(d_hiz),
        .sync(sync_raw), .cm_rom(cm_rom), .cm_ram(cm_ram)
    );

    wire       oe_raw   = d_oe;
    wire       hiz_raw  = d_hiz;
    wire [3:0] dout_raw = d_out;
    assign sync = sync_raw;
// Observation layer (from keplertech_demo/sec/observe.vh, plus d_hiz).
// d_out is observed only while d_oe = 1: c4004_top drives the D pads with
// `D = d_oe ? d_out : 'z`, so d_out has no effect on the pins otherwise.
// `MASK_X1`: bus not observed during X1 (not used unless a run says so).



wire obs_mask = 1'b0;

assign bus_oe  = oe_raw  & ~obs_mask;
assign bus_hiz = hiz_raw | obs_mask;
assign bus_d   = bus_oe ? dout_raw : 4'd0;
endmodule
// Environment model shared by both SEC wrappers (instantiated identically in
// wrap_c4004.v and wrap_ref.v). It sits between the top-level d_in/test pins
// and the CPU core, and recovers the 4004 state phase from the design's OWN
// sync output (sync is high in X3, so the next clock is A1).
//
// Always on:
//   d_in and test are registered on the rising edge before reaching the core.
//   Kepler-Formal SEC evaluates negedge state with the same input vector as the
//   preceding posedge, while the C4004 routes d_in combinationally into its
//   negedge-clocked SRAM address (SRC pre-read in M2). Registering the inputs
//   removes that mixed-edge combinational input path; it is a pure 1-clock
//   delay applied identically to both designs.
//
// `STACK_ENV`: restrict programs to stack-legal ones. Tracks subroutine depth
//   (0..3). In M1/M2 of an instruction's FIRST word, a JMS at depth 3 or a BBL
//   at depth 0 is replaced by NOP (00). Two-word instructions: JCN(1), FIM(2,
//   OPA[0]=0), FIN(3, OPA[0]=0), JUN(4), JMS(5), ISZ(7).
//
// `TEST_CONST`: TEST is captured while rst is high and then held, so it is a
//   free but constant value per run.
module env_filter (
    input  wire       clk,
    input  wire       rst,
    input  wire [3:0] d_in,
    input  wire       test,
    input  wire       sync,      // this design's own SYNC (high in X3)
    output wire [3:0] d_core,
    output wire       test_core,
    output reg  [2:0] ph,        // recovered phase: 0=A1 ... 3=M1 4=M2 5=X1 ... 7=X3
    output reg        ph_valid   // 1 once the first SYNC has been seen
);
    // ---- input registers (always) ----
    reg [3:0] d_q;
    reg       t_q;
    always @(posedge clk) begin
        if (rst) begin
            d_q <= 4'd0;
            t_q <= test;
        end else begin
            d_q <= d_in;



            t_q <= test;

        end
    end





    assign test_core = t_q;


    // ---- phase recovery from the design's own sync ----
    always @(posedge clk) begin
        if (rst) begin
            ph       <= 3'd0;
            ph_valid <= 1'b0;
        end else if (sync) begin
            ph       <= 3'd0;
            ph_valid <= 1'b1;
        end else begin
            ph       <= ph + 3'd1;
        end
    end


















































    assign d_core = d_q;

endmodule
// SEC copy of rtl/sram128_model.v (SRAM128, 32 x 4-bit macro model).
// Identical behaviour, but WITHOUT the `initial` block that zeroes mem[] and DO:
// contents start unknown, as in silicon. The C4004 clears all 32 words itself
// during its post-reset init sequence, so nothing may rely on initial zeros.
//   - posedge CE with WR = 0: DO <= mem[addr]
//   - posedge CE with WR = 1: mem[addr] <= DIN (DO holds)
module SRAM128 (
    input  wire       CE,
    input  wire       WR,
    input  wire       A0,
    input  wire       A1,
    input  wire       A2,
    input  wire       A3,
    input  wire       A4,
    input  wire [3:0] DIN,
    output reg  [3:0] DO
);
    reg [3:0] mem [0:31];
    wire [4:0] addr = {A4, A3, A2, A1, A0};

    always @(posedge CE) begin
        if (WR) begin
            mem[addr] <= DIN;
        end else begin
            DO <= mem[addr];
        end
    end
endmodule
// Intel 4004 Compatible 4-bit CPU Core (c4004_core) - Trial 04 Optimized
// Target: OpenSUSI TR-1um (1um CMOS, 5.0V, 2-Metal M1/M2)
// Optimizations in Trial 04 (Fusion of Trial 02 & Trial 03 + Non-Stop Shadow PC):
//   1. Dedicated 4-bit Pipeline Register (reg_buf[3:0]):
//      - Latches early X1 SRAM128 index register pre-reads (ADD, SUB, LD, XCH, INC, ISZ, SRC)
//        WITHOUT ever overwriting buf_l (Shadow PC_L)!
//      - Latches 2-cycle M1/M2 operand nibbles for FIM/FIN so buf_l/buf_m/buf_h (Shadow PC)
//        are NEVER clobbered by any register, ALU, or FIM instruction!
//      - Thus pc_init_n stays 0 continuously across FIM, LD, ADD, SUB, INC, ISZ, XCH, SRC, JCN,
//        JUN, JMS, and JIN (only asserted on BBL stack return or FIN cycle 2 PC_H restore).
//   2. Full-Cycle FF-to-FF CLA ALU in X2:
//      - In X2, the 4-bit CLA operates purely FF-to-FF from reg_buf[3:0] and acc[3:0] over a
//        full clock period (1.0 * T_clk), achieving >= 15 MHz internal ALU timing!


module c4004_core (
    input  wire       clk,
    input  wire       rst,
    input  wire       test,
    input  wire [3:0] d_in,
    output reg  [3:0] d_out,
    output reg        d_oe,
    output wire       d_hiz,
    output wire       sync,
    output wire       cm_rom,
    output wire [3:0] cm_ram
);

    wire rst_n = ~rst;
    assign d_hiz = ~d_oe;

    // Core Architectural & Pipeline Flip-Flops (41 FFs total)
    reg [2:0] phase;      // 0:A1, 1:A2, 2:A3, 3:M1, 4:M2, 5:X1, 6:X2, 7:X3
    reg [3:0] acc;        // 4-bit Accumulator
    reg       cy;         // 1-bit Carry Flag
    reg [3:0] opr;        // 4-bit Instruction OPR
    reg [3:0] opa;        // 4-bit Instruction OPA
    reg [1:0] sp;         // 2-bit Stack Pointer (0..3)
    reg       two_cycle;  // 2nd cycle of 2-byte instruction
    reg       pc_init_n;  // 0: drive PC directly from buf_l/m/h FFs, 1: read PC from SRAM128
    reg       isz_taken;  // ISZ branch taken flag from cycle 1
    reg       pc_cout;    // Carry between A1->A2->A3 serial 4-bit PC incrementer
    reg [3:0] buf_l;      // Dedicated Shadow PC_L (or FIN address low nibble)
    reg [3:0] buf_m;      // Dedicated Shadow PC_M (or FIN address mid nibble)
    reg [3:0] buf_h;      // Dedicated Shadow PC_H
    reg [3:0] reg_buf;    // Dedicated 4-bit Pipeline Reg: X1 SRAM pre-read latch & 2nd-byte operand buffer!
    reg [2:0] dcl_reg;    // 3-bit Command Control Register (DCL)
    reg       wpm_half;   // 1-bit WPM half-byte toggle

    // Instruction decode helpers
    wire fin_cyc       = two_cycle & (opr == 4'h3);
    wire is_fim_fin_c2 = two_cycle & (opr[3:1] == 3'b001);
    wire is_src        = (~two_cycle) & (opr == 4'h2) & opa[0];
    wire is_jin        = (~two_cycle) & (opr == 4'h3) & opa[0];
    wire is_fin_c1     = (~two_cycle) & (opr == 4'h3) & (~opa[0]);
    wire is_io_ram     = (~two_cycle) & (opr == 4'hE);
    wire is_io_wr      = is_io_ram & (~opa[3]);
    wire is_reg_op     = (~two_cycle) & (opr == 4'h6 || opr == 4'h7 || (opr >= 4'h8 && opr <= 4'hB));

    // JCN branch condition evaluation
    wire jcn_raw   = (opa[2] & (acc == 4'd0)) | (opa[1] & cy) | (opa[0] & (~test));
    wire jcn_taken = opa[3] ^ jcn_raw;

    // -------------------------------------------------------------------------
    // SRAM128 Interface (A[4]=0: R0..R15, A[4]=1: Stack 4x3 nibbles)
    //   - Phase M1 (3): Write PC_L + 1 to stack
    //   - Phase M2 (4): Write R[2P] early for FIM/FIN cycle 2, else Write PC_M + c to stack
    //   - Phase X1 (5): Write R[2P+1] early for FIM/FIN cycle 2, or Pre-Read R[opa]/R0/R[2P] into reg_buf!
    //   - Phase X2 (6): Write PC_M + c (for FIM c2) or Write PC_H + c to stack
    //   - Phase X3 (7): Write PC_H + c (for FIM c2), Writeback R[opa] (from reg_buf), or Read R1/R[2P+1]
    // -------------------------------------------------------------------------
    reg  [4:0] sram_addr;
    reg  [3:0] sram_din;
    reg        sram_we;
    reg        sram_ce_en;
    wire [3:0] sram_do;

    always @(*) begin
        sram_addr  = {1'b1, sp, 2'b00};
        sram_din   = buf_l;
        sram_we    = 1'b0;
        sram_ce_en = 1'b1;

        case (phase)
            3'd0: begin // A1: Read PC_L (only needed when pc_init_n == 1)
                sram_addr  = {1'b1, sp, 2'b00};
                sram_we    = 1'b0;
                sram_ce_en = pc_init_n & (~fin_cyc);
            end
            3'd1: begin // A2: Read PC_M (only needed when pc_init_n == 1)
                sram_addr  = {1'b1, sp, 2'b01};
                sram_we    = 1'b0;
                sram_ce_en = pc_init_n & (~fin_cyc);
            end
            3'd2: begin // A3: Read PC_H (only needed when pc_init_n == 1)
                sram_addr  = {1'b1, sp, 2'b10};
                sram_we    = 1'b0;
                sram_ce_en = pc_init_n;
            end
            3'd3: begin // M1: Write PC_L + 1
                sram_addr  = {1'b1, sp, 2'b00};
                sram_din   = buf_l;
                sram_we    = ~fin_cyc;
                sram_ce_en = ~fin_cyc;
            end
            3'd4: begin // M2: Write R[2P] early for FIM/FIN cycle 2 (keeping buf_l/m/h untouched!), else Write PC_M + c
                if (is_fim_fin_c2) begin
                    sram_addr  = {1'b0, opa[3:1], 1'b0}; // Write R[2P] from reg_buf!
                    sram_din   = reg_buf;
                    sram_we    = 1'b1;
                    sram_ce_en = 1'b1;
                end else begin
                    sram_addr  = {1'b1, sp, 2'b01};
                    sram_din   = buf_m;
                    sram_we    = ~fin_cyc;
                    sram_ce_en = ~fin_cyc;
                end
            end
            3'd5: begin // X1: Write R[2P+1] early for FIM/FIN cycle 2, or Pre-Read Index Reg into reg_buf!
                if (is_fim_fin_c2) begin
                    sram_addr  = {1'b0, opa[3:1], 1'b1}; // Write R[2P+1] from reg_buf!
                    sram_din   = reg_buf;
                    sram_we    = 1'b1;
                    sram_ce_en = 1'b1;
                end else if (is_fin_c1) begin
                    sram_addr  = 5'h00; // Pre-read R0 for FIN
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b1;
                end else if (is_src || is_jin) begin
                    sram_addr  = {1'b0, opa[3:1], 1'b0}; // Pre-read R[2P] for SRC / JIN
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b1;
                end else if (is_reg_op) begin
                    sram_addr  = {1'b0, opa}; // Pre-read R[opa] into reg_buf!
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b1;
                end else begin
                    sram_addr  = {1'b1, sp, 2'b10}; // Write PC_H + c for non-reg ops
                    sram_din   = buf_h;
                    sram_we    = ~fin_cyc;
                    sram_ce_en = ~fin_cyc;
                end
            end
            3'd6: begin // X2: Write PC_M + c (for FIM c2) or Write PC_H + c (for reg_op / SRC / FIN / JIN)
                if (two_cycle && (opr == 4'h2)) begin
                    sram_addr  = {1'b1, sp, 2'b01}; // Write PC_M + c to stack for FIM c2
                    sram_din   = buf_m;
                    sram_we    = 1'b1;
                    sram_ce_en = 1'b1;
                end else if (is_reg_op || is_src || is_jin || is_fin_c1) begin
                    sram_addr  = {1'b1, sp, 2'b10}; // Write PC_H + c to stack in X2
                    sram_din   = buf_h;
                    sram_we    = 1'b1;
                    sram_ce_en = 1'b1;
                end else begin
                    sram_addr  = {1'b1, sp, 2'b10};
                    sram_din   = buf_h;
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b0;
                end
            end
            3'd7: begin // X3: Write PC_H + c (for FIM c2), Register Writeback (from reg_buf), or 2nd Nibble Read
                if (two_cycle && (opr == 4'h2)) begin
                    sram_addr  = {1'b1, sp, 2'b10}; // Write PC_H + c to stack for FIM c2
                    sram_din   = buf_h;
                    sram_we    = 1'b1;
                    sram_ce_en = 1'b1;
                end else if (is_fin_c1) begin
                    sram_addr  = 5'h01; // Read R1
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b1;
                end else if (is_src || is_jin) begin
                    sram_addr  = {1'b0, opa[3:1], 1'b1}; // Read R[2P+1]
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b1;
                end else if ((~two_cycle) && (opr == 4'h6 || opr == 4'h7 || opr == 4'hB)) begin
                    sram_addr  = {1'b0, opa}; // Write R[opa] from reg_buf for INC, ISZ, XCH!
                    sram_din   = reg_buf;
                    sram_we    = 1'b1;
                    sram_ce_en = 1'b1;
                end else begin
                    sram_addr  = {1'b0, opa};
                    sram_we    = 1'b0;
                    sram_ce_en = 1'b0;
                end
            end
        endcase
    end

    wire sram_ce = (~clk) & rst_n & sram_ce_en;

    SRAM128 u_sram (
        .A0  (sram_addr[0]),
        .A1  (sram_addr[1]),
        .A2  (sram_addr[2]),
        .A3  (sram_addr[3]),
        .A4  (sram_addr[4]),
        .CE  (sram_ce),
        .WR  (sram_we),
        .DIN (sram_din),
        .DO  (sram_do)
    );

    // -------------------------------------------------------------------------
    // Unified 4-bit Carry-Lookahead Adder:
    //   - Phase A1..A3 (0..2): Serial PC+1 incrementer
    //   - Phase X2     (6):    Full-cycle FF-to-FF ALU using reg_buf (latched from SRAM in X1!)
    // -------------------------------------------------------------------------
    wire [3:0] cur_pc_l = (fin_cyc | (~pc_init_n)) ? buf_l : sram_do;
    wire [3:0] cur_pc_m = (fin_cyc | (~pc_init_n)) ? buf_m : sram_do;
    wire [3:0] cur_pc_h = (~pc_init_n)             ? buf_h : sram_do;

    wire       daa_adj  = (acc > 4'd9) | cy;
    reg  [3:0] alu_a;
    reg  [3:0] alu_b;
    reg        alu_cin;

    always @(*) begin
        if (!phase[2]) begin
            // A1, A2, A3: PC serial incrementer
            alu_a   = (phase[1:0] == 2'd0) ? cur_pc_l :
                      (phase[1:0] == 2'd1) ? cur_pc_m : cur_pc_h;
            alu_b   = 4'd0;
            alu_cin = (phase[1:0] == 2'd0) ? 1'b1 : pc_cout;
        end else if (opr[3:1] == 3'b011) begin
            // INC (4'h6) or ISZ (4'h7): reg_buf + 1 (pure FF-to-FF!)
            alu_a   = reg_buf;
            alu_b   = 4'd0;
            alu_cin = 1'b1;
        end else if (opr == 4'h8) begin
            // ADD: acc + reg_buf + cy (pure FF-to-FF!)
            alu_a   = acc;
            alu_b   = reg_buf;
            alu_cin = cy;
        end else if (opr == 4'h9) begin
            // SUB: acc + ~reg_buf + ~cy (pure FF-to-FF!)
            alu_a   = acc;
            alu_b   = ~reg_buf;
            alu_cin = ~cy;
        end else if (opr == 4'hE) begin
            // ADM (opa=B) or SBM (opa=8)
            alu_a   = acc;
            alu_b   = opa[0] ? d_in : (~d_in);
            alu_cin = opa[0] ? cy   : (~cy);
        end else begin
            // OPR == 4'hF: IAC (opa=2), DAC (opa=8), DAA (opa=B)
            alu_a   = acc;
            alu_b   = (opa == 4'h8) ? 4'hF :
                      (opa == 4'hB && daa_adj) ? 4'd6 : 4'd0;
            alu_cin = (opa == 4'h2);
        end
    end

    // Explicit Carry-Lookahead Logic
    wire [3:0] cla_p = alu_a ^ alu_b;
    wire [3:0] cla_g = alu_a & alu_b;
    wire       c1    = cla_g[0] | (cla_p[0] & alu_cin);
    wire       c2    = cla_g[1] | (cla_p[1] & cla_g[0]) | (cla_p[1] & cla_p[0] & alu_cin);
    wire       c3    = cla_g[2] | (cla_p[2] & c2);
    wire       c4    = cla_g[3] | (cla_p[3] & c3);
    wire [4:0] alu_res = {c4, cla_p ^ {c3, c2, c1, alu_cin}};

    // -------------------------------------------------------------------------
    // Data Bus Output Multiplexer (d_out, d_oe) & External Control Signals
    //   Notice: In X2 (phase 6), SRC outputs reg_buf directly from FF at posedge clk!
    // -------------------------------------------------------------------------
    always @(*) begin
        d_oe  = 1'b0;
        d_out = 4'd0;
        case (phase)
            3'd0: begin
                d_oe  = 1'b1;
                d_out = cur_pc_l;
            end
            3'd1: begin
                d_oe  = 1'b1;
                d_out = cur_pc_m;
            end
            3'd2: begin
                d_oe  = 1'b1;
                d_out = cur_pc_h;
            end
            3'd6: begin
                if (is_src) begin
                    d_oe  = 1'b1;
                    d_out = reg_buf; // Driven directly from reg_buf FF at posedge clk!
                end else if (is_io_wr) begin
                    d_oe  = 1'b1;
                    d_out = acc;
                end
            end
            3'd7: begin
                if (is_src) begin
                    d_oe  = 1'b1;
                    d_out = sram_do;
                end
            end
            default: begin
                d_oe  = 1'b0;
                d_out = 4'd0;
            end
        endcase
    end

    assign sync   = (phase == 3'd7);
    assign cm_rom = ~((phase == 3'd2) || ((phase == 3'd6) && is_io_ram && (opa == 4'h2 || opa == 4'h3 || opa == 4'hA)));

    // DCL bank decode (active-low during A3 and during X2 for SRC/IO_RAM/DCL)
    wire cm_ram_active = (phase == 3'd2) ||
                         ((phase == 3'd6) && (~two_cycle) && (is_src || is_io_ram || (opr == 4'hF && opa == 4'hD)));
    wire [3:0] bank_sel = {dcl_reg[2], dcl_reg[1], dcl_reg[0], (dcl_reg == 3'b000)};
    assign cm_ram = cm_ram_active ? (~bank_sel) : 4'b1111;

    // -------------------------------------------------------------------------
    // Sequential State & Pipelined Execution Logic
    // -------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            phase     <= 3'd0;
            acc       <= 4'd0;
            cy        <= 1'b0;
            opr       <= 4'd0;
            opa       <= 4'd0;
            sp        <= 2'd0;
            two_cycle <= 1'b0;
            pc_init_n <= 1'b0;
            isz_taken <= 1'b0;
            pc_cout   <= 1'b0;
            buf_l     <= 4'd0;
            buf_m     <= 4'd0;
            buf_h     <= 4'd0;
            reg_buf   <= 4'd0;
            dcl_reg   <= 3'd0;
            wpm_half  <= 1'b0;
        end else begin
            phase <= phase + 3'd1;

            case (phase)
                3'd0: begin // A1
                    if (!fin_cyc) begin
                        buf_l   <= alu_res[3:0];
                        pc_cout <= alu_res[4];
                    end
                end
                3'd1: begin // A2
                    if (!fin_cyc) begin
                        buf_m   <= alu_res[3:0];
                        pc_cout <= alu_res[4];
                    end
                end
                3'd2: begin // A3
                    if (!fin_cyc) begin
                        buf_h     <= alu_res[3:0];
                        pc_init_n <= 1'b0; // buf_l/m/h always hold valid PC+1!
                    end
                end
                3'd3: begin // M1
                    if (!two_cycle) begin
                        opr <= d_in;
                    end else begin
                        reg_buf <= d_in; // Latch 2nd-byte high nibble into reg_buf (buf_l/m/h stay untouched!)
                        if (!is_fim_fin_c2) begin
                            pc_init_n <= 1'b1; // Only JCN/ISZ/JUN/JMS might clobber buf_m in M2
                        end
                    end
                end
                3'd4: begin // M2
                    if (!two_cycle) begin
                        opa <= d_in;
                    end else if (is_fim_fin_c2) begin
                        reg_buf <= d_in; // Latch 2nd-byte low nibble into reg_buf (R[2P] written in M2!)
                    end else begin
                        buf_m <= d_in;   // 2nd-byte low nibble for JCN/ISZ/JUN/JMS
                    end
                end
                3'd5: begin // X1: Latch Pre-Read SRAM128 Index Reg into reg_buf (buf_l stays untouched!)
                    if (is_fin_c1 || is_jin) begin
                        buf_m <= sram_do; // R0 (for FIN) or R[2P] (for JIN)
                    end else if (is_reg_op || is_src) begin
                        reg_buf <= sram_do; // R[opa] (for ALU/LD/XCH/INC/ISZ) or R[2P] (for SRC) -> buf_l untouched!
                    end
                end
                3'd6: begin // X2: Main Execution Phase (Full Clock Period from reg_buf FF!)
                    if (!two_cycle) begin
                        case (opr)
                            4'h6: begin // INC Rr
                                reg_buf <= alu_res[3:0]; // buf_l stays untouched!
                            end
                            4'h7: begin // ISZ Rr (1st cycle)
                                reg_buf   <= alu_res[3:0]; // buf_l stays untouched!
                                isz_taken <= (reg_buf != 4'hF);
                            end
                            4'h8, 4'h9: begin // ADD, SUB (FF-to-FF from reg_buf!)
                                acc <= alu_res[3:0];
                                cy  <= alu_res[4];
                            end
                            4'hA: begin // LD Rr (FF-to-FF from reg_buf!)
                                acc <= reg_buf;
                            end
                            4'hB: begin // XCH Rr (FF-to-FF swap between acc and reg_buf!)
                                acc     <= reg_buf;
                                reg_buf <= acc;
                            end
                            4'hC: begin // BBL d
                                sp        <= sp - 2'd1;
                                acc       <= opa;
                                pc_init_n <= 1'b1; // Read return PC from SRAM128 at new SP
                            end
                            4'hD: begin // LDM d
                                acc <= opa;
                            end
                            4'hE: begin // I/O and RAM instructions
                                case (opa)
                                    4'h3: wpm_half <= ~wpm_half;
                                    4'h8, 4'hB: begin // SBM, ADM
                                        acc <= alu_res[3:0];
                                        cy  <= alu_res[4];
                                    end
                                    4'h9, 4'hA, 4'hC, 4'hD, 4'hE, 4'hF: begin // RDM, RDR, RD0..RD3
                                        acc <= d_in;
                                    end
                                    default: ;
                                endcase
                            end
                            4'hF: begin // Accumulator Group
                                case (opa)
                                    4'h0: begin // CLB
                                        acc <= 4'd0;
                                        cy  <= 1'b0;
                                    end
                                    4'h1: cy  <= 1'b0; // CLC
                                    4'h2, 4'h8: begin // IAC, DAC
                                        acc <= alu_res[3:0];
                                        cy  <= alu_res[4];
                                    end
                                    4'h3: cy  <= ~cy;  // CMC
                                    4'h4: acc <= ~acc; // CMA
                                    4'h5: begin // RAL
                                        acc <= {acc[2:0], cy};
                                        cy  <= acc[3];
                                    end
                                    4'h6: begin // RAR
                                        acc <= {cy, acc[3:1]};
                                        cy  <= acc[0];
                                    end
                                    4'h7: begin // TCC
                                        acc <= {3'd0, cy};
                                        cy  <= 1'b0;
                                    end
                                    4'h9: begin // TCS
                                        acc <= cy ? 4'd10 : 4'd9;
                                        cy  <= 1'b0;
                                    end
                                    4'hA: cy  <= 1'b1; // STC
                                    4'hB: begin // DAA
                                        acc <= alu_res[3:0];
                                        if (alu_res[4]) cy <= 1'b1;
                                    end
                                    4'hC: begin // KBP
                                        case (acc)
                                            4'b0000: acc <= 4'd0;
                                            4'b0001: acc <= 4'd1;
                                            4'b0010: acc <= 4'd2;
                                            4'b0100: acc <= 4'd3;
                                            4'b1000: acc <= 4'd4;
                                            default: acc <= 4'hF;
                                        endcase
                                    end
                                    4'hD: dcl_reg <= acc[2:0]; // DCL
                                    default: ;
                                endcase
                            end
                            default: ;
                        endcase
                    end
                end
                3'd7: begin // X3: 2nd Nibble Read/Write & Two-Cycle Branch Commit
                    if (!two_cycle) begin
                        if ((opr == 4'h1) ||
                            (opr == 4'h2 && !opa[0]) ||
                            (opr == 4'h3 && !opa[0]) ||
                            (opr == 4'h4) ||
                            (opr == 4'h5) ||
                            (opr == 4'h7)) begin
                            two_cycle <= 1'b1;
                        end
                        if (opr == 4'h3) begin // FIN (read R1) or JIN (read R[2P+1])
                            buf_l <= sram_do;
                            if (opa[0]) begin // JIN jumps to {PC_H, R[2P], R[2P+1]}
                                pc_init_n <= 1'b0;
                            end else begin
                                pc_init_n <= 1'b1; // FIN cycle 2 A3 reads PC_H from SRAM128
                            end
                        end
                    end else begin
                        two_cycle <= 1'b0;
                        case (opr)
                            4'h1: begin // JCN cycle 2
                                if (jcn_taken) begin
                                    buf_l     <= buf_m;
                                    buf_m     <= reg_buf;
                                    pc_init_n <= 1'b0;
                                end
                            end
                            4'h4: begin // JUN cycle 2
                                buf_l     <= buf_m;
                                buf_m     <= reg_buf;
                                buf_h     <= opa;
                                pc_init_n <= 1'b0;
                            end
                            4'h5: begin // JMS cycle 2
                                sp        <= sp + 2'd1;
                                buf_l     <= buf_m;
                                buf_m     <= reg_buf;
                                buf_h     <= opa;
                                pc_init_n <= 1'b0;
                            end
                            4'h7: begin // ISZ cycle 2
                                if (isz_taken) begin
                                    buf_l     <= buf_m;
                                    buf_m     <= reg_buf;
                                    pc_init_n <= 1'b0;
                                end
                            end
                            default: ;
                        endcase
                    end
                end
            endcase
        end
    end

endmodule
