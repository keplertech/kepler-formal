module top(
  input clock,
  input reset,
  input io_interrupts_debug,
  output io_decode_0_fp_illegal
);
  wire csr_io_decode_0_fp_illegal;
  wire reset_n;
  wire s1d;
  wire s1q;
  wire s1qn;
  wire s2d;
  wire s2q;
  wire s2qn;
  wire s3d;
  wire s3q;
  wire s3qn;
  wire edited_io_decode_0_fp_illegal;

  CSRFile dut(
    .clock(clock),
    .reset(reset),
    .io_ungated_clock(1'b0),
    .io_interrupts_debug(io_interrupts_debug),
    .io_interrupts_mtip(1'b0),
    .io_interrupts_msip(1'b0),
    .io_interrupts_meip(1'b0),
    .io_rw_addr(12'b0),
    .io_rw_cmd(3'b0),
    .io_rw_wdata(32'b0),
    .io_decode_0_csr(12'b0),
    .io_exception(1'b0),
    .io_retire(1'b0),
    .io_cause(32'b0),
    .io_pc(32'b0),
    .io_tval(32'b0),
    .io_inst_0(32'b0),
    .io_decode_0_fp_illegal(csr_io_decode_0_fp_illegal)
  );

  INV_X1 inv0(.A(reset), .ZN(reset_n));
  AND2_X1 g1(.A1(io_interrupts_debug), .A2(reset_n), .ZN(s1d));
  DFF_X1 ff1(.D(s1d), .CK(clock), .Q(s1q), .QN(s1qn));
  AND2_X1 g2(.A1(s1q), .A2(reset_n), .ZN(s2d));
  DFF_X1 ff2(.D(s2d), .CK(clock), .Q(s2q), .QN(s2qn));
  AND2_X1 g3(.A1(s2q), .A2(reset_n), .ZN(s3d));
  DFF_X1 ff3(.D(s3d), .CK(clock), .Q(s3q), .QN(s3qn));
  XOR2_X1 x1(
    .A(csr_io_decode_0_fp_illegal),
    .B(s3q),
    .Z(edited_io_decode_0_fp_illegal)
  );

  assign io_decode_0_fp_illegal = edited_io_decode_0_fp_illegal;
endmodule
