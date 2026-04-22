module top(
  input clock,
  input reset,
  input io_interrupts_debug,
  output io_decode_0_fp_illegal
);
  wire csr_io_decode_0_fp_illegal;

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

  assign io_decode_0_fp_illegal = csr_io_decode_0_fp_illegal;
endmodule
