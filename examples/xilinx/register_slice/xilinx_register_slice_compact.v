// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

module xilinx_register_slice(clk, rst, ce, d, x, q);
  input clk;
  wire clk;
  input rst;
  wire rst;
  input ce;
  wire ce;
  input d;
  wire d;
  input x;
  wire x;
  output q;
  wire q;
  wire lut_d;
  XOR2 explicit_xor(.A(d), .B(x), .Y(lut_d));
  FDRE #(
    .INIT(1'hx)
  ) _09_ (
    .C(clk),
    .CE(ce),
    .D(lut_d),
    .Q(q),
    .R(rst)
  );
endmodule
