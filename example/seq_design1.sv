module top(
    input logic clk,
    input logic rst,
    input logic d,
    output logic q
);
  always_ff @(posedge clk)
  if (rst) begin
    q <= 1'b0;
  end else begin
    q <= d;
  end
endmodule
