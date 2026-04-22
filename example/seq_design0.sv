module top(
    input  logic clk,
    input  logic rst,
    input  logic d,
    output logic q
);

  logic [3:0] cycle;   // enough bits to count to 11

  always_ff @(posedge clk) begin
    if (rst) begin
      cycle <= 4'd0;
      q     <= 1'b0;
    end else begin
      cycle <= cycle + 4'd1;

      if (cycle == 4'd11)
        q <= 1'b1;      // force q = 1 on cycle 11
      else
        q <= d;         // normal behavior
    end
  end

endmodule

