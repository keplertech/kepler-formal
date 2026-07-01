//-----------------------------------------------------------------------------
// Copyright 1997 VAutomation Inc. Nashua NH USA (603) 882-2282.
// This software is provided AS-IS and free of charge with the restriction that
// this copyright notice remain in all copies of the Source Code at all times.
// Visit HTTP://www.vautomation.com for more information on our cores.
//
// File: vcrc32_8.v     Ethernet/FDDI/1394 8-bit parallel 32-bit CRC generator.
// Description:
//      This module does an 8-bit parallel CRC generation.  The polynomial
//      is that specified for IEEE 802.3 (ethernet) LANs and other standards.
//
//
// I. Functionality:
//
// 1. The module does an 8-bit parallel CRC generation.
// 2. The module provides a synchronous 8-bit per clock load and unload
//    function.
// 3. The polynomial is that specified for 802.3 LANs and other standards.
//
// The polynomial computed is:
//
// G(x)=X**32+X**26+X**23+X**22+X**16+X**12+X**11+X**10+X**8+X**7+X**5+X**4+X
// **2+X+1
//
// II. Module I/O
//
// Inputs: clk, clken, reset, load, compute, data_in[7:0]
//
// outputs: crc_ok, data_out[7:0], crc[31:0]
//
//
// III.Truth Table:
//
// clken  reset  compute   load    | data_out
// ------------------------------------------
//  0      X      X        X       | No change
//  1      0      0        0       | No change
//  1      1      X        X       | 0xFFFF (all ones)
//  1      0      X        1       | load and shift 1 byte of CRC
//  1      0      1        0       | Compute CRC
//
//  0      0      1        1       | unload 4 byte CRC (NOT IMPLEMENTED)
//
// Loading and unloading of the 32-bit CRC register is done one byte at a time
// by asserting load and clken.  The data on data_in is shifted into the the
// LSB of the CRC register. The MSB of the CRC register is available on
// data_out.
//
//-----------------------------------------------------------------------------
module main(clk,clken,reset,load,compute,data_in,data_out,crc_ok,crc);
    input         clk;          // everything clocks on rising edge
    input         clken;        // Clock enable
    input         reset;        // RESET
    input         load;         // Load the CRC register
    input         compute;      // compute the next 8-bits of CRC
    input [7:0]   data_in;      // data in
    output [7:0]  data_out;
    output        crc_ok;       // the last CRC check cycle was OK
    output [31:0] crc;          // direct acces to the CRC value

    parameter CRC_INITIAL_VALUE = 32'b11111111111111111111111111111111;
    // 0xC704DD7B
    parameter CRC_REMAINDER = 32'b11000111000001001101110101111011;

    reg [31:0]    crc;

    wire [31:0]   new_crc;      // D inputs to the flops

  function [31:0] parallel_crc;
        // Because the equations are long I am keeping the name of
        // the incoming CRC and the XOR vector short.
        input [31:0] c;         // crc_in
        input [7:0] d;          // data_in
    begin: _parallel_crc
        reg [31:24] x;

        // The first thing that a parallel CRC needs to do is to develop
        // the vector formed by XORing the input vector by the current CRC.
        // This vector is then used during the CRC calculation.

        x = c[31:24] ^ d;

        // The parellel CRC is a function of the X vector and the current CRC.
        parallel_crc[31] = x[29] ^ c[23];
        parallel_crc[30] = x[28] ^ x[31] ^ c[22];
        parallel_crc[29] = x[27] ^ x[30] ^ x[31] ^ c[21];
        parallel_crc[28] = x[26] ^ x[29] ^ x[30] ^ c[20];
        parallel_crc[27] = x[31] ^ x[25] ^ x[28] ^ x[29] ^ c[19];
        parallel_crc[26] = x[30] ^ x[24] ^ x[27] ^ x[28] ^ c[18];
        parallel_crc[25] = x[26] ^ x[27] ^ c[17];
        parallel_crc[24] = x[31] ^ x[25] ^ x[26] ^ c[16];
        parallel_crc[23] = x[30] ^ x[24] ^ x[25] ^ c[15];
        parallel_crc[22] = x[24] ^ c[14];
        parallel_crc[21] = x[29] ^ c[13];
        parallel_crc[20] = x[28] ^ c[12];
        parallel_crc[19] = x[27] ^ x[31] ^ c[11];
        parallel_crc[18] = x[26] ^ x[30] ^ x[31] ^ c[10];
        parallel_crc[17] = x[25] ^ x[29] ^ x[30] ^ c[9];
        parallel_crc[16] = x[24] ^ x[28] ^ x[29] ^ c[8];
        parallel_crc[15] = x[27] ^ x[28] ^ x[29] ^ x[31] ^ c[7];
        parallel_crc[14] = x[26] ^ x[27] ^ x[28] ^ x[30] ^ x[31] ^ c[6];
        parallel_crc[13] = x[31] ^ x[25] ^ x[26] ^ x[27] ^ x[29] ^ x[30]
                           ^ c[5];
        parallel_crc[12] = x[30] ^ x[24] ^ x[25] ^ x[26] ^ x[28] ^ x[29]
                           ^ c[4];
        parallel_crc[11] = x[24] ^ x[25] ^ x[27] ^ x[28] ^ c[3];
        parallel_crc[10] = x[24] ^ x[26] ^ x[27] ^ x[29] ^ c[2];
        parallel_crc[9] = x[25] ^ x[26] ^ x[28] ^ x[29] ^ c[1];
        parallel_crc[8] = x[24] ^ x[25] ^ x[27] ^ x[28] ^ c[0];
        parallel_crc[7] = x[24] ^ x[26] ^ x[27] ^ x[29] ^ x[31];
        parallel_crc[6] = x[25] ^ x[26] ^ x[28] ^ x[29] ^ x[30] ^ x[31];
        parallel_crc[5] = x[31] ^ x[30] ^ x[29] ^ x[28] ^ x[27] ^ x[25]
                          ^ x[24];
        parallel_crc[4] = x[30] ^ x[28] ^ x[27] ^ x[26] ^ x[24];
        parallel_crc[3] = x[31] ^ x[25] ^ x[26] ^ x[27];
        parallel_crc[2] = x[30] ^ x[24] ^ x[31] ^ x[25] ^ x[26];
        parallel_crc[1] = x[30] ^ x[24] ^ x[31] ^ x[25];
        parallel_crc[0] = x[30] ^ x[24];
    end // block: _parallel_crc
    endfunction // parallel_crc

    initial crc = CRC_INITIAL_VALUE;

    always @ (posedge clk) begin
    // Create the 32 Flip flops (with clock enable flops)
        if (clken)
          crc = new_crc;
    end

    // Concurrent signal assignment for the crc generation.
    assign new_crc = reset ? CRC_INITIAL_VALUE :
                     load ? {crc[23:0], data_in} :
                     compute ? parallel_crc(crc,data_in) :
                     crc;

    assign data_out = ~crc[31:24];

    assign crc_ok = crc == CRC_REMAINDER;
    // This is a 32-bit wide AND function, so proper attention should be paid
    // when synthesizing to achieve good results.  If there are cycles
    // available pipeling this gate would be appropriate.

   assert property (crc[31:0] != b10101010101010101010101010101010);

endmodule // vcrc32_8
