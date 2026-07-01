/*###############################
## Author:  Rajdeep Mukherjee ##
## Date:    8/3/2017          ##
## Task:    CRC Generator     ##
###############################*/

#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

#define CRC_INITIAL_VALUE 0xFFFFFFFF;
#define CRC_REMAINDER 0xC704DD7B;

struct state_elements_main {
  unsigned int crc;
};
struct state_elements_main smain;

unsigned int parallel_crces(
// Because the equations are long I am keeping the name of
// the incoming CRC and the XOR vector short.
unsigned int c1,                // crc_in
unsigned char d)                // data_in
{
unsigned char x1;

// The first thing that a parallel CRC needs to do is to develop
// the vector formed by XORing the input vector by the current CRC.
// This vector is then used during the CRC calculation.

x1 = ((c1>>24)&0xFF) ^ d;

// decompose c1 into c array
_Bool c[32];
for(unsigned i=0;i<=31;i++)
  c[i]=((c1>>i)&0x1);

_Bool x[32];
// fill x from index 24 to 31
x[24] = x1&0x1;
x[25] = x1&0x2;
x[26] = x1&0x4;
x[27] = x1&0x8;
x[28] = x1&0x10;
x[29] = x1&0x20;
x[30] = x1&0x40;
x[31] = x1&0x80;

_Bool parallel_crc[32];

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
parallel_crc[13] = x[31] ^ x[25] ^ x[26] ^ x[27] ^ x[29] ^ x[30] ^ c[5];
parallel_crc[12] = x[30] ^ x[24] ^ x[25] ^ x[26] ^ x[28] ^ x[29] ^ c[4];
parallel_crc[11] = x[24] ^ x[25] ^ x[27] ^ x[28] ^ c[3];
parallel_crc[10] = x[24] ^ x[26] ^ x[27] ^ x[29] ^ c[2];
parallel_crc[9] = x[25] ^ x[26] ^ x[28] ^ x[29] ^ c[1];
parallel_crc[8] = x[24] ^ x[25] ^ x[27] ^ x[28] ^ c[0];
parallel_crc[7] = x[24] ^ x[26] ^ x[27] ^ x[29] ^ x[31];
parallel_crc[6] = x[25] ^ x[26] ^ x[28] ^ x[29] ^ x[30] ^ x[31];
parallel_crc[5] = x[31] ^ x[30] ^ x[29] ^ x[28] ^ x[27] ^ x[25] ^ x[24];
parallel_crc[4] = x[30] ^ x[28] ^ x[27] ^ x[26] ^ x[24];
parallel_crc[3] = x[31] ^ x[25] ^ x[26] ^ x[27];
parallel_crc[2] = x[30] ^ x[24] ^ x[31] ^ x[25] ^ x[26];
parallel_crc[1] = x[30] ^ x[24] ^ x[31] ^ x[25];
parallel_crc[0] = x[30] ^ x[24];

// conver the value in the array to a number
unsigned int n=0;
for(unsigned int i=0;i<=31;i++)
 n=n+2^i*parallel_crc[i];

 return n;
} // endfunction parallel_crc

void initial() {
  smain.crc = CRC_INITIAL_VALUE;
}

void crc_main(
    _Bool         clk,          // everything clocks on rising edge
    _Bool         clken,        // Clock enable
    _Bool         reset,        // RESET
    _Bool         load,         // Load the CRC register
    _Bool         compute,      // compute the next 8-bits of CRC
    unsigned char data_in,      // data in
    unsigned char  data_out,
    _Bool         *crc_ok,      // the last CRC check cycle was OK
    unsigned int *crc           // direct acces to the CRC value
    )
{
  unsigned int new_crc;

  // Concurrent signal assignment for the crc generation.
  new_crc = reset ? 0xFFFFFFFF :
    load ? (((smain.crc&0xFFFFFF) << 7) | data_in) :
    compute ? parallel_crces(smain.crc,data_in) :
    smain.crc;

  data_out = ~((smain.crc>>24) & 0xFF);


  // clocked block
  // Create the 32 Flip flops (with clock enable flops)
  if (clken)
   smain.crc = new_crc;

  // update the output
  *crc = smain.crc;

  *crc_ok = (smain.crc == 0xC704DD7B);

  // This is a 32-bit wide AND function, so proper attention should be paid
  // when synthesizing to achieve good results.  If there are cycles
  // available pipeling this gate would be appropriate.
}

int main()
{
  _Bool         clk;            // everything clocks on rising edge
  _Bool         clken;  // Clock enable
  _Bool         reset;  // RESET
  _Bool         load;   // Load the CRC register
  _Bool         compute;        // compute the next 8-bits of CRC
  unsigned char data_in;        // data in
  unsigned char  data_out;
  _Bool           crc_ok;       // the last CRC check cycle was OK
  unsigned int crc;             // direct acces to the CRC value


  unsigned int  new_crc;        // D inputs to the flops

  initial();
  //while(1) {
    crc_main(clk,clken,reset,load,compute,data_in,data_out,&crc_ok,&crc);

    // property
    assert(crc != 0xAAAAAAAA);
  //}
}
