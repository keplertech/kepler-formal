#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

struct state_elements_main{
unsigned int dout;
unsigned int inr;
};
struct state_elements_main  smain;

void rotate(_Bool clock, unsigned char amount, unsigned int din, unsigned int *dout)
{
  unsigned int dout_old;
  unsigned int inr_old;
  unsigned int tmp0;
  unsigned int tmp1;
  unsigned int tmp2;
  unsigned int tmp3;
  unsigned int tmp4;
  unsigned int tmp5;

  // clocked block
  smain.inr = din;

  // continuous assignment
  tmp0 = smain.inr;
  tmp1 = (amount&0x1) ? (tmp0 >> 1)|(tmp0 << (32 - 1)) : tmp0;
  tmp2 = (amount&0x2) ? (tmp1 >> 2)|(tmp1 << (32 - 2)) : tmp1;
  tmp3 = (amount&0x4) ? (tmp2 >> 4)|(tmp2 << (32 - 4)) : tmp2;
  tmp4 = (amount&0x8) ? (tmp3 >> 8)|(tmp3 << (32 - 8)) : tmp3;
  tmp5 = (amount&0x10) ? (tmp4 >> 16)|(tmp4 << (32 - 16)) : tmp4;

  // clocked
  smain.dout = tmp5;
  *dout = tmp5;
}

void initial_main(){
  smain.dout = 0;
  smain.inr = 0;
}

void main() {
  _Bool clock;
  unsigned char amount;
  unsigned int din;
  unsigned int dout;
  initial_main();
  //while(1) {
    rotate(clock, amount, din, &dout);
    // # FAIL:
    assert(dout != 0xAAAAAAAA);
  //}
}
