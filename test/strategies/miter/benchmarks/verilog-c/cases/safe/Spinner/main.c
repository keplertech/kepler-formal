#include <stdio.h>
#include <assert.h>

struct state_elements_main{
  unsigned int inr;
  _Bool spl;
  unsigned int dout;
};
struct state_elements_main  smain;

void spinner(_Bool clock, _Bool spin, unsigned char amount, unsigned int din, unsigned int *dout)
{
  unsigned int tmp0;
  unsigned int tmp1;
  unsigned int tmp2;
  unsigned int tmp3;
  unsigned int tmp4;
  unsigned int tmp5;

  tmp0 = smain.inr;
  tmp1 = ((amount & 0x1) ? ((tmp0 >> 1) & 2147483647) : (tmp0 & 2147483647)) | ((tmp1 >> 31)<<31);
  tmp1 = (tmp1 & 0x7FFFFFFF) | (((amount&0x1) ? (tmp0&0x1) : ((tmp0>>31) & 1)) << 31);

  tmp2 = (((amount>>1) & 0x1) ? ((tmp1 >> 2) & 1073741823) : (tmp1 & 1073741823)) | ((tmp2 >> 30)<<30);
  tmp2 = (((((amount>>1) & 0x1) ? (tmp1 & 0x3) : ((tmp1>>30) & 0x3))&0x00000003)<<30) | tmp2&0x3FFFFFFF;

  tmp3 = (((amount>>2) & 0x1) ? ((tmp2 >> 4) & 0xFFFFFFF) : (tmp2 & 0xFFFFFFF)) | ((tmp3 >> 28)<<28);
  tmp3 = (((((amount>>2) & 0x1) ? (tmp2 & 0xF) : ((tmp2>>28) & 0xF))&0x0000000F)<<28) | tmp3&0x0FFFFFFF;

  tmp4 = (((amount>>3) & 0x1) ? ((tmp3 >> 8) & 0x0FFFFFF) : (tmp3 & 0x0FFFFFF)) | ((tmp4 >> 24)<<24);
  tmp4 = (((((amount>>3) & 0x1) ? (tmp3 & 0xFF) : ((tmp3>>24) & 0xFF))&0x000000FF)<<24) | tmp3&0x00FFFFFF;


  tmp5 = (((amount>>4) & 0x1) ? ((tmp4 >> 16) & 0xFFFF) : (tmp4 & 0xFFFF)) | ((tmp4 >> 16)<<16);
  tmp5 = (((((amount>>4) & 0x1) ? (tmp4 & 0xFFFF) : ((tmp4>>16) & 0xFFFF))&0x0000FFFF)<<16) | tmp4&0x0000FFFF;

  if(smain.spl)
    smain.inr = smain.dout;
  else
    smain.inr = din;

  smain.dout = tmp5;
  smain.spl = spin;

  // pass the output
  *dout = smain.dout;

  // PASS
  assert(!((smain.inr&0xFFFFFFFF!=0) && (smain.dout&0xFFFFFFFF!=0)) || (smain.spl==0 || smain.inr&0xFFFFFFFF!=0));
  // PASS
  //assert(!(inr[31:0]==0 && dout[31:0]==0) ||  (spl==0 || inr[31:0]==0));
  assert(!((smain.inr&0xFFFFFFFF==0) && (smain.dout&0xFFFFFFFF==0)) || (smain.spl==0 || smain.inr&0xFFFFFFFF==0));
}

void main() {
  _Bool clock;
  _Bool spin;
  unsigned char amount;
  unsigned int din;
  unsigned int dout;
  while(1) {
    spinner(clock, spin, amount, din, &dout);
  }
}
