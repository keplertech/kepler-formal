#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

struct state_elements_main{
  unsigned char b[8];
  unsigned char freg;
  unsigned char treg;
};
struct state_elements_main smain;

void initial_main(){
  smain.treg = 0;
  smain.freg = 0;
}

void four_by_four(_Bool clock, unsigned char from, unsigned char to)
{
  _Bool valid;
  _Bool parity;
  _Bool permutation;
  _Bool oddInversions;

  smain.freg = from;
  smain.treg = to;


  valid = (smain.b[smain.treg] == 0) &&
    (((smain.treg&0x3 == smain.freg&0x3) && !(smain.treg&0x4 == smain.freg&0x4)) ||
     ((smain.treg&0x4 == smain.freg&0x4) &&
      (((smain.treg&0x3 == 0) && (smain.freg&0x3 == 1)) ||
       ((smain.treg&0x3 == 1) && (smain.freg&0x1 == 0)) ||
       ((smain.treg&0x3 == 2) && (smain.freg&0x1 == 1)) ||
       ((smain.treg&0x3 == 3) && (smain.freg&0x3 == 2)) )));

  parity =
    (((smain.b[0] & 5) == 1) | ((smain.b[0] & 5) == 4)) ^
    (((smain.b[1] & 5) == 0) | ((smain.b[1] & 5) == 5)) ^
    (((smain.b[2] & 5) == 1) | ((smain.b[2] & 5) == 4)) ^
    (((smain.b[3] & 5) == 0) | ((smain.b[3] & 5) == 5)) ^
    (((smain.b[4] & 5) == 1) | ((smain.b[4] & 5) == 4)) ^
    (((smain.b[5] & 5) == 0) | ((smain.b[5] & 5) == 5)) ^
    (((smain.b[6] & 5) == 1) | ((smain.b[6] & 5) == 4)) ^
    (((smain.b[7] & 5) == 0) | ((smain.b[7] & 5) == 5));

  permutation = smain.b[0]!=smain.b[1] && smain.b[0]!=smain.b[2] && smain.b[0]!=smain.b[3] &&
    smain.b[0]!=smain.b[4] && smain.b[0]!=smain.b[5] && smain.b[0]!=smain.b[6] && smain.b[0]!=smain.b[7] &&
    smain.b[1]!=smain.b[2] && smain.b[1]!=smain.b[3] && smain.b[1]!=smain.b[4] && smain.b[1]!=smain.b[5] &&
    smain.b[1]!=smain.b[6] && smain.b[1]!=smain.b[7] && smain.b[2]!=smain.b[3] && smain.b[2]!=smain.b[4] &&
    smain.b[2]!=smain.b[5] && smain.b[2]!=smain.b[6] && smain.b[2]!=smain.b[7] && smain.b[3]!=smain.b[4] &&
    smain.b[3]!=smain.b[5] && smain.b[3]!=smain.b[6] && smain.b[3]!=smain.b[7] && smain.b[4]!=smain.b[5] &&
    smain.b[4]!=smain.b[6] && smain.b[4]!=smain.b[7] && smain.b[5]!=smain.b[6] && smain.b[5]!=smain.b[7] &&
    smain.b[6]!=smain.b[7];

  oddInversions = (smain.b[4]==0 || smain.b[5]==0 || smain.b[6]==0 || smain.b[7]==0) ^
    (smain.b[0]>smain.b[1] && smain.b[1]!=0) ^ (smain.b[0]>smain.b[2] && smain.b[2]!=0) ^
    (smain.b[0]>smain.b[3] && smain.b[3]!=0) ^ (smain.b[0]>smain.b[4] && smain.b[4]!=0) ^
    (smain.b[0]>smain.b[5] && smain.b[5]!=0) ^ (smain.b[0]>smain.b[6] && smain.b[6]!=0) ^
    (smain.b[0]>smain.b[7] && smain.b[7]!=0) ^ (smain.b[1]>smain.b[2] && smain.b[2]!=0) ^
    (smain.b[1]>smain.b[3] && smain.b[3]!=0) ^ (smain.b[1]>smain.b[4] && smain.b[4]!=0) ^
    (smain.b[1]>smain.b[5] && smain.b[5]!=0) ^ (smain.b[1]>smain.b[6] && smain.b[6]!=0) ^
    (smain.b[1]>smain.b[7] && smain.b[7]!=0) ^ (smain.b[2]>smain.b[3] && smain.b[3]!=0) ^
    (smain.b[2]>smain.b[4] && smain.b[4]!=0) ^ (smain.b[2]>smain.b[5] && smain.b[5]!=0) ^
    (smain.b[2]>smain.b[6] && smain.b[6]!=0) ^ (smain.b[2]>smain.b[7] && smain.b[7]!=0) ^
    (smain.b[3]>smain.b[4] && smain.b[4]!=0) ^ (smain.b[3]>smain.b[5] && smain.b[5]!=0) ^
    (smain.b[3]>smain.b[6] && smain.b[6]!=0) ^ (smain.b[3]>smain.b[7] && smain.b[7]!=0) ^
    (smain.b[4]>smain.b[5] && smain.b[5]!=0) ^ (smain.b[4]>smain.b[6] && smain.b[6]!=0) ^
    (smain.b[4]>smain.b[7] && smain.b[7]!=0) ^ (smain.b[5]>smain.b[6] && smain.b[6]!=0) ^
    (smain.b[5]>smain.b[7] && smain.b[7]!=0) ^ (smain.b[6]>smain.b[7] && smain.b[7]!=0);

  if(valid)
  {
    smain.b[smain.treg] = smain.b[smain.freg];
    smain.b[smain.freg] = 0;
  }

  //#PASS:
  assert(!permutation==1 || parity==0);
  /*
  //#PASS:
  assert property (permutation==1 |-> ##[0:100] permutation==1);
  //#PASS:
  assert property (permutation==1 && oddInversions==0 |-> ##[0:100] oddInversions==0);
  //#PASS:
  assert property (parity==0 |-> ##[0:100] parity==0);
  */
}

void main() {
  _Bool clock;
  unsigned char from;
  unsigned char to;
  initial_main();
  //while(1) {
    four_by_four(clock, from, to);
  //}
}
