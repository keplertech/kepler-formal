#include <stdio.h>
#define TRUE 1
#define FALSE 0

struct state_elements_main {
  _Bool ack0;
  _Bool ack1;
  _Bool req0;
  _Bool req1;
  _Bool robin;
};
struct state_elements_main  smain;

void rrobin(_Bool clock, _Bool ir0, _Bool ir1, _Bool *ack0, _Bool *ack1)
{
  if(!smain.req0)
    smain.ack0 = FALSE;
  else if(!smain.req1)
    smain.ack0 = TRUE;
  else if(!smain.ack0 && !smain.ack1)
    smain.ack0 = !smain.robin;
  else
    smain.ack0 = !smain.ack0;

  if(!smain.req1)
    smain.ack1 = FALSE;
  else if(!smain.req0)
    smain.ack1 = TRUE;
  else if(!smain.ack0 && !smain.ack1)
    smain.ack1 = smain.robin;
  else
    smain.ack1 = !smain.ack1;

  if(smain.req0 && smain.req1 && !smain.ack0 && !smain.ack1)
    smain.robin = !smain.robin;

  smain.req0 = ir0;
  smain.req1 = ir1;

  // update the acknowledgement
  *ack0 = smain.ack0;
  *ack1 = smain.ack1;
}

void initial_main(_Bool ir0, _Bool ir1) {
  smain.ack0 = FALSE;
  smain.ack1 = FALSE;
  smain.robin = FALSE;
  smain.req0 = ir0;
  smain.req1 = ir1;
}

void main() {
  _Bool clock;
  _Bool ir0;
  _Bool ir1;
  _Bool ack0;
  _Bool ack1;
  _Bool clock;
  _Bool ir0;
  _Bool ir1;
  _Bool ack0;
  _Bool ack1;
  initial_main(ir0, ir1);
  //while(1) {
    // All property passes
    /*if(smain.req0==0 && smain.req1==0)
    {
      rrobin(clock, ir0, ir1, &ack0, &ack1);
      assert(ack0==0 && ack1==0);
    }
    if(smain.req0==1 && smain.req1==0) {
      rrobin(clock, ir0, ir1, &ack0, &ack1);
      assert(ack0==1);
    }
    if(smain.req0==0 && smain.req1==1) {
      rrobin(clock, ir0, ir1, &ack0, &ack1);
      assert(ack1==1);
    }
    if(smain.req0==1 && smain.ack1==1) {
      rrobin(clock, ir0, ir1, &ack0, &ack1);
      assert(ack0==1);
    }*/
    if(smain.req1==1 && smain.ack0==1) {
      rrobin(clock, ir0, ir1, &ack0, &ack1);
      assert(ack1==1);
    }
  //}
}
