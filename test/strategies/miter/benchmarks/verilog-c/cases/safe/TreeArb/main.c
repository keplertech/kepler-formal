#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

int I1 = 0;
int I2 = 1;
int R1 = 2;
int R2 = 3;
int A1 = 4;
int A2 = 5;

struct state_elements_cell{
  _Bool state;
};

void initial_cell(struct state_elements_cell *scell)
{
   scell->state=I1;
}

void cell(struct state_elements_cell *scell, _Bool clk, _Bool ack0, _Bool req1, _Bool req2, _Bool *req0, _Bool *ack1, _Bool *ack2)
{
  if(scell->state == I1) {
   if(req1)
     scell->state = R1;
   else if(req2)
     scell->state = R2;
  }
  else if(scell->state == R1) {
   if(ack0)
     scell->state = A1;
  }
  else if(scell->state == A1) {
   if(!req1)
     scell->state = I2;
  }
  else if(scell->state == I2) {
    if(req2)
      scell->state = R2;
    else if(req1)
      scell->state = R1;
  }
  else if(scell->state == R2) {
    if(ack0)
      scell->state = A2;
  }
  else if(scell->state == A2) {
    if(!req2)
      scell->state = I1;
  }

  *req0 = (scell->state == R1) || (scell->state == R2) || (scell->state == A1) || (scell->state == A2);
  *ack1 = (scell->state == A1);
  *ack2 = (scell->state == A2);
}

struct state_elements_proc{
  _Bool state;
};

void initial_proc(struct state_elements_cell *sproc)
{
   sproc->state=0;
}

void proc(struct state_elements_proc *sproc, _Bool clk, _Bool ack, _Bool choice, _Bool *req)
{
  int idle = 0;
  int requesting = 1;
  int locking = 2;
  if(sproc->state == idle) {
   if(choice)
     sproc->state = requesting;
  }
  else if(sproc->state == requesting) {
    if(ack)
      sproc->state = locking;
  }
  else if(sproc->state == locking) {
    if(choice)
      sproc->state = idle;
  }

  *req = (sproc->state == requesting) || (sproc->state == locking);
}

void main()
{
  _Bool clk;
  _Bool a3;
  _Bool a2[2];
  _Bool a1[4];
  _Bool a0[8];
  _Bool r3;
  _Bool r2[2];
  _Bool r1[4];
  _Bool r0[8];
  _Bool choice[8];


  struct state_elements_cell c0;
  struct state_elements_cell c1;
  struct state_elements_cell c2;
  struct state_elements_cell c3;
  struct state_elements_cell c4;
  struct state_elements_cell c5;
  struct state_elements_cell c10;

  initial_cell(&c10);
  initial_cell(&c0);
  initial_cell(&c1);
  initial_cell(&c2);
  initial_cell(&c3);
  initial_cell(&c4);
  initial_cell(&c5);

  cell(&c10, clk, r2[0], r2[1], a2[0], &a2[1], &r3, &a3);
  a3 = r3;
  cell(&c0, clk, r1[0], r1[1], a1[0], &a1[1], &r2[0], &a2[0]);
  cell(&c1, clk, r1[2], r1[3], a1[2], &a1[3], &r2[1], &a2[1]);
  cell(&c2, clk, r0[0], r0[1], a0[0], &a0[1], &r1[0], &a1[0]);
  cell(&c3, clk, r0[2], r0[3], a0[2], &a0[3], &r1[1], &a1[1]);
  cell(&c4, clk, r0[4], r0[5], a0[4], &a0[5], &r1[2], &a1[2]);
  cell(&c5, clk, r0[6], r0[7], a0[6], &a0[7], &r1[3], &a1[3]);

  struct state_elements_proc p0;
  struct state_elements_proc p1;
  struct state_elements_proc p2;
  struct state_elements_proc p3;
  struct state_elements_proc p4;
  struct state_elements_proc p5;
  struct state_elements_proc p6;
  struct state_elements_proc p7;

  initial_proc(&p0);
  initial_proc(&p1);
  initial_proc(&p2);
  initial_proc(&p3);
  initial_proc(&p4);
  initial_proc(&p5);
  initial_proc(&p6);
  initial_proc(&p7);

  proc(&p0, clk, a0[0], r0[0], &choice[0]);
  proc(&p1, clk, a0[1], r0[1], &choice[1]);
  proc(&p2, clk, a0[2], r0[2], &choice[2]);
  proc(&p3, clk, a0[3], r0[3], &choice[3]);
  proc(&p4, clk, a0[4], r0[4], &choice[4]);
  proc(&p5, clk, a0[5], r0[5], &choice[5]);
  proc(&p6, clk, a0[6], r0[6], &choice[6]);
  proc(&p7, clk, a0[7], r0[7], &choice[7]);

  assert(!((p0.state==2 || p1.state==2 || p2.state==2 || p3.state==2) &&
 (p4.state==2 || p5.state==2 || p6.state==2 || p7.state==2) ||
 (p0.state==2 || p1.state==2 || p4.state==2 || p5.state==2) &&
 (p2.state==2 || p3.state==2 || p6.state==2 || p7.state==2) ||
 (p0.state==2 || p2.state==2 || p4.state==2 || p6.state==2) &&
 (p1.state==2 || p3.state==2 || p5.state==2 || p7.state==2)));
}
