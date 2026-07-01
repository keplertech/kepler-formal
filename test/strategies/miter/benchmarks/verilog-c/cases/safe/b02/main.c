#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0
_Bool nondet_bool();

int A=0;
int B=1;
int C=2;
int D=3;
int E=4;
int F=5;
int G=6;

struct state_elements_main {
  _Bool U;
  unsigned char stato;
};
struct state_elements_main smain;

void initial() {
  smain.U = 0;
  smain.stato = A;
}

void b02(
    _Bool clock,
    _Bool  LINEA, _Bool *U)
{
  // clocked block
  switch(smain.stato) {
    case 0:
              smain.stato = B;
              smain.U=0;
              break;

    case 1:
              if (LINEA==0)
                smain.stato = C;
              else
                smain.stato = F;
              smain.U = 0;
              break;

    case 2:
              if (LINEA==0)
                smain.stato = D;
              else
                smain.stato = G;
              smain.U = 0;
              break;

    case 3:
              smain.stato = E;
              smain.U=0;
              break;

    case 4:
              smain.stato = B;
              smain.U=1;
              break;

    case 5:
              smain.stato = G;
              smain.U=0;
              break;

    case 6:
              if (LINEA==0)
                smain.stato = E;
              else
                smain.stato = A;
              smain.U = 0;
              break;

  }
  *U=smain.U;
}

void main()
{
  _Bool LINEA=nondet_bool();
  _Bool U=nondet_bool;
  _Bool  clock;
  initial();
  assert(smain.stato==A);

  while(1) {
    b02(clock,LINEA,&U);
    if(smain.stato==D) {
      b02(clock,LINEA,&U);
      b02(clock,LINEA,&U);
      assert(smain.stato==B);
    }
    /*if(smain.U==1) {
      b02(clock,LINEA,&U);
      assert(smain.U==0);
    }
    if(smain.U==1) {
      assert(smain.stato==B);
    }*/
  }
}
