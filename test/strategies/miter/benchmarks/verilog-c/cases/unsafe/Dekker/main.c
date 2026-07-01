#include <stdio.h>
#define TRUE 1
#define FALSE 0

_Bool nondet_bool();

struct state_elements_main{
  _Bool c[2];
  _Bool turn;
  _Bool self;
  unsigned char pc[2];
};

int L0 = 0;
int L1 = 1;
int L2 = 2;
int L3 = 3;
int L4 = 4;
int L5 = 5;
int L6 = 6;


void initial(struct state_elements_main *smain, _Bool clock, _Bool select, _Bool pause) {
   smain->pc[0] = L0 & 1;
   smain->pc[1] = (L0 & 1);
   smain->c[0] = TRUE;
   smain->c[1] = TRUE;
   smain->self = select;
}

 void dekker(struct state_elements_main *smain, _Bool clock, _Bool select, _Bool pause)
 {
    smain->self = select;
    if(smain->pc[smain->self] == L0) {
      if(!pause)
        smain->pc[smain->self] = L1;
    }
    else
      if(smain->pc[smain->self] == L1)
      {
        smain->c[smain->self] = FALSE;
        smain->pc[smain->self] = L2;
      }
      else
        if(smain->pc[smain->self] == L2) {
          if(smain->c[!smain->self] == 1)
            smain->pc[smain->self] = L5;
          else
            smain->pc[smain->self] = L3;
        }
        else
          if(smain->pc[smain->self] == L3) {
            if(smain->turn == smain->self)
              smain->pc[smain->self] = L2;
            else
            {
              smain->c[smain->self] = TRUE;
              smain->pc[smain->self] = L4;
            }
          }
          else
            if(smain->pc[smain->self] == L4) {
              if(smain->turn == smain->self)
              {
                smain->c[smain->self] = FALSE;
                smain->pc[smain->self] = L2;
              }
            }
            else
              if(smain->pc[smain->self] == L5) {
                if(!pause)
                  smain->pc[smain->self] = L6;
             }
              else
                if(smain->pc[smain->self] == L6)
                {
                  smain->c[smain->self] = TRUE;
                  smain->turn = !smain->self;
                  smain->pc[smain->self] = L0;
                }

  }

void main() {
 _Bool clock;
 _Bool select;
 _Bool pause;
 struct state_elements_main smain;

 initial(&smain,clock, select, pause);
 select = nondet_bool();
 dekker(&smain,clock, select, pause);
 // Mutual exclusion property, both the processes can not be in the critical
 // section at the same time
 assert((smain.pc[0] == L5 && smain.pc[1] == L5));
 // fair scheduling
 //assert((smain.self == 1 || smain.self == 0));
}
