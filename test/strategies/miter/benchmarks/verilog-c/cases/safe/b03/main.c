#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

int U1=0;
int U2=1;
int U3=2;
int U4=3;
int INIT=0;
int ANALISI_REQ=1;
int ASSEGNA=2;

struct state_elements_main {
  unsigned char  GRANT_O;
  unsigned char  stato;
  unsigned char  coda0, coda1, coda2, coda3;
  _Bool          ru1, ru2, ru3, ru4, fu1, fu2, fu3, fu4;
  unsigned char grant;
};
struct state_elements_main smain;

void initial() {
  smain.stato = INIT;
  smain.coda0 = 0;
  smain.coda1 = 0;
  smain.coda2 = 0;
  smain.coda3 = 0;
  smain.ru1 = 0;
  smain.fu1 = 0;
  smain.ru2 = 0;
  smain.fu2 = 0;
  smain.ru3 = 0;
  smain.fu3 = 0;
  smain.ru4 = 0;
  smain.fu4 = 0;
  smain.grant = 0;
  smain.GRANT_O = 0;
}

void b03(
    _Bool clock,
    _Bool REQUEST1,
    _Bool REQUEST2,
    _Bool REQUEST3,
    _Bool REQUEST4,
    unsigned char *GRANT_O
    )
{
  // clocked block
  switch (smain.stato) {
    case 0: {
              smain.ru1 = REQUEST1;
              smain.ru2 = REQUEST2;
              smain.ru3 = REQUEST3;
              smain.ru4 = REQUEST4;
              smain.stato = ANALISI_REQ;
              break;
            }
    case 1: {
              smain.GRANT_O = smain.grant;
              if (smain.ru1) {
                if (!smain.fu1) {
                  smain.coda3 = smain.coda2;
                  smain.coda2 = smain.coda1;
                  smain.coda1 = smain.coda0;
                  smain.coda0 = U1;
                }
              }
              else if (smain.ru2) {
                if (!smain.fu2) {
                  smain.coda3 = smain.coda2;
                  smain.coda2 = smain.coda1;
                  smain.coda1 = smain.coda0;
                  smain.coda0 = U2;
                }
              }
              else if (smain.ru3) {
                if (!smain.fu3) {
                  smain.coda3 = smain.coda2;
                  smain.coda2 = smain.coda1;
                  smain.coda1 = smain.coda0;
                  smain.coda0 = U3;
                }
              }
              else if (smain.ru4) {
                if (!smain.fu4) {
                  smain.coda3 = smain.coda2;
                  smain.coda2 = smain.coda1;
                  smain.coda1 = smain.coda0;
                  smain.coda0 = U4;
                }
              }
              smain.fu1 = smain.ru1;
              smain.fu2 = smain.ru2;
              smain.fu3 = smain.ru3;
              smain.fu4 = smain.ru4;
              smain.stato = ASSEGNA;
              break;
            }
    case 2: {
              if (smain.fu1 | smain.fu2 | smain.fu3 | smain.fu4) {
                switch (smain.coda0) {
                  case 4:       smain.grant = 8;
                  case 2: smain.grant = 4;
                  case 1: smain.grant = 2;
                  case 7: smain.grant = 1;
                  default: smain.grant = 0;
                }
                smain.coda0 = smain.coda1;
                smain.coda1 = smain.coda2;
                smain.coda2 = smain.coda3;
                smain.coda3 = 0;
              }
              smain.ru1 = REQUEST1;
              smain.ru2 = REQUEST2;
              smain.ru3 = REQUEST3;
              smain.ru4 = REQUEST4;
              smain.stato = ANALISI_REQ;
              break;
            }

  } // end switch
  *GRANT_O = smain.GRANT_O;
}

void main()
{
  _Bool    clock;
  _Bool          REQUEST1, REQUEST2, REQUEST3, REQUEST4;
  unsigned char GRANT_O;
  initial();
  while(1) {
    b03(clock,REQUEST1, REQUEST2, REQUEST3, REQUEST4, &GRANT_O);
    assert(GRANT_O==0 || GRANT_O==8 || GRANT_O==4 || GRANT_O==2 || GRANT_O==1);
  }
}
