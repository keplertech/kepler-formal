#include <stdio.h>
#include <assert.h>

_Bool nondet_bool();

int S_RESET = 0;
int S_START = 1;
int S_LOAD_X = 2;
int S_UPDATE_MAR = 3;
int S_LOAD_Y = 4;
int S_CALC_RETTA = 5;
int S_INCREMENTA = 6;

struct state_elements_b07 {
    unsigned char        PUNTI_RETTA;
    _Bool    stato;
    unsigned char        cont, x;
    unsigned char        y, t;
    unsigned char        mar;
};
struct state_elements_b07 sb;

unsigned char mems(unsigned char a)
{
  unsigned char mem;
  if      (a ==  0) mem = 1;
  else if (a ==  1) mem = 255;
  else if (a ==  2) mem = 0;
  else if (a ==  3) mem = 0;
  else if (a ==  4) mem = 0;
  else if (a ==  5) mem = 2;
  else if (a ==  6) mem = 0;
  else if (a ==  7) mem = 0;
  else if (a ==  8) mem = 0;
  else if (a ==  9) mem = 2;
  else if (a == 10) mem = 255;
  else if (a == 11) mem = 5;
  else if (a == 12) mem = 0;
  else if (a == 13) mem = 2;
  else if (a == 14) mem = 0;
  else              mem = 2;
  return mem;
} // mem

void initial()
{
        sb.stato = S_RESET;
        sb.PUNTI_RETTA = 0;
        sb.cont = 0;
        sb.mar = 0;
        sb.x = 0;
        sb.y = 0;
        sb.t = 0;
}

void b07(
    unsigned char *PUNTI_RETTA,
    _Bool        START,
    _Bool        clock)
{
   unsigned char mem_mar;
   mem_mar = mems(sb.mar);

  // clocked block
   switch (sb.stato) {
     case S_RESET: {
                     sb.stato = S_START;
                     break;
                   }
     case S_START: {
                     if (START) {
                       sb.cont = 0;
                       sb.mar = 0;
                       sb.stato = S_LOAD_X;
                     } else {
                       sb.stato = S_START;
                       sb.PUNTI_RETTA = 0;
                     }
                     break;
                   }
     case S_LOAD_X: {
                      sb.x = mem_mar;
                      sb.stato = S_UPDATE_MAR;
                      break;
                    }
     case S_UPDATE_MAR: {
                          sb.mar = sb.mar + 1;
                          sb.t = (((sb.x&0x3F)<<1) | 0)&0x7F;
                          sb.stato = S_LOAD_Y;
                          break;
                        }
     case S_LOAD_Y: {
                      sb.y = mem_mar&0x7F;
                      sb.x =  (((0<<7) | (sb.x&0x7F))&0xFF) + (((0<<7) | (sb.t&0x7F))&0xFF);
                      sb.stato = S_CALC_RETTA;
                      break;
                    }
S_CALC_RETTA: {
                sb.x = (((0<<7) | (sb.x&0x7F))&0xFF) + (((0<<7) | (sb.t&0x7F))&0xFF);
                sb.stato = S_INCREMENTA;
                break;
              }
S_INCREMENTA: {
                if (sb.mar != 15) {
                  if (sb.x == 2)
                    sb.cont = sb.cont + 1;
                  sb.mar = sb.mar + 1;
                  sb.stato = S_LOAD_X;
                } else {
                  if (START == 0) {
                    if (sb.x == 2)
                      sb.PUNTI_RETTA = sb.cont + 1;
                    else
                      sb.PUNTI_RETTA = sb.cont;
                    sb.stato = S_START;
                  } else {
                    sb.stato = S_INCREMENTA;
                  }
                }
                break;
              }
   }
   *PUNTI_RETTA = sb.PUNTI_RETTA;
}

void main()
{
    unsigned char PUNTI_RETTA;
    _Bool        START;
    _Bool        clock;
    initial();
    while(1) {
      START=nondet_bool();
      b07(&PUNTI_RETTA, START, clock);
      assert ((sb.x&0xFF)!=148);
    }
}
