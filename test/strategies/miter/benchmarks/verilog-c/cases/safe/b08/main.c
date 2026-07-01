#include <stdio.h>
#include <assert.h>

_Bool nondet_bool();
unsigned char nondet_char();

int start_st = 0;
int init = 1;
int loop_st = 2;
int the_end = 3;

struct state_elements_b08 {
    unsigned char        O;
    _Bool    STATO;
    unsigned char        IN_R;
    unsigned char        MAR;
    unsigned char        OUT_R;
};
struct state_elements_b08 sb;

unsigned char ROMS(unsigned char a)
{
  unsigned int ROM;
  if      (a ==  0) ROM = 0x7F97A;
  else if (a ==  1) ROM = 0x39A62;
  else if (a ==  2) ROM = 0xA8FFF;
  else if (a ==  3) ROM = 0xFF6BA;
  else if (a ==  4) ROM = 0xFFF6E;
  else if (a ==  5) ROM = 0xFFBA8;
  else if (a ==  6) ROM = 0xCA75B;
  else              ROM = 0x2FFF4;
  return ROM;
} // ROM

void initial()
{
        sb.STATO = start_st;
        sb.MAR = 0;
        sb.IN_R = 0;
        sb.OUT_R = 0;
        sb.O = 0;
}

void b08(
    _Bool        CLOCK,
    _Bool        START,
    unsigned char  I,
    unsigned char *O)
{

  unsigned char          ROM_1, ROM_2;
  unsigned char          ROM_OR;
  unsigned int  ROM_OUT;


        // clocked block
  switch (sb.STATO) {
    case start_st: {
                     if (START) sb.STATO = init;
                     break;
                   }
    case init: {
                 sb.IN_R  = I;
                 sb.OUT_R = 0;
                 sb.MAR   = 0;
                 sb.STATO = loop_st;
                 break;
               }
    case loop_st: {
                    if (((ROM_2 & ~sb.IN_R) | (ROM_1 & sb.IN_R) | (ROM_2 & ROM_1)) ==
                        0xFF) {
                      sb.OUT_R = (sb.OUT_R | ROM_OR);
                    }
                    sb.STATO = the_end;
                    break;
                  }
    case the_end: {
                    if (sb.MAR != 7) {
                      sb.MAR = sb.MAR + 1;
                      sb.STATO = loop_st;
                    } else if (!START) {
                      sb.O = sb.OUT_R;
                      sb.STATO = start_st;
                    }
                    break;
                  }
  }

  ROM_OUT = ROMS(sb.MAR);
  ROM_1 = ((ROM_OUT>>12)&0xFF);
  ROM_2 = ((ROM_OUT>>4)&0xFF);
  ROM_OR = (ROM_OUT&0xF);
  assert((ROM_OR&0xF)!=0);
  *O=sb.O;
}

void main()
{
    _Bool        CLOCK;
    _Bool        START;
    unsigned char  I;
    unsigned char O;
    initial();
    while(1) {
      START=nondet_bool();
      I=nondet_char();
      b08(CLOCK, START, I, &O);
    }
}
