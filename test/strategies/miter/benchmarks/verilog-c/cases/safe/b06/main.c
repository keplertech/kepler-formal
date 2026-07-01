#include <stdio.h>
#include <assert.h>

struct state_elements_b06 {
    unsigned char CC_MUX, USCITE;
    _Bool        ENABLE_COUNT, ACKOUT;
    _Bool    state;
};
struct state_elements_b06 sb;

int      cc_nop   = 1;
int      cc_enin  = 1;
int      cc_intr  = 2;
int      cc_ackin = 3;
int      out_norm = 1;

int s_init=0;
int s_wait=1;
int s_enin=2;
int s_enin_w=3;
int s_intr=4;
int s_intr_1=5;
int s_intr_w=6;

void initial()
{
  sb.state = s_init;
  sb.CC_MUX = 0;
  sb.ENABLE_COUNT = 0;
  sb.ACKOUT = 0;
  sb.USCITE = 0;
}

void b06(
    unsigned char *CC_MUX,
    _Bool        EQL,
    unsigned char *USCITE,
    _Bool        clock,
    _Bool  *ENABLE_COUNT,
    _Bool  *ACKOUT,
    _Bool  CONT_EQL
    )
{
  // clocked block
        if (CONT_EQL) {
            *ACKOUT = 0;
            *ENABLE_COUNT = 0;
  } else {
            *ACKOUT = 1;
            *ENABLE_COUNT = 1;
  }

        switch (sb.state)
  {

          case s_init: {
                   sb.CC_MUX = cc_enin;
                   sb.USCITE = out_norm;
                   sb.state = s_wait;
                   break;
                 }
    case s_wait: {
                   if (EQL) {
                     sb.USCITE = 0;
                     sb.CC_MUX = cc_ackin;
                     sb.state = s_enin;
                   } else {
                     sb.USCITE = out_norm;
                     sb.CC_MUX = cc_intr;
                     sb.state = s_intr_1;
                   }
                   break;
                 }
    case s_intr_1: {
                     if (EQL) {
                       sb.USCITE = 0;
                       sb.CC_MUX = cc_ackin;
                       sb.state = s_intr;
                     } else {
                       sb.USCITE = out_norm;
                       sb.CC_MUX = cc_enin;
                       sb.state = s_wait;
                     }
                     break;
                   }
    case s_enin: {
                   if (EQL) {
                     sb.USCITE = 0;
                     sb.CC_MUX = cc_ackin;
                     sb.state = s_enin;
                   } else {
                     sb.USCITE = 1;
                     sb.ACKOUT = 1;
                     sb.ENABLE_COUNT = 1;
                     sb.CC_MUX = cc_enin;
                     sb.state = s_enin_w;
                   }
                   break;
                 }
    case s_enin_w: {
                     if (EQL) {
                       sb.USCITE = 1;
                       sb.CC_MUX = cc_enin;
                       sb.state = s_enin_w;
                     } else {
                       sb.USCITE = out_norm;
                       sb.CC_MUX = cc_enin;
                       sb.state = s_wait;
                     }
                     break;
                   }
    case s_intr: {
                   if (EQL) {
                     sb.USCITE = 0;
                     sb.CC_MUX = cc_ackin;
                     sb.state = s_intr;
                   } else {
                     sb.USCITE = 3;
                     sb.CC_MUX = cc_intr;
                     sb.state = s_intr_w;
                   }
                   break;
                 }
    case s_intr_w: {
                     if (EQL) {
                       sb.USCITE = 3;
                       sb.CC_MUX = cc_intr;
                       sb.state = s_intr_w;
                     } else {
                       sb.USCITE = out_norm;
                       sb.CC_MUX = cc_enin;
                       sb.state = s_wait;
                     }
                   }
                   break;
  }
  *USCITE = sb.USCITE;
  *CC_MUX = sb.CC_MUX;
  *ENABLE_COUNT = sb.ENABLE_COUNT;
  *ACKOUT = sb.ACKOUT;
}

void main()
{
  unsigned char CC_MUX;
  _Bool  EQL;
  unsigned char USCITE;
  _Bool  clock;
  _Bool  ENABLE_COUNT;
  _Bool  ACKOUT;
  _Bool  CONT_EQL;
  initial();
  while (1) {
    b06(&CC_MUX, EQL, &USCITE, clock, &ENABLE_COUNT, &ACKOUT, CONT_EQL);
    assert (ENABLE_COUNT==ACKOUT);
    assert ((USCITE&0x3)!=2);
  }
}
