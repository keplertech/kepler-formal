#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

struct state_elements_control {
  unsigned char prev;
  _Bool         Lcmd;
  _Bool         Ucmd;
  _Bool         Ncmd;
  _Bool         Ccmd;
};
struct state_elements_control scontrol;

void initial_control()
{
  scontrol.Lcmd = 0;
  scontrol.Ucmd = 0;
  scontrol.Ncmd = 1;
  scontrol.Ccmd = 0;
  scontrol.prev = 0;
}

void control(
    _Bool       clock,
    _Bool       reset,
    unsigned char       in,
    _Bool       *Lcmd,
    _Bool       *Ucmd,
    _Bool       *Ncmd,
    _Bool       *Ccmd)
{

  _Bool load;

  // clocked block
  if (reset)
    scontrol.prev = 0;
  else
    scontrol.prev = in;

  load = (scontrol.prev == 0x1b); // escape

  if (reset) {
    scontrol.Ncmd = 1;
    scontrol.Lcmd = 0; scontrol.Ucmd = 0; scontrol.Ccmd = 0;
  }
  else if (load) {
    switch (in) {
      case 0x4c: { // L
                   scontrol.Lcmd = 1; scontrol.Ucmd = 0; scontrol.Ncmd = 0; scontrol.Ccmd = 0;
                   break;
                 }
      case 0x55: { // U
                   scontrol.Lcmd = 0; scontrol.Ucmd = 1; scontrol.Ncmd = 0; scontrol.Ccmd = 0;
                   break;
                 }
      case 0x4e: { // N
                   scontrol.Lcmd = 0; scontrol.Ucmd = 0; scontrol.Ncmd = 1; scontrol.Ccmd = 0;
                   break;
                 }
      case 0x43: { // C
                   scontrol.Lcmd = 0; scontrol.Ucmd = 0; scontrol.Ncmd = 0; scontrol.Ccmd = 1;
                   break;
                 }
      default: {
                 scontrol.Lcmd = 0; scontrol.Ucmd = 0; scontrol.Ncmd = 0; scontrol.Ccmd = 0;
                 break;
               }
    }
  } // if (load)
  *Lcmd = scontrol.Lcmd;
  *Ucmd = scontrol.Ucmd;
  *Ncmd = scontrol.Ncmd;
  *Ccmd = scontrol.Ccmd;

} // control

// Utility functions
_Bool isUpper(
        unsigned char   in)
{
  _Bool isupper;
  isupper = ~((in>>5)&1);
  return isupper;
}

unsigned char toLower(
        unsigned char   in)
{
  unsigned char tolower;
        if (isUpper(in))
            tolower = in + 0x20;
        else
            tolower = in;
  return tolower;
}

unsigned char toUpper(
        unsigned char   in)
{
        unsigned char toupper;
        if (!isUpper(in))
            toupper = in - 0x20;
        else
            toupper = in;
  return toupper;
}

unsigned char changeCase(
        unsigned char   in)
{
        unsigned char changecase;
  return changecase;
  if (isUpper(in))
            changecase = in + 0x20;
        else
            changecase = in - 0x20;
  return changecase;
}

void transform(
    unsigned char in,
    _Bool        Lcmd,
    _Bool        Ucmd,
    _Bool        Ncmd,
    _Bool        Ccmd,
    unsigned char *out)
{
    *out = Lcmd ? toLower(in) :
           Ucmd ? toUpper(in) :
           Ncmd ? in :
           Ccmd ? changeCase(in) : 0x00;
}

struct state_elements_main {
  unsigned char regIn;
  unsigned char dataOut;
};
struct state_elements_main smain;

void initial() {
        smain.regIn = 0;
        smain.dataOut = 0;
}

void vlunc(
    _Bool        clock,
    _Bool        reset,
    unsigned char        dataIn,
    unsigned char *dataOut)
{
    unsigned char transformed;
    _Bool        Lcmd, Ucmd, Ncmd, Ccmd;

    control(clock,reset,smain.regIn,&Lcmd,&Ucmd,&Ncmd,&Ccmd);

    transform(smain.regIn,Lcmd,Ucmd,Ncmd,Ccmd,&transformed);

    // clocked block
        if (reset) {
            smain.dataOut = 0;
            smain.regIn = 0;
  } else {
            smain.dataOut = transformed;
            smain.regIn = dataIn;
  } // else: !if(reset)

  *dataOut=smain.dataOut;
  //#PASS: The command lines are 1-hot encoded.
  assert(!(Lcmd==1 && Ucmd==1 || Lcmd==1 && Ccmd==1 || Lcmd==1 && Ncmd==1 || Ucmd==1 && Ccmd==1 || Ucmd==1 && Ncmd==1 || Ccmd==1 && Ncmd==1));
}

void main()
{
  _Bool  clock;
   _Bool         reset;
  unsigned char  dataIn;
  unsigned char dataOut;
  initial();
  initial_control();
  while(1)
    vlunc(clock,reset,dataIn,&dataOut);
}
