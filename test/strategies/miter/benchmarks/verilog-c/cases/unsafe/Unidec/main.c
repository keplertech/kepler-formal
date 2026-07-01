#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

struct state_elements_main {
  unsigned short int word;
  _Bool         found;
  _Bool         init;
};
struct state_elements_main smain;

unsigned short int codes(unsigned char sel)
{
  unsigned short int code;
  if(sel&0x7 == 0)
    code = 8;

  else
    if(sel&0x7 == 1)
      code = 10;

    else
      if(sel&0x7 == 2)
        code = 0x58;

      else
        if(sel&0x7 == 3)
          code = 0x248;

        else
          if(sel&0x7 == 4)
            code = 0x2C1;

          else
            if(sel&0x7 == 5)
              code = 0x263;

            else
              if(sel&0x7 == 6)
                code = 0xC689;

              else
                if(sel&0x7 == 7)
                  code = 0x8;

  return code;
}

unsigned short int suffixes(unsigned short int word, unsigned char sel)
{
  unsigned short int suffix;
  if(sel&0x3 == 0)
    suffix = ((0x0&0x7)<<13) | ((smain.word>>3)&0x1FFF);

  else
    if(sel&0x3 == 1)
    suffix = ((0x0&0x3F)<<10) | ((smain.word>>6)&0x3FF);


    else
      if(sel&0x3 == 2)
    suffix = ((0x0&0x1FF)<<7) | ((smain.word>>9)&0x7F);

      else
        if(sel&0x3 == 3)
    suffix = ((0x0&0xFFF)<<4) | ((smain.word>>12)&0xF);

  return suffix;
}

unsigned short int prefixes(unsigned short int word, unsigned char sel)
{
  unsigned short int prefix;
  if(sel&0x3 == 0) {
    if(((word>>4)&0xFFF)==0)
      prefix=0x7FFF;
    else
      prefix=(0x1FFF<<2) | (word&0x7);
  }
  if(sel&0x3 == 1) {
    if(((word>>7)&0x1FF)==0)
      prefix=0x7FFF;
    else
      prefix=(0x3FF<<5) | (word&0x3F);
  }
  if(sel&0x3 == 2) {
    if(((word>>10)&0x3F)==0)
      prefix=0x7FFF;
    else
      prefix=(0x7F<<8) | (word&0x1FF);
  }
  if(sel&0x3 == 3) {
    if(((word>>13)&0x7)==0)
      prefix=0x7FFF;
    else
      prefix=(0xF<<11) | (word&0x7FF);
  }
  return prefix;
}

void initial(unsigned char sel1)
{
  smain.word = codes(sel1);
  smain.found = 0;
  smain.init = 1;
}

void unidec(
    _Bool       clk,
    unsigned char sel1,
    unsigned char sel2,
    _Bool *found)
{
  unsigned short int other;
  other = codes(sel1);

  smain.init = 0;
  if (other == prefixes(smain.word,sel2)) {
    // There is a code word that is a prefix of the current word.
    // Make the suffix of the current word the next word.
    smain.word = suffixes(smain.word,sel2);
  }
  else if (prefixes(other,sel2) == smain.word) {
    // The current word is a prefix of another code word.
    // Make the suffix of the other word the next word.
    smain.word = suffixes(other,sel2);
  }
  else {
    // Neither applies.  Go to trap state.
    smain.word = 0;
  }

  // clocked block
  smain.found = !smain.init && smain.word == other;
  *found = !smain.init && smain.word == other;
}

int main()
{
  _Bool clk;
  unsigned char sel1;
  unsigned char sel2;
  _Bool found;
  //while(1) {
    unidec(clk, sel1, sel2, &found);

    //#PASS: The trap state lives up to its name.
    //assert((smain.word!=0) || (smain.word==0));

    // FAIL
    assert(found==1);
  //}
}
