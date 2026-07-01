#include <stdio.h>
#include <assert.h>

struct state_elements_barrel {
  unsigned char b0, r0;
  unsigned char b1, r1;
  unsigned char b2, r2;
  unsigned char b3, r3;
};
struct state_elements_barrel sb;

_Bool fvalid(
    unsigned char b0, unsigned char r0,
    unsigned char b1, unsigned char r1,
    unsigned char b2, unsigned char r2,
    unsigned char b3, unsigned char r3
    )
{
  _Bool valid;
  valid =
    ((b0&0x3 != r0&0x3)|| ( b1&0x3 == r1&0x3)) &&
    ((b0&0x3 != r1&0x3)|| ( b1&0x3 == r2&0x3)) &&
    ((b0&0x3 != r2&0x3)|| ( b1&0x3 == r3&0x3)) &&
    ((b0&0x3 != r3&0x3)|| ( b1&0x3 == r0&0x3)) &&
    ((b1&0x3 != r0&0x3)|| ( b2&0x3 == r1&0x3)) &&
    ((b1&0x3 != r1&0x3)|| ( b2&0x3 == r2&0x3)) &&
    ((b1&0x3 != r2&0x3)|| ( b2&0x3 == r3&0x3)) &&
    ((b1&0x3 != r3&0x3)|| ( b2&0x3 == r0&0x3)) &&
    ((b2&0x3 != r0&0x3)|| ( b3&0x3 == r1&0x3)) &&
    ((b2&0x3 != r1&0x3)|| ( b3&0x3 == r2&0x3)) &&
    ((b2&0x3 != r2&0x3)|| ( b3&0x3 == r3&0x3)) &&
    ((b2&0x3 != r3&0x3)|| ( b3&0x3 == r0&0x3)) &&
    ((b3&0x3 != r0&0x3)|| ( b0&0x3 == r1&0x3)) &&
    ((b3&0x3 != r1&0x3)|| ( b0&0x3 == r2&0x3)) &&
    ((b3&0x3 != r2&0x3)|| ( b0&0x3 == r3&0x3)) &&
    ((b3&0x3 != r3&0x3) || ( b0&0x3 == r0&0x3));
  return valid;
}

void initial()
{
  sb.b0 = 0&0x3; sb.r0 = 0&0x3;
  sb.b1 = 0&0x3; sb.r1 = 0&0x3;
  sb.b2 = 0&0x3; sb.r2 = 0&0x3;
  sb.b3 = 0&0x3; sb.r3 = 0&0x3;
  // Then remap invalid initial states to valid one.
  if (!fvalid(
        sb.b0&0x3, sb.r0&0x3,
        sb.b1&0x3, sb.r1&0x3,
        sb.b2&0x3, sb.r2&0x3,
        sb.b3&0x3, sb.r3&0x3
        ))
  {
    sb.b0 = 0; sb.r0 = 0;
    sb.b1 = 0; sb.r1 = 0;
    sb.b2 = 0; sb.r2 = 0;
    sb.b3 = 0; sb.r3 = 0;
  }
}

int barrel(_Bool clk)
{
  /*unsigned char b0_old=sb.b0;
  unsigned char b1_old=sb.b1;
  unsigned char b2_old=sb.b2;
  unsigned char b3_old=sb.b3;

  // always block
  sb.b0 = b1_old;
  sb.b1 = b2_old;
  sb.b2 = b3_old;
  sb.b3 = b0_old;*/

  sb.b0 = sb.b1;
  sb.b1 = sb.b2;
  sb.b2 = sb.b3;
  sb.b3 = sb.b0;
} // barrel

int main()
{
  _Bool clk;
  initial();
  barrel(clk);
  barrel(clk);
  assert
    ((((sb.b0&0x3 != sb.r0&0x3) || (sb.b1&0x3 == sb.r1&0x3)) &&
     ((sb.b0&0x3 != sb.r1&0x3) || ( sb.b1&0x3 == sb.r2&0x3)) &&
     ((sb.b0&0x3 != sb.r2&0x3) || ( sb.b1&0x3 == sb.r3&0x3)) &&
     ((sb.b0&0x3 != sb.r3&0x3) || ( sb.b1&0x3 == sb.r0&0x3)) &&
     ((sb.b1&0x3 != sb.r0&0x3) || ( sb.b2&0x3 == sb.r1&0x3)) &&
     ((sb.b1&0x3 != sb.r1&0x3) || ( sb.b2&0x3 == sb.r2&0x3)) &&
     ((sb.b1&0x3 != sb.r2&0x3) || ( sb.b2&0x3 == sb.r3&0x3)) &&
     ((sb.b1&0x3 != sb.r3&0x3) || ( sb.b2&0x3 == sb.r0&0x3)) &&
     ((sb.b2&0x3 != sb.r0&0x3) || ( sb.b3&0x3 == sb.r1&0x3)) &&
     ((sb.b2&0x3 != sb.r1&0x3) || ( sb.b3&0x3 == sb.r2&0x3)) &&
     ((sb.b2&0x3 != sb.r2&0x3) || ( sb.b3&0x3 == sb.r3&0x3)) &&
     ((sb.b2&0x3 != sb.r3&0x3) || ( sb.b3&0x3 == sb.r0&0x3)) &&
     ((sb.b3&0x3 != sb.r0&0x3) || ( sb.b0&0x3 == sb.r1&0x3)) &&
     ((sb.b3&0x3 != sb.r1&0x3) || ( sb.b0&0x3 == sb.r2&0x3)) &&
     ((sb.b3&0x3 != sb.r2&0x3) || ( sb.b0&0x3 == sb.r3&0x3)) &&
     ((sb.b3&0x3 != sb.r3&0x3) || ( sb.b0&0x3 == sb.r0&0x3))));
}
