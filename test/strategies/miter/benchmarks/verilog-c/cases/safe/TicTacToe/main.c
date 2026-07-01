#include <stdio.h>
#define TRUE 1
#define FALSE 0

unsigned char nondet_uchar();

struct state_elements_main{
 unsigned char b[9];
 unsigned char turn;
 unsigned char move;
};
struct state_elements_main  smain;

int EMPTY = 0;
int X = 1;
int O = 2;

void initial() {
  int i;
  for(i = 0; i < 9; i = i + 1)
    smain.b[i] = EMPTY;
  smain.turn = X;
  smain.move = 0;
}

void tictactoe(_Bool clock, unsigned char imove, _Bool *winX, _Bool *winO, _Bool *finished)
{
  smain.move = imove < 9 ? imove : 0;
  if(!finished && smain.b[smain.move] == EMPTY)
  {
    smain.b[smain.move] = smain.turn;
    smain.turn = (unsigned int)smain.turn == X ? (unsigned char)O : (unsigned char)X;
  }

  *winX = (unsigned int)smain.b[0] == X && ((unsigned int)smain.b[1] == X && (unsigned int)smain.b[2] == X || (unsigned int)smain.b[3] == X && (unsigned int)smain.b[6] == X) || (unsigned int)smain.b[8] == X && ((unsigned int)smain.b[7] == X && (unsigned int)smain.b[6] == X || (unsigned int)smain.b[5] == X && (unsigned int)smain.b[2] == X) || (unsigned int)smain.b[4] == X && ((unsigned int)smain.b[0] == X && (unsigned int)smain.b[8] == X || (unsigned int)smain.b[2] == X && (unsigned int)smain.b[6] == X || (unsigned int)smain.b[1] == X && (unsigned int)smain.b[7] == X || (unsigned int)smain.b[3] == X && (unsigned int)smain.b[5] == X);

  *winO = (unsigned int)smain.b[0] == O && ((unsigned int)smain.b[1] == O && (unsigned int)smain.b[2] == O || (unsigned int)smain.b[3] == O && (unsigned int)smain.b[6] == O) || (unsigned int)smain.b[8] == O && ((unsigned int)smain.b[7] == O && (unsigned int)smain.b[6] == O || (unsigned int)smain.b[5] == O && (unsigned int)smain.b[2] == O) || (unsigned int)smain.b[4] == O && ((unsigned int)smain.b[0] == O && (unsigned int)smain.b[8] == O || (unsigned int)smain.b[2] == O && (unsigned int)smain.b[6] == O || (unsigned int)smain.b[1] == O && (unsigned int)smain.b[7] == O || (unsigned int)smain.b[3] == O && (unsigned int)smain.b[5] == O);

  *finished = *winX || *winO || (unsigned int)smain.b[0] != EMPTY && (unsigned int)smain.b[1] != EMPTY && (unsigned int)smain.b[2] != EMPTY && (unsigned int)smain.b[3] != EMPTY && (unsigned int)smain.b[4] != EMPTY && (unsigned int)smain.b[5] != EMPTY && (unsigned int)smain.b[6] != EMPTY && (unsigned int)smain.b[7] != EMPTY && (unsigned int)smain.b[8] != EMPTY;

}

void main() {
  _Bool clock;
  unsigned char imove;
  _Bool winX;
  _Bool winO;
  _Bool finished;
  initial();
  //while(1) {
    imove = nondet_uchar();
    tictactoe(clock, imove, &winX, &winO, &finished);

   /*
   assert((((smain.move >> 3)&1)==0) || ((smain.move&3)==0));
   assert(!(winX==1 && winO==1));
   assert(!(smain.b[0]==X && smain.b[1]==X && smain.b[2]==X) || (winX==1)); // top row
   assert(!(smain.b[0]==O && smain.b[1]==O && smain.b[2]==O) || (winO==1));
   assert(!(smain.b[3]==X && smain.b[4]==X && smain.b[5]==X) || (winX==1)); // middle row
   assert(!(smain.b[3]==O && smain.b[4]==O && smain.b[5]==O) || (winO==1));
   assert(!(smain.b[6]==X && smain.b[7]==X && smain.b[8]==X) || (winX==1)); // bottom row
   */

   assert(!(smain.b[6]==O && smain.b[7]==O && smain.b[8]==O) || (winO==1));
  //}
}
