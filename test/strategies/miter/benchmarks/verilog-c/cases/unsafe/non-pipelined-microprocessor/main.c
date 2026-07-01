#include <stdio.h>
#define TRUE 1
#define FALSE 0

#define ADD 0
#define JMP 1
#define AND 2
#define XOR 3

#define INSTR 11
#define OPCODE 2
#define DATA_ADDR 2
#define PROG_ADDR 2
#define DATA 4


struct state_elements_memory{
unsigned char m2;
unsigned char m3;
unsigned char m4;
unsigned char m5;
unsigned char m6;
unsigned char m7;
unsigned char m0;
unsigned char m1;
};
struct state_elements_memory  smemory;

void memory(_Bool clk, unsigned char opcode, unsigned char readLoc1, unsigned char readLoc2, unsigned char writeLoc, unsigned char data, unsigned char *memOut1, unsigned char *memOut2)
{
  unsigned char m0_old;
  unsigned char m1_old;
  unsigned char m2_old;
  unsigned char m3_old;
  unsigned char m4_old;
  unsigned char m5_old;
  unsigned char m6_old;
  unsigned char m7_old;
  if((unsigned int)opcode != 1)
  {
    if((unsigned int)writeLoc == 0)
    {
      smemory.m0 = data;
    }

    if((unsigned int)writeLoc == 1)
    {
      smemory.m1 = data;
    }

    if((unsigned int)writeLoc == 2)
    {
      smemory.m2 = data;
    }

    if((unsigned int)writeLoc == 3)
    {
      smemory.m3 = data;
    }

    if((unsigned int)writeLoc == 4)
    {
      smemory.m4 = data;
    }

    if((unsigned int)writeLoc == 5)
    {
      smemory.m5 = data;
    }

    if((unsigned int)writeLoc == 6)
    {
      smemory.m6 = data;
    }

    if((unsigned int)writeLoc == 7)
    {
      smemory.m7 = data;
    }

  }

  *memOut1 = (readLoc1 == 0) ? smemory.m0 :
    (readLoc1 == 1) ? smemory.m1 :
    (readLoc1 == 2) ? smemory.m2 :
    (readLoc1 == 3) ? smemory.m3 :
    (readLoc1 == 4) ? smemory.m5 :
    (readLoc1 == 5) ? smemory.m5 :
    (readLoc1 == 6) ? smemory.m6 : smemory.m7;

  *memOut2 = (readLoc2 == 0) ? smemory.m0 :
    (readLoc2 == 1) ? smemory.m1 :
    (readLoc2 == 2) ? smemory.m2 :
    (readLoc2 == 3) ? smemory.m3 :
    (readLoc2 == 4) ? smemory.m5 :
    (readLoc2 == 5) ? smemory.m5 :
    (readLoc2 == 6) ? smemory.m6 : smemory.m7;
}

void program(_Bool clk, unsigned char progCntr, unsigned short int *instruction)
{
  unsigned short int instr0;
  unsigned short int instr1;
  unsigned short int instr2;
  unsigned short int instr3;
  unsigned short int instr4;
  unsigned short int instr5;
  unsigned short int instr6;
  unsigned short int instr7;
  instr0 = 576;
  instr1 = 1152;
  instr2 = 1728;
  instr3 = 2304;
  instr4 = 505;
  instr5 = 0;
  instr6 = 0;
  instr7 = 0;

  *instruction =  (progCntr == 0) ? instr0 :
    (progCntr == 1) ? instr1 :
    (progCntr == 2) ? instr2 :
    (progCntr == 3) ? instr3 :
    (progCntr == 4) ? instr4 :
    (progCntr == 5) ? instr5 :
    (progCntr == 6) ? instr6 : instr7;
}

void decodeOpcd(_Bool clk, unsigned short int instruction, unsigned char *opcode)
{
  *opcode = instruction & 7;
}

void decodeLoc1(_Bool clk, unsigned short int instruction, unsigned char *readLoc1)
{
  *readLoc1 = ((instruction >> 3) & 7);
}

void decodeLoc2(_Bool clk, unsigned short int instruction, unsigned char *readLoc2)
{
  *readLoc2 = ((instruction >> 6) & 7);
}

void decodeLoc3(_Bool clk, unsigned short int instruction, unsigned char *writeLoc)
{
  *writeLoc = ((instruction >> 9) & 7);
}

void alu(_Bool clk, unsigned char opcode, unsigned char operand1, unsigned char operand2, unsigned char *aluOut)
{
  *aluOut = (opcode == 0) ? (operand1 + operand2) :
    (opcode == 3) ? (operand1^operand2) :
    (opcode == 2) ? (operand1 & operand2) : 0;
}

struct state_elements_pc{
  unsigned char progCntr;
};
struct state_elements_pc  spc;

void pc(_Bool clk, unsigned char opcode, unsigned char operand1, unsigned char operand2, unsigned char *progCntr)
{
  unsigned char progCntr_old;
  if((unsigned int)opcode == 1 && (unsigned int)operand1 == 0)
  {
    spc.progCntr = operand2 & 7;
    *progCntr = operand2 & 7;
  }
  else
  {
    spc.progCntr = spc.progCntr + 1;
    *progCntr = *progCntr + 1;
  }
}

int main()
{
  _Bool clk;
  unsigned char aluOut;
  unsigned char memOut1;
  unsigned char memOut2;
  unsigned char readLoc1;
  unsigned char readLoc2;
  unsigned char writeLoc;
  unsigned char progCntr;
  unsigned short int instruction;
  unsigned char opcode;
  //while(1) {
   memory(clk, opcode, readLoc1, readLoc2, writeLoc, aluOut, &memOut1, &memOut2);
   program(clk, progCntr, &instruction);
   decodeOpcd(clk, instruction, &opcode);
   decodeLoc1(clk, instruction, &readLoc1);
   decodeLoc2(clk, instruction, &readLoc2);
   decodeLoc3(clk, instruction, &writeLoc);
   alu(clk, opcode, memOut1, memOut2, &aluOut);
   pc(clk, opcode, memOut1, memOut2, &progCntr);
   assert(!((opcode == JMP) && (memOut1 == 0)) && (progCntr == (memOut2 & 7)));
 //}
}
