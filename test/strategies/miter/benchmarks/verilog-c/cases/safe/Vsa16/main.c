#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

// Control states.
int      IF = 0;                // instruction fetch
int  ID = 1;            // instruction decode
int  EX = 2;            // execution
int  MEM = 3;           // memory access
int  WB = 4;            // write-back

// Opcodes.
int      LW  = 0;
int  SW    = 1;
int  BEQZ  = 2;
int  ALUop = 3;
int  ADDI  = 4;
int  SUBI  = 5;

// ALU function codes.
int         ADD = 0;
int               SUB = 1;
int               AND = 2;
int               OR  = 3;
int               XOR = 4;
int               SRL = 5;

struct state_elements_alu16 {
  unsigned short int    Registers[4];   // regiter file
  unsigned short int      PC;                   // program counter
  unsigned short int      NPC;                  // next program counter
  unsigned short int      IR;                   // instruction register
  unsigned short int      A;                    // first ALU operand
  unsigned short int      B;                    // second ALU operand
  unsigned short int      ALUOutput;            // ALU output register
  _Bool           Cond;                 // comparison result
  unsigned short int      LMD;                  // load memory data register
  unsigned char   State;                // state for non-pipelined unit
};
struct state_elements_alu16 sa;

void initial_alu16()
{
        for (unsigned i = 0; i < 4; i = i + 1)
            sa.Registers[i] = 0;
        sa.PC = 0;
        sa.NPC = 0;
        sa.IR = 0;
        sa.A = 0;
        sa.B = 0;
        sa.ALUOutput = 0;
        sa.Cond = 0;
        sa.LMD = 0;
        sa.State = IF;
}

void alu16(
  _Bool clk,
  unsigned short int *PC,
  unsigned short int instruction,
  unsigned short int *ALUOutput,
  unsigned short int datain,
  unsigned short int *dataout,
  _Bool *wr)
{
    // Interesting fields of the instruction register.
    unsigned char         opcode = sa.IR&0x7;
    unsigned char         adFld1 = ((sa.IR>>3)&0x3);
    unsigned char         adFld2 = ((sa.IR>>3)&0x3);
    unsigned char         adFld3 = ((sa.IR>>7)&0x3);
    unsigned short int immFld = ((sa.IR>>7)&0x1FF);
    unsigned char       funFld = ((sa.IR>>9)&0x7F);


  // Decoding of the instruction type.
  _Bool memRef = (opcode == LW || opcode == SW);
  _Bool regRegALU = (opcode == ALUop);
  _Bool regImmALU = (opcode == ADDI || opcode == SUBI);
  _Bool branch = (opcode == BEQZ);
  // Immediate operand with sign extension.
  unsigned short int    Imm =  ((((immFld>>8)&0x1) & 0xFF) << 7) | (immFld & 0xFF);

  // clocked block
        switch(sa.State) {
         case IF: {
              sa.NPC = sa.PC + 2;
              sa.IR = instruction;
              break;
        } // case: IF
         case ID: {
              sa.A = sa.Registers[adFld1];
              sa.B = sa.Registers[adFld2];
              break;
        } // case: ID
         case EX: {
              if (memRef) {
                *ALUOutput = sa.A + Imm;
              } // if (memRef)
              else if (regRegALU) {
                if (funFld == ADD)
                  *ALUOutput = sa.A + sa.B;
                else if (funFld == SUB)
                  *ALUOutput = sa.A - sa.B;
                else if (funFld == AND)
                  *ALUOutput = sa.A & sa.B;
                else if (funFld == OR)
                  *ALUOutput = sa.A | sa.B;
                else if (funFld == XOR)
                  *ALUOutput = sa.A ^ sa.B;
                else if (funFld == SRL)
                  *ALUOutput = ((0&0x0000)<<15) | ((sa.A>>1)&0x7FFF);
              } // if (regRegALU)
              else if (regImmALU) {
                if (opcode == ADDI)
                  *ALUOutput = sa.A + Imm;
                else if (opcode == SUBI)
                  *ALUOutput = sa.A - Imm;
              } // if (regImmALU)
              else if (branch) {
                *ALUOutput =  ((0x0<<11) | (sa.NPC&0xfff)) + ((Imm<<1) | 0);
                sa.Cond = (sa.A == 0);
              } // if (branch)
            break;
            } // case: EX
   case MEM: {
               if (memRef) {
                 if (opcode == LW)
                   sa.LMD = datain;
               } // if (memRef)
               if (branch) {
                 if (sa.Cond)
                   sa.PC = *ALUOutput&0xFFF;
                 else
                   sa.PC = sa.NPC;
               }
               else
                 sa.PC = sa.NPC;
             } // case: MEM
   case WB: {
              if (regRegALU) {
                if (adFld3 != 0)
                  sa.Registers[adFld3] = *ALUOutput;
              }
              else if (regImmALU) {
                if (adFld2 != 0)
                  sa.Registers[adFld2] = *ALUOutput;
              }
              else if (opcode == LW) {
                if (adFld2 != 0)
                  sa.Registers[adFld2] = sa.LMD;
              }
            } // case: WB
  } // case (State)

  // State update.
        if (sa.State == 4)
            sa.State = 0;
        else
            sa.State = sa.State + 1;

  // Combinational outputs.
  *dataout = sa.B;
  *wr = (sa.State == MEM) & (opcode == SW);

  // update the output
  *PC = sa.PC;

  // #PASS: The program counter is always aligned on a half-word boundary.
  //assert((sa.PC&0x1)==0 && (sa.NPC&0x1)==0);

  // #PASS: Register R0 is never written.
  assert(sa.Registers[0]==0);

  // #PASS: The state is between 0 and 4.
  //assert((((sa.State>>2)&0x1)==0) || ((sa.State&0x3)==0));
}

int main()
{
  _Bool clk;
  unsigned short int PC;
  unsigned short int instruction;
  unsigned short int ALUOutput;
  unsigned short int datain;
  unsigned short int dataout;
  _Bool wr;

  initial_alu16();
  //while(1)
  alu16(clk,&PC,instruction,&ALUOutput,datain,&dataout,&wr);
  alu16(clk,&PC,instruction,&ALUOutput,datain,&dataout,&wr);
  alu16(clk,&PC,instruction,&ALUOutput,datain,&dataout,&wr);
  alu16(clk,&PC,instruction,&ALUOutput,datain,&dataout,&wr);
  alu16(clk,&PC,instruction,&ALUOutput,datain,&dataout,&wr);
}
