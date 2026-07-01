#include <stdio.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0

// parameter
int RxInit_c = 0x3FF;

struct state_elements_control {
    unsigned char        rxBuf;
    unsigned char        txBuf;
    _Bool        shiftLoad;
    unsigned char freqDiv;
};
struct state_elements_control scontrol;

void initial_control()
{
  scontrol.rxBuf = 0x80;
  scontrol.txBuf = 0x01;
  scontrol.shiftLoad = 1;
  scontrol.freqDiv = 0;
}

void control(_Bool clock, _Bool reset, _Bool ld, unsigned char dataIn, _Bool *enable,
    unsigned char *parallelOut, unsigned char parallelIn, _Bool *shiftLoad, _Bool txEmpty,
    _Bool dataRdy, _Bool bitClock, _Bool *ok)
{
  *enable = (scontrol.freqDiv == 7);

  // clocked block
  if (reset == 0) {
    scontrol.shiftLoad = 1;
    scontrol.freqDiv = 0;
    if (ld == 1)
      scontrol.txBuf = dataIn;
  }
  else {
    if (dataRdy == 1) {
      scontrol.rxBuf = parallelIn;
    }
    if (*enable == 1 && txEmpty == 1) {
      if (scontrol.shiftLoad == 1) {
        scontrol.shiftLoad = 0;
      } // if (shifLoad == 1)
      else {
        scontrol.shiftLoad = 1;
      } // else: !if(shifLoad == 1)
    }    // if (txEmpty)
    else if (ld == 1) {
      scontrol.txBuf = dataIn;
    }
    scontrol.freqDiv = ((scontrol.freqDiv + 1)&0xF);
  } // else: !if(reset == 0)

  // continuous assignment statement
  *ok = (scontrol.rxBuf == scontrol.txBuf);
  *parallelOut = scontrol.txBuf;
}

struct state_elements_rx {
    unsigned short int   RxReg;         // the receive register
    unsigned char        Count16;       // to divide by 16
    _Bool        RxMT;          // receive register empty
    _Bool        RxIn;          // registered serial input
};
struct state_elements_rx srx;

void initial_rx()
{
        srx.RxReg = RxInit_c;
        srx.Count16 = 0;
        srx.RxMT = 1;
        srx.RxIn = 0;
} // initial begin

void UartRx(_Bool clk, _Bool ResetF, _Bool Serial_InT, _Bool *DataRdyT, unsigned char *DataOuT, _Bool *BitClkT)
{
        // clocked block
        // Clock serial input into RxIn
        srx.RxIn = Serial_InT;

        // Check if a data word is received
        if (DataRdyT == 1) {
           srx.RxMT = 1;
  }

  // Reset
        if (ResetF == 0) {
            srx.Count16 = 0;            // reset divide by 16 counter
            srx.RxMT = 1;                       // new message starting
            srx.RxReg = RxInit_c;
  }
        // Start bit
        else if (srx.RxMT == 1 && srx.RxIn == 0) {
            srx.Count16 = 0;            // reset divide by 16 counter
            srx.RxMT = 0;                       // new message starting
            srx.RxReg = RxInit_c;
  }

        // If in a receive transaction mode
        // if @ mid bit clock then clock data into register
        else if (srx.Count16 == 7 && srx.RxMT == 0) {   // mid clock
            srx.RxReg = (srx.RxReg&0x00) | ((srx.RxReg>>1)&0x1FF);
            srx.RxReg = (srx.RxReg&0x0FF) | (srx.RxIn<<9);
            srx.Count16 = srx.Count16 + 1;
  }

        // Normal counter increment modulo 16
        else {
            srx.Count16 = srx.Count16 + 1;
  }

  // continuous assignment
  *DataRdyT = ((srx.RxMT == 0) && (((srx.RxReg>>9)&0x1) == 1) && ((srx.RxReg&0x1) == 0));

  *BitClkT = (srx.Count16 == 9);

  *DataOuT = ((srx.RxReg>>1)&0xFF);
} // UartRx

struct state_elements_tx {
  unsigned short int    XmitReg;        // the transmit register
  unsigned char Count;          // # of serial bits sent
};
struct state_elements_tx stx;

void initial_tx()
{
        stx.Count = 0;
        stx.XmitReg = 0x3FF;
} // initial begin

void UartXmt(_Bool Shift_LdF, _Bool ClkEnbT, _Bool Clk, unsigned char DataT, _Bool ResetF, _Bool *Serial_OuT, _Bool *XmitMT)
{
  // clocked block
  if (ResetF == 0) {
    stx.XmitReg = 0x3FF;
    stx.Count   = 9;
  }
  else if (ClkEnbT == 1 && Shift_LdF == 0 && ResetF == 1) {
    stx.XmitReg   = (stx.XmitReg&0x0FF) | (1<<9);               // stop bit(s)
    stx.XmitReg = (stx.XmitReg&0x001) | ((DataT<<1)&0xff);      // payload
    stx.XmitReg   = stx.XmitReg&0x0 | 0;                // start bit
    stx.Count        = 0;
  }
        else if (ClkEnbT == 1 && Shift_LdF == 1 && ResetF == 1) {
          stx.XmitReg = (stx.XmitReg&0x00) | ((stx.XmitReg>>1)&0x1FF);
          stx.XmitReg   = (stx.XmitReg&0x0FF) | 1<<9;
          if (stx.Count != 9)
            stx.Count = stx.Count + 1;
  }

  *Serial_OuT = stx.XmitReg&0x1;
  *XmitMT = (stx.Count == 9);

} // UartXmt

struct state_elements_nullModem{
  _Bool rst;
  _Bool ld;
};
struct state_elements_nullModem  snullModem;

void initial_main()
{
        snullModem.rst = 1;
        snullModem.ld = 0;
}

void main()
{
  _Bool clock;
  _Bool reset;
  _Bool load;
  unsigned char dataIn;
  _Bool ok;
  unsigned char parallelOut;
  unsigned char parallelIn;
  _Bool shiftLoad;
  _Bool enable;
  _Bool serialOut;
  _Bool txEmpty;
  _Bool serialIn;
  _Bool dataRdy;
  _Bool bitClock;
  _Bool rst_old;
  _Bool ld_old;

  serialIn = serialOut;

  // clocked block
  snullModem.rst = reset;
  snullModem.ld = load;

  initial_main();
  initial_control();
  initial_tx();
  initial_rx();

  control(clock,snullModem.rst,snullModem.ld,dataIn,&enable,&parallelOut,parallelIn,shiftLoad,txEmpty,dataRdy,bitClock,&ok);
  UartXmt(shiftLoad,enable,clock,parallelOut,snullModem.rst,&serialOut,&txEmpty);
  UartRx(clock,snullModem.rst,serialIn,&dataRdy,&parallelIn,&bitClock);

  assert(snullModem.rst == 0 || snullModem.ld == 1);
}
