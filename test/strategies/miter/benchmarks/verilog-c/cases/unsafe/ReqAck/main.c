#include <stdio.h>
#define TRUE 1
#define FALSE 0

int idle = 0;
int starting = 1;
int     working = 2;
int     done = 3;

// This module is a timer.
struct state_elements_slave {
  unsigned short int count;
};
struct state_elements_slave sslave;

void initial_slave()
{
  sslave.count = 0;
}

void slave(_Bool clock, _Bool start, _Bool *ready)
{
        if (start)
          sslave.count = 0;
        else
          sslave.count = sslave.count + 1;

  *ready = (sslave.count == 0x3FF);
}

struct state_elements_reqAck {
  unsigned char state;
};
struct state_elements_reqAck sreqAck;

void initial_reqAck()
{
  sreqAck.state = idle;
}

void reqAck(_Bool clock, _Bool req, _Bool *ack)
{
  _Bool start;
  _Bool ready;

  slave (clock,start,&ready);

  if(sreqAck.state == idle)
  {
    if(req)
      sreqAck.state = starting;
    else
      sreqAck.state = idle;
  }
  if(sreqAck.state == idle)
  {
    sreqAck.state = working;
  }
  if(sreqAck.state == working)
  {
    if (ready)
      sreqAck.state = done;
    else
      sreqAck.state = working;
  }
  if(sreqAck.state == done)
  {
    sreqAck.state = idle;
  }

  ack = (sreqAck.state == done);
  start = (sreqAck.state == starting);
}

struct state_elements_main {
  _Bool req;
};
struct state_elements_main smain;

void initial_main()
{
  smain.req =  0;
}

void main()
{
  _Bool clock, ack, nd;

  nd=1;
  smain.req=nd;
  if(smain.req==1) {
    reqAck(clock, smain.req, &ack);
    assert(ack==1);
  }
}
