#include <stdio.h>
#define TRUE 1
#define FALSE 0

#define IDLE 0
#define READY 1
#define BUSY 2

#define NO_REQ 0
#define REQ 1
#define HAVE_TOKEN 2


int A = 0;
int B = 1;
int C = 2;
int X = 3;

struct state_elements_controller{
  _Bool ack;
  unsigned char state;
  _Bool pass_token;
};
struct state_elements_controller  scontroller;

void controller_initial()
{
  scontroller.state = IDLE;
  scontroller.ack = 0;
  scontroller.pass_token = 1;
}

void controller(_Bool clk, _Bool req, _Bool *ack, unsigned char sel, _Bool *pass_token, unsigned char id)
{
  _Bool is_selected;
  is_selected = (sel&0x3 == id&0x3);

  // clocked block
  if(scontroller.state == IDLE) {
    if(is_selected)
      if(req)
      {
        scontroller.state = READY;
        scontroller.pass_token = FALSE;
      }
      else
        scontroller.pass_token = TRUE;
    else
      scontroller.pass_token = FALSE;
  }
  else if(scontroller.state == READY)
  {
    scontroller.state = BUSY;
    scontroller.ack = TRUE;
  }
  else if(scontroller.state == BUSY) {
    if(!req)
    {
      scontroller.state = IDLE;
      scontroller.ack = FALSE;
      scontroller.pass_token = TRUE;
    }
  }
  // pass the output variable
  *ack = scontroller.ack;
  *pass_token = scontroller.pass_token;

}

struct state_elements_arbiter{
  unsigned char state;
};
struct state_elements_arbiter sarbiter;

void arbiter_initial()
{
  sarbiter.state = A;
}

void arbiter(_Bool clk, unsigned char *sel, _Bool active)
{
  if(active) {
    if(sarbiter.state == A)
      sarbiter.state = B;
    else if(sarbiter.state == B)
      sarbiter.state = C;
    if(sarbiter.state == C)
      sarbiter.state = A;
  }
  // continuous assignment
  *sel = active ? sarbiter.state: X;
}

struct state_elements_client{
  unsigned char state;
  _Bool req;
};
struct state_elements_client sclient;

void client_initial()
{
  sclient.state = NO_REQ;
  sclient.req = 0;
}
void client(_Bool clk, _Bool *req, _Bool ack)
{
  _Bool rand_choice = 0;
  if(sclient.state == NO_REQ) {
    if (rand_choice) {
      sclient.req = 1;
      sclient.state = REQ;
    }
  }
  else if(sclient.state == REQ) {
    if (ack)
      sclient.state = HAVE_TOKEN;
  }
  else if(sclient.state == HAVE_TOKEN)
  {
    if (rand_choice) {
      sclient.req = 0;
      sclient.state = NO_REQ;
    }
  }
  // pass the output variable
  *req = sclient.req;
}

void main()
{
  _Bool clk, ackA, ackB, ackC;
  unsigned char sel;
  _Bool active;
  _Bool reqA, reqB, reqC;
  _Bool pass_tokenA, pass_tokenB, pass_tokenC;

  // initial functions
  client_initial();
  controller_initial();
  arbiter_initial();

  //while(1) {
    client(clk, &reqA, ackA);
    client(clk, &reqB, ackB);
    client(clk, &reqC, ackC);

    controller(clk, reqA, &ackA, sel, &pass_tokenA, A);
    controller(clk, reqB, &ackB, sel, &pass_tokenB, B);
    controller(clk, reqC, &ackC, sel, &pass_tokenC, C);

    arbiter(clk, &sel, active);

    active = pass_tokenA || pass_tokenB || pass_tokenC;

    assert(!(ackA == 1 && ackB == 1 || ackB == 1 && ackC == 1 || ackC == 1 && ackA ==1));
  //}
}
