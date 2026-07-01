#include <stdio.h>
#define TRUE 1
#define FALSE 0

unsigned char nondet_uchar();

// parameters
  int BITS = 4;
  int NONE = 0;
  int NICKEL = 1;
  int DIME = 2;
  int QUARTER = 3;
  int ACCEPTING = 0;
  int CHANGE = 1;
  int REFUND = 2;
  int BEVERAGE = 3;

struct state_elements_vending{
  unsigned char t5;
  unsigned char t10;
  unsigned char t25;
  unsigned char l5;
  unsigned char l10;
  unsigned char l25;
  unsigned char change;
  unsigned char state;
};
struct state_elements_vending svending;

void initial_vending(_Bool clock, unsigned char deposit, unsigned char *change, _Bool *beverage, _Bool *enable)
{
  svending.t5 = 0;
  svending.t10 = 0;
  svending.t25 = 0;
  svending.l5 = 0;
  svending.l10 = 0;
  svending.l25 = 0;
  svending.state = ACCEPTING;
  svending.change = NONE;
  *change = NONE;
}

void vending(_Bool clock, unsigned char deposit, unsigned char *change, _Bool *beverage, _Bool *enable)
{
  unsigned char total;

  total = ((0<<3) || (svending.l5&7)) + ((0<<3) || ((svending.l10 & 3)<<1) || 0) + ((0<<3) || ((svending.l25 & 1) << 2) || ((0&1)<<1) || (svending.l25 & 1));
  if(svending.state == ACCEPTING)
  {
      if((unsigned int)total >= 5)
      {
          svending.change = deposit & 3;
          *change = deposit & 3;
          svending.state = CHANGE;
      }
      else
      {
          if(deposit == NICKEL)
          {
            if(svending.t5 == 15)
            {
                svending.change = NICKEL;
                *change = NICKEL;
            }
            else
            {
                svending.change = NONE;
                *change = NONE;
              svending.t5 = (svending.t5 + 1)&15;
              svending.l5 = (svending.l5 + 1)&7;
            }
          }
          else
            if(deposit == DIME)
            {
              if(svending.t10 == 15)
              {
                  svending.change = DIME;
                  *change = DIME;
              }
              else
              {
                  svending.change = NONE;
                  *change = NONE;
                svending.t10 = (svending.t10 + 1)&15;
                svending.l10 = (svending.l10 + 1)&3;
              }
            }
            else
              if(deposit == QUARTER)
              {
                if(svending.t25 == 15)
                {
                    svending.change = QUARTER;
                    *change = QUARTER;
                }
                else
                {
                    svending.change = NONE;
                    *change = NONE;
                  svending.t25 = (svending.t25 + 1)&15;
                  svending.l25 = (svending.l25 + 1)&1;
                }
              }
              else
                if(deposit == NONE)
                {
                    svending.change = NONE;
                    *change = NONE;
                }

      } // end of else
    } // end of ACCEPTING

    else
      if(svending.state == CHANGE)
      {
        if((unsigned int)total == 5)
        {
            svending.change = NONE;
            *change = NONE;
          svending.state = BEVERAGE;
        }
        else
          if((unsigned int)total == 6)
          {
            if((unsigned int)svending.t5 > 0)
            {
                svending.change = NICKEL;
                *change = NICKEL;
              svending.t5 = (svending.t5 - 1)&15;
              svending.state = BEVERAGE;
            }
            else
            {
              {
                svending.change = NONE;
                *change = NONE;
              }
              svending.state = REFUND;
            }
          }
          else
          {
            if((unsigned int)svending.l10 > 0)
            {
                svending.change = DIME;
                *change = DIME;
              svending.t10 = (svending.t10 - 1)&15;
              svending.l10 = (svending.l10 - 1)&3;
            }
            else
            {
                svending.change = NICKEL;
                *change = NICKEL;
              svending.t5 = (svending.t5 - 1)&15;
              svending.l5 = (svending.l5 - 1)&7;
            }
          }
      }
      else
        if(svending.state == BEVERAGE)
        {
            svending.change = NONE;
            *change = NONE;
          svending.l5 = 0;
          svending.l10 = 0;
          svending.l25 = 0;
          svending.state = ACCEPTING;
        }
        else
          if(svending.state == REFUND)
          {
            if((unsigned int)svending.l10 > 0)
            {
              svending.l10 = (svending.l10 - 1)&3;
              svending.t10 = (svending.t10 - 1)&15;
                svending.change = DIME;
                *change = DIME;
            }
            else
            {
              svending.state = ACCEPTING;
                svending.change = NONE;
                *change = NONE;
            }
          }


  *beverage = (svending.state == BEVERAGE);
  *enable = (svending.state == ACCEPTING && ((total&15) < 5));


 // Specifications
 // Invariant properties
 //assert((((total >> 3)&1)==0) || (((total >> 1)&3)==0));
 //assert((!(total&15==8)) || (!((svending.l5 & 7) ==0)));
 //assert((!(total&15==7)) || ( (!((svending.l5&7)==0)) || (!((svending.l10 & 3)==0))));
 //#PASS: No more than 45c can be deposited during a transaction.
 //assert((((total>>3)&1)==0) || (((total>>1)&3)==0));


 //#PASS: If in the CHANGE state the total is 35c or more, then either a
 //#      nickel or a dime has been deposited in the current transaction.
 //assert((!(svending.state==CHANGE && ((total&15)==7) || (((total>>3)&1) == 1))) || (!(((svending.l5 & 7) == 0) && ((svending.l10 & 3) ==0))));

 //#PASS: In the REFUND state we have no quarters from this transaction,
 //#      and no nickels at all.
 //assert((!(svending.state==REFUND)) || (((svending.l25 & 1)==0) && ((svending.l5 & 7)==0) && ((svending.t5 & 15)==0)));

 //#PASS: On entry to REFUND, we have exactly three dimes from this transaction.
 // assert((!svending.state==REFUND) || ##1(!v.state==REFUND || v.l10[1:0]==3));

 /*#PASS: In the BEVERAGE state we have exactly 25c from this transaction.
 #      However, total is not up to date if we borrowed a nickel from the
 #      total count to give change out of three dimes.  Hence, total may
 #      read either 25c or 30c.*/
 assert((!(svending.state==BEVERAGE)) || ((((total>>2)&3)==1) && ((total&3)==1) || ((total&3)==2)));

}

struct state_elements_monitor{
  unsigned char balance;
};
struct state_elements_monitor smonitor;

void initial_monitor(_Bool clock, unsigned char deposit, _Bool beverage, unsigned char change, unsigned char *balance)
{
  smonitor.balance = 0;
  *balance = 0;
}

void monitor(_Bool clock, unsigned char deposit, _Bool beverage, unsigned char change, unsigned char *balance)
{
  unsigned char deposit;
  unsigned char change;
  unsigned char valD;
  unsigned char valC;
  unsigned char valB;
  valD = ((deposit == QUARTER) ? 5 : ((deposit == DIME) ? 2 : ((deposit == NICKEL) ? 1 : 0)));
  valC = ((change == QUARTER) ? 5 : (change == DIME ? 2 : (change == NICKEL ? 1 : 0)));
  valB = (beverage == 1 ? 5 : 0);
  smonitor.balance = (smonitor.balance + valD) - (valC + valB);
  *balance = (smonitor.balance + valD) - (valC + valB);
}

unsigned char deposit;

void main()
{

  _Bool clock;
  unsigned char balance;
  unsigned char change;
  unsigned char nd;
  _Bool beverage;
  _Bool enable;

  initial_vending(clock, deposit, &change, &beverage, &enable);
  initial_monitor(clock, deposit, beverage, change, &balance);

  deposit = nondet_uchar() & 3;
  vending(clock, deposit, &change, &beverage, &enable);
  monitor(clock, deposit, beverage, change, &balance);

  if(enable)
    deposit = nd & 3;
  else
    deposit = NONE;
}
