#define TRUE 1
#define FALSE 0

struct state_elements_Traffic_Light{
  unsigned char Light_Sign;
  unsigned int Counter;
};
struct state_elements_Traffic_Light  sTraffic_Light;

 // parameters
  int RED_LIGHT = 0;
  int GREEN_LIGHT = 1;
  int YELLOW_LIGHT = 2;

  unsigned int RED_count = 60;
  unsigned int GREEN_count = 40;
  unsigned int YELLOW_count = 5;

void Traffic_Light(_Bool reset, _Bool clk, unsigned int *time_left)
{
  unsigned char Light_Sign_old;
  unsigned int Counter_old;

  // assignment statements
  Light_Sign_old = sTraffic_Light.Light_Sign;
  Counter_old = sTraffic_Light.Counter;

  if(!reset)
  {
    sTraffic_Light.Light_Sign = RED_LIGHT;
    sTraffic_Light.Counter = 0;
  }

  else
  {
    if(Light_Sign_old == RED_LIGHT)
        sTraffic_Light.Light_Sign = (Counter_old == 0) ? GREEN_LIGHT : RED_LIGHT;
    else
     if(Light_Sign_old == GREEN_LIGHT)
        sTraffic_Light.Light_Sign = (Counter_old == 0) ? YELLOW_LIGHT : GREEN_LIGHT;
    else
     if(Light_Sign_old == YELLOW_LIGHT)
        sTraffic_Light.Light_Sign = (Counter_old == 0) ? RED_LIGHT : YELLOW_LIGHT;

    if(Light_Sign_old == RED_LIGHT)
       sTraffic_Light.Counter = (Counter_old == 0) ? GREEN_count : Counter_old - 1;
    else
      if(Light_Sign_old == GREEN_LIGHT)
       sTraffic_Light.Counter = (Counter_old == 0) ? YELLOW_count : Counter_old - 1;

    else
      if(Light_Sign_old == YELLOW_LIGHT)
       sTraffic_Light.Counter = (Counter_old == 0) ? RED_count : Counter_old - 1;
  }

  *time_left = sTraffic_Light.Counter;
}

int main() {
  _Bool reset;
  _Bool clk;
  unsigned int time_left;
  unsigned char temp;
  // do reset
  Traffic_Light(0, clk, &time_left);

  //while(1) {
    Traffic_Light(1, clk, &time_left);
    temp = sTraffic_Light.Light_Sign;
    //assert((time_left != 0xffffffff));
    // This property failes after 42 clock cycles when GREEN_count reaches 40 down to 0
    //assert((sTraffic_Light.Light_Sign != 2));
    assert((sTraffic_Light.Light_Sign != 3));
  //}
  return 1;
}
