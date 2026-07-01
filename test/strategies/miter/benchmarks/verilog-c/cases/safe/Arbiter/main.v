//dified by: Rajdeep Mukherjee <rajdeep.mukherjee@cs.ox.ac.ik>


//typedef enum {A, B, C, X} selection;
//typedef enum {IDLE, READY, BUSY} controller_state;
//typedef enum {NO_REQ, REQ, HAVE_TOKEN} client_state;

module main(clk, ackA, ackB, ackC);
input clk;
output ackA, ackB, ackC;

parameter A =0;
parameter B =1;
parameter C =2;
parameter X =3;

parameter IDLE = 0;
parameter READY = 1;
parameter BUSY = 2;

parameter NO_REQ = 0;
parameter REQ = 1;
parameter HAVE_TOKEN = 2;

wire [1:0] sel;
wire active;

assign active = pass_tokenA || pass_tokenB || pass_tokenC;

controller controllerA(clk, reqA, ackA, sel, pass_tokenA, A);
controller controllerB(clk, reqB, ackB, sel, pass_tokenB, B);
controller controllerC(clk, reqC, ackC, sel, pass_tokenC, C);
arbiter arbiter(clk, sel, active);

client clientA(clk, reqA, ackA);
client clientB(clk, reqB, ackB);
client clientC(clk, reqC, ackC);

assert property1: !(ackA == 1 && ackB == 1 || ackB == 1 && ackC == 1 || ackC == 1 && ackA ==1);
endmodule

module controller(clk, req, ack, sel, pass_token, id);
input clk, req, sel, id;
output ack, pass_token;

parameter A =0;
parameter B =1;
parameter C =2;
parameter X =3;

parameter IDLE = 0;
parameter READY = 1;
parameter BUSY = 2;

parameter NO_REQ = 0;
parameter REQ = 1;
parameter HAVE_TOKEN = 2;

wire [1:0] sel, id;
reg ack, pass_token;
reg [1:0] state;

initial state = IDLE;
initial ack = 0;
initial pass_token = 1;

wire is_selected;
assign is_selected = (sel == id);

always @(posedge clk) begin
  case(state)
    IDLE:
      if (is_selected)
        if (req)
          begin
          state = READY;
          pass_token = 0; /* dropping off this line causes a safety bug */
          end
        else
          pass_token = 1;
      else
        pass_token = 0;
    READY:
      begin
      state = BUSY;
      ack = 1;
      end
    BUSY:
      if (!req)
        begin
        state = IDLE;
        ack = 0;
        pass_token = 1;
        end
  endcase
end
endmodule

module arbiter(clk, sel, active);
input clk, active;
output sel;


parameter A =0;
parameter B =1;
parameter C =2;
parameter X =3;

parameter IDLE = 0;
parameter READY = 1;
parameter BUSY = 2;

parameter NO_REQ = 0;
parameter REQ = 1;
parameter HAVE_TOKEN = 2;
wire [1:0] sel;
reg [1:0] state;

initial state = A;

assign sel = active ? state: X;

always @(posedge clk) begin
  if (active)
    case(state)
      A:
        state = B;
      B:
        state = C;
      C:
        state = A;
    endcase
end
endmodule

module client(clk, req, ack);
input clk, ack;
output req;

parameter A =0;
parameter B =1;
parameter C =2;
parameter X =3;

parameter IDLE = 0;
parameter READY = 1;
parameter BUSY = 2;

parameter NO_REQ = 0;
parameter REQ = 1;
parameter HAVE_TOKEN = 2;
reg req;
reg [1:0] state;

// RM: This can be nondeterministically set to 0 or 1
wire rand_choice = 1'b0;

initial req = 0;
initial state = NO_REQ;

always @(posedge clk) begin
  case(state)
    NO_REQ:
      if (rand_choice)
        begin
        req = 1;
        state = REQ;
        end
    REQ:
      if (ack)
        state = HAVE_TOKEN;
    HAVE_TOKEN:
      if (rand_choice)
        begin
        req = 0;
        state = NO_REQ;
        end
  endcase
end
endmodule
