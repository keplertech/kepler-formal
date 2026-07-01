// Author: Fabio Somenzi <Fabio@Colorado.EDU>

module main (clock,reset,dataIn,dataOut);
    input        clock;
    inout        reset;
    input [7:0]  dataIn;
    output [7:0] dataOut;
    reg [7:0]    dataOut;

    reg [7:0]    regIn;
    wire [7:0]   transformed;
    wire         Lcmd, Ucmd, Ncmd, Ccmd;

    initial begin
        regIn   = 0;
        dataOut = 0;
    end

    control     c(clock,reset,regIn,Lcmd,Ucmd,Ncmd,Ccmd);

    transform   t(regIn,Lcmd,Ucmd,Ncmd,Ccmd,transformed);

    always @ (posedge clock) begin
        if (reset) begin
            dataOut = 0;
            regIn = 0;
        end else begin
            dataOut = transformed;
            regIn = dataIn;
        end // else: !if(reset)
    end

//#PASS: The command lines are 1-hot encoded.
assert property (!(Lcmd==1 && Ucmd==1 || Lcmd==1 && Ccmd==1 || Lcmd==1 && Ncmd==1 || Ucmd==1 && Ccmd==1
   || Ucmd==1 && Ncmd==1 || Ccmd==1 && Ncmd==1));
/*#PASS: If the output is a lowercase character and Lcmd is asserted,
# then Lcmd must be deasserted or the output continues to be a lowercase
# character.*/

assert property ((Lcmd==1 && dataOut[5]==1) |-> ##1 (Lcmd==0)); // Fails at fourth cycle

endmodule // lunc


module control (clock,reset,in,Lcmd,Ucmd,Ncmd,Ccmd);
    input       clock;
    input       reset;
    input [7:0] in;
    output      Lcmd;
    output      Ucmd;
    output      Ncmd;
    output      Ccmd;

    reg         Lcmd;
    reg         Ucmd;
    reg         Ncmd;
    reg         Ccmd;

    wire        load;
    reg [7:0]   prev;

    initial begin
        Lcmd = 0;
        Ucmd = 0;
        Ncmd = 1;
        Ccmd = 0;
        prev = 0;
    end // initial begin

    always @ (posedge clock)
      if (reset)
        prev = 0;
      else
        prev = in;

    assign load = (prev == 8'h1b); // escape

    always @ (posedge clock)
      if (reset) begin
          Ncmd = 1;
          Lcmd = 0; Ucmd = 0; Ccmd = 0;
      end
      else if (load) begin
          case (in)
            8'h4c: begin // L
                Lcmd = 1; Ucmd = 0; Ncmd = 0; Ccmd = 0;
            end
            8'h55: begin // U
                Lcmd = 0; Ucmd = 1; Ncmd = 0; Ccmd = 0;
            end
            8'h4e: begin // N
                Lcmd = 0; Ucmd = 0; Ncmd = 1; Ccmd = 0;
            end
            8'h43: begin // C
                Lcmd = 0; Ucmd = 0; Ncmd = 0; Ccmd = 1;
            end
            default: begin
                Lcmd = 1'bx; Ucmd = 1'bx; Ncmd = 1'bx; Ccmd = 1'bx;
            end
          endcase
      end // if (load)

endmodule // control


module transform (in,Lcmd,Ucmd,Ncmd,Ccmd,out);
    input [7:0]  in;
    input        Lcmd;
    input        Ucmd;
    input        Ncmd;
    input        Ccmd;
    output [7:0] out;

    assign out = Lcmd ? toLower(in) :
           Ucmd ? toUpper(in) :
           Ncmd ? in :
           Ccmd ? changeCase(in) : 8'hxx;

    function [7:0] toLower;
        input [7:0]   in;
    begin: _toLower
        if (isUpper(in))
            toLower = in + 8'h20;
        else
            toLower = in;
    end
    endfunction // toLower

    function [7:0] toUpper;
        input [7:0]   in;
    begin: _toUpper
        if (!isUpper(in))
            toUpper = in - 8'h20;
        else
            toUpper = in;
    end
    endfunction // toUpper

    function [7:0] changeCase;
        input [7:0]   in;
    begin: _changeCase
        if (isUpper(in))
            changeCase = in + 8'h20;
        else
            changeCase = in - 8'h20;
    end
    endfunction // changeCase

    function isUpper;
        input [7:0]   in;
    begin: _isUpper
        isUpper = ~in[5];
    end
    endfunction // isUpper

endmodule // transform
