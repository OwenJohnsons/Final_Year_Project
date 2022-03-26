       Real Function RIPOLY( Argmnt, Coeff, Order )
C
       Integer Argmnt, Order
       Double Precision Coeff( 0: Order )
       Double Precision Dummy, Arg
       Integer I
C
       Arg = Argmnt
       Dummy = 0.0
       Do 100 I = Order, 0, -1
         Dummy = Dummy * Arg +  Coeff( I )
100     Continue
C
       RIPOLY = Dummy

       End

