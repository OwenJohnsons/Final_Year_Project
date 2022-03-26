       Subroutine FIND_FREE( f_Unit )
C
C       Subroutine to find a free unit number
C
       Integer f_Unit
C
       Logical There
C
       f_Unit = 11
       Inquire( Unit = f_Unit, Opened = There )
       Do While( There )
         f_Unit = f_Unit + 1
         Inquire( Unit = f_Unit, Opened = There )
       End Do
C
       End
