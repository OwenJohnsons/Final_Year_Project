      Subroutine SETDT( Type )

      use get_header
C
       Integer Type
C
       Integer Type2, K, n_Dross
C
       Character*80 Header( -1: 6 )
       Data Header/'DATATYPE=                    Z /DIMENSIONLESS',
     1       'DATATYPE=                    U /UNKNOWN',
     2       'DATATYPE=                    C /COUNTS',
     3       'DATATYPE=                    E /ERRORS (COUNTS)',
     4       'DATATYPE=                    L /F-LAMBDA',
     5       'DATATYPE=                    F /ERROR (F-LAMBDA)',
     6       'DATATYPE=                    N /F-NU',
     7       'DATATYPE=                    G /ERROR (F-NU)' /
C                            Find where to put the datatype
       Type2 = Type
       If( Type2 .lt. -1 .or. Type2 .gt. 6 ) Type2 = 0
       n_Dross = LNDRSS( )
       Do K = 1, n_Dross
         If( Dross( K )( 1: 8 ) .eq. 'DATATYPE' ) Goto 100
       End Do
       If( Dross( n_Dross ) .eq. 'END' ) Then
         If( n_Dross .le. 35 ) Then
           Dross( n_Dross + 1 ) = 'END'
           K = n_Dross
         Else
           K = 0
         End If
       Else
         If( n_Dross .le. 34 ) Then
           Dross( n_Dross + 2 ) = 'END'
           K = n_Dross + 1
         Else If( n_Dross .eq. 35 ) Then
           K = n_Dross + 1
         Else
           K = 0
         End If
       End If
C                            Insert the datatype header
100       Continue
       If( K .ne. 0 ) Then
         Dross( K ) = Header( Type2 )
       End If
C
       End
