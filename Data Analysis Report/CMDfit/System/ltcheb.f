      Subroutine LTCHEB( Y, Pix, Coef, Order, RMS )
      
       ! New version of ltcheb, where order on input is the maximum
       ! order to be fitted.

       Integer, Parameter :: maxM=9

       Integer, intent(in) :: Pix
       Real, intent(in), dimension(pix) :: y
       Double Precision, intent(out), dimension( 0: maxM ) :: Coef
       integer, intent(inout) :: order
       real, intent(out) :: RMS       
       
       Double Precision Add, Div, Mul
       Double Precision Gamma( 0: maxM ), Bj( 0: maxM )
       Double Precision Ssq1( 0: maxM ), E( 0: 1 )
       Double Precision Xi, Sx, Pw, BetaJ, Sum, Diff
       Double Precision Alpha( 0: maxM, 0: maxM ), Temp
       Real Test, ftRat( maxM )
       Integer Limit, I, J, J1, J2, K, Skip
C
       Double Precision LPPOLY
C
       Double Precision C( 0: maxM, -1: maxM )
c       ...C( I, J ) is the factor of X**I in Jth Tchebyceff polynomial
       Common / cheb / C
C
       Real Ftest( 9 )
       Data Ftest / 161.0, 18.5, 10.1, 7.71, 6.61,
     1                                          5.99, 5.59, 5.32, 5.12 /
C
       Limit = Min( maxM, Pix - 2, order)
       If( Limit .le. 0 ) Then
         Write( *, '('' Number of pixels is too small: '',I5)' ) Pix
         Return
       End If
       E( 1 ) = 1.0
       C( 0, -1 ) = 0.0
       C( 0, 0 ) = 1.0
       Mul = Real( Pix ) - 1.0
       If( Mul .eq. 0.0 ) Mul = 1.0
       Div = 1.0 / Mul
       Add = Mul - 1.0
       Skip = Pix / 100
       If( Skip .le. 0 ) Skip = 1
C
c       Section 1.  Calculate the polynomials for this data
c
       Do J = 0, Limit
         Gamma( J ) = 0.0
         Sx = 0.0
         Do I = 1, Pix, Skip
           Xi = ( I + Add ) * Div
           Pw = LPPOLY( J, Xi ) ** 2
           Gamma( J ) = Gamma( J ) + Pw
           Sx = Sx + Pw * Xi
         End Do
         If( J .le. 0 ) Then
           BetaJ = 0.0
         Else
           BetaJ = Gamma( J ) / Gamma( J - 1 )
         End If
         E( 0 ) = -Sx / Gamma( J )
         J1 = J - 1
         J2 = J + 1
         If( J2 .le. maxM ) Then
           Call LPRDCT( C( 0, J ), J, E, 1, C( 0, J2 ) )
           Do K = 0, J1
             C( K, J2 ) = C( K, J2 ) - C( K, J1 ) * BetaJ
           End Do
         End If
       End Do
c
c       Section 2.  Calculate the coefficients Bj
c
       Do J = 0, Limit
         Sum = 0.0
         Do I = 1, Pix, Skip
           Xi = ( I + Add ) * Div
           Sum = Sum + Y( I ) * LPPOLY( J, Xi )
         End Do
         Bj( J ) = Sum / Gamma( J )
       End Do
c
c       Section 3.  Calculate residuals and test for best order fit
c
       Do J = 0, Limit
         Ssq1( J ) = 0.0
       End Do
       Do I = 1, Pix, Skip
         Xi = ( I + Add ) * Div
         Sum = 0.0
         Do J = 0, Limit
           Sum = Sum + Bj( J ) * LPPOLY( J, Xi )
           Diff = Y( I ) - Sum
           Ssq1( J ) = Ssq1( J ) + Diff * Diff
         End Do
       End Do
       RMS = RMS * RMS
       Do J = 0, Limit
         K = Pix / Skip - J - 1
         Ssq1( J ) = Ssq1( J ) / K
         If( Ssq1( J ) .le. RMS ) Then
           Order = J
           Goto 350
         End If
         If( J .ge. 1 ) Then
           Test = ( Ssq1( J - 1 ) - Ssq1( J ) ) / Ssq1( J ) * K
           If( K .ge. 10 ) Then
             ftRat( J ) = Test / 5.0
           Else
             ftRat( J ) = Test / Ftest( K )
           End If
         End If
       End Do
       Do J = 1, Limit
         If( ftRat( J ) .ge. 1.0 ) Order = J
       End Do
       If( Order .ne. Limit ) Then
         Write( *, 337 ) Order
337       Format
     1    ( ' @ The fit does not improve with orders higher than ', I1 )
       Else
         Write( *, '('' * Failed to find a good fit !!!'')' )
       End If
350     Continue
       RMS = Sqrt( RMS )
c
c       Section 5. Calculate coefficients Coeff and associated errors
c
       Do J = 0, Order
         Coef( J ) = 0.0
         Do K = J, Order
           Coef( J ) = Coef( J ) + Bj( K ) * C( J, K )
         End Do
       End Do
       Do J = Order + 1, Limit
         Coef( J ) = 0.0
       End Do
C
C       convert back to the real data set
C
       Call LCNVRT( Alpha, maxM, -Add, Mul )
       Do J = 0, Order
         Temp = 0.0
         Do K = J, Order
           Temp = Temp + Alpha( K, J ) * Coef( K )
         End Do
         Coef( J ) = Temp
       End Do
C
C     END
C
       End
C
C
C
       Subroutine LPRDCT( A, M, B, N, P )
c
c       Makes the product of two polynomials
c       A( J ) * X ** J ( J = 0, M ) and B( K ) * X ** K ( K = 0, N )
c       and returns the result P( L ) * X ** L ( L = 0, M + N )
c                                       Koji, 10.1.84
C
       Integer M, N
       Double Precision A( 0:M ), B( 0:N ), P( 0:M+N )
C
       Integer J, K, L
       Integer MinA, MaxA
C
       Do L = 0, M + N
         P( L ) = 0.0
         MinA = Max( 0, L - N )
         MaxA = Min( L, M )
         Do J = MinA, MaxA
           K = L - J
           P( L ) = P( L ) + A( J ) * B( K )
         End Do
       End Do
C
       End
C
C
C
       Double Precision Function LPPOLY( J, X )
c
c       Tchebycheff polynomials with coefficient C
C
       Integer MaxM
       Parameter( MaxM = 9 )
C
       Integer J
       Double Precision X
C
       Integer K
       Double Precision Dummy
C
       Double Precision C( 0: maxM, -1: maxM )
       Common / cheb / C
C
       If( J .le. -1 ) Then
         LPPOLY = 0.0
       Else If( J .eq. 0 ) Then
         LPPOLY = 1.0
       Else
         If( X .eq. 0.0 ) Then
           LPPOLY = C( 0, J )
         Else
           Dummy = 0.0
           Do K = J, 0, -1
             Dummy = Dummy * X + C( K, J )
           End Do
           LPPOLY = Dummy
         End If
       End If
C
       End
C
C
C
       Subroutine LCNVRT( Alpha, Order, A0, A1 )
C
       Integer Order
       Double Precision Alpha( 0:Order, 0:Order )
       Double Precision A0, A1
c
       Double Precision B0, B1
       Integer I, J
C
       Do I = 0, Order
         Do J = 0, Order
           Alpha( I, J ) = 0.0
         End Do
       End Do
       B1 = 1.0 / A1
       B0 = -A0 * B1
       Alpha( 0, 0 ) = 1.0
       Do I = 1, Order
         Do J = 0, I - 1
           Alpha( I, J ) = Alpha( I, J ) + B0 * Alpha( I - 1, J )
           Alpha( I, J + 1 ) = Alpha( I, J + 1 )
     1                                          + B1 * Alpha( I - 1, J )
         End Do
       End Do
C
       End
