       Subroutine UNIX_EXTEND( o_Name, n_Name, Flag )
C
C       Looks through the file name as given, and expands any
C       environment name (e.g., $ARKINFO) into full pathname.
C
       Character*( * ) o_Name, n_Name
       Integer Flag
C
       Integer L, M, o_Len, n_Len, Colon, Tilde, Slash, Dollar
       Character One, Work*7
C
        Flag = 1
       o_Len = Len( o_Name )
       n_Len = Len( n_Name )
       M = 0
       Colon = 0
       Tilde = 0
       Slash = 0
        Dollar = 0
       Do L = 1, o_Len
         One = o_Name( L: L )
         If( One .eq. ':' ) Then
           Colon = L
           M = L
         Else If( One .eq. '~' ) Then
           Tilde = L
           M = L
         Else If( One .eq. '/' ) Then
           If( Slash .eq. 0 .and.
     1             ( Colon.gt.0 .or. Tilde.gt.0 .or. Dollar.gt.0) ) Then
             Slash = L
           End If
           M = L
          Else If( One .eq. '$' ) Then
            Dollar = L
            M = L
         Else If( One .gt. ' ' ) Then
           M = L
         End If
       End Do
       o_Len = M
C
       If( Colon.gt.0 .or. Tilde.gt.0 .or. Dollar.gt.0) Then
         If(Colon.eq.8 .or. Colon.eq.7) Then
            Work=o_name(1:Colon-1)
            call upcase(Work)
            If (Work.eq.'ARKTEST' .or.
     &          Work.eq.'ARKDOCS' .or.
     &          Work.eq.'ARKDATA' .or.
     &          Work.eq.'FOSDAT') Then
c             Env name found, which matches an ARK one.
             Call ARKGSM( Work, n_Name )
             Slash = Colon
            End if
         Else If( Tilde .gt. 0 ) Then
c                                   ! ~ used
           If( Slash .ne. Tilde + 1 ) Then
             Flag = -11
             Write( *, 120 )
 120         Format( ' ERROR:: Cannot expand other user''s home' )
             Return
           Else If( Tilde .ne. 1 ) Then
             Flag = 0
             Write( *, 125 )
125             Format( ' ERROR:: UNIX file name syntax error' )
           End If
           Call ARKGSM( 'HOME', n_Name )
          Else If( Dollar .gt. 0) Then
            Call ARKGSM( o_Name( Dollar + 1: Slash - 1 ), n_Name )
         End If
         M = 0
         Do L = 1, n_Len
           If( n_Name( L: L ) .gt. ' ' ) M = L
         End Do
         If( M .eq. 0 ) Then
           Flag = -12
           Write(*,*) ' ERROR:: Non-existent environment name ', work
           Return
         Else If( M + o_Len - Slash + 1 .gt. n_Len ) Then
           Flag = -13
           Write( *, 140 )
140           Format( ' ERROR:: Running out of file name buffer' )
           Return
         End If
         n_Name = n_Name( 1: M ) //'/'//o_Name( Slash+1: o_Len )
       Else
         If( o_Len .gt. n_Len ) Then
           Flag = -13
           Write( *, 140 )
           Return
         End If
         n_Name = o_Name
       End If
C
       End
    
