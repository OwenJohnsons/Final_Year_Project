       Subroutine SETBUG( )
C
C       Call this at the start of a program.  Then in any subroutines where
C       you have debug statements include in the declarations;
C
       Logical Debug
       Common / bug / Debug
C
C       Then you can write statements like;
C       if (debug) print*, '@ Variable is ', i
C       These statements will only be used if the VMS symbol debug has been
C       set to YES or yes.
C
c       Original by Timn, better version by Koji, ARK Version 3 Version by
c       Timn.
c
       logical Set, arkgsm
       Character*3 String
C
       Set = arkgsm( 'ARK_DEBUG', String )
       If( Set ) THen
         If( String .eq. 'YES' .or. String .eq. 'yes' ) Then
           Debug = .True.
         Else If( String .eq. 'NO' .or. String .eq. 'no' ) Then
           Debug = .False.
         Else
           Write( *, 100 )
100        Format( ' ERROR:: Illegal value of DCL symbol ARK_DEBUG; ',
     1             'use ''yes'' or ''no''.' )
         End If
       Else
         Debug = .False.
       End If
C
       End
