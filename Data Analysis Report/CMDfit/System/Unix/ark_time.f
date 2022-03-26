       subroutine ARK_TIME( hh, mm, ss )
C
C       Get time of day as 3 integers
c
       integer hh, mm, ss
c
       integer array( 3 )
c
       call ITIME( array )
       hh = array( 1 )
       mm = array( 2 )
       ss = array( 3 )
c
       end
