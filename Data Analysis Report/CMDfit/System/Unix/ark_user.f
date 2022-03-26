       subroutine ARK_USER( usrnam )
C
C       Tries to guess user's login name from his/her home directory
C
       character*( * ) usrnam
C
       character*64 hname
       integer l, slash, last
C
       call GETENV( 'HOME', hname )
       slash = 0
       last = 0
       do l = 1, 64
         if( hname( l: l ) .gt. ' ' ) then
           last = l
           if( hname( l: l ) .eq. '/' ) then
             slash = l
           end if
         end if
       end do
       if( slash .eq. 0 ) then
         usrnam = hname( 1: last )
       else
         usrnam = hname( slash + 1: last )
       end if
C
       end
