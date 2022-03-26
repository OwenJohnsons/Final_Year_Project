       Subroutine WILDF( fildef, extdef, lastf, lstfil )
                                                  
       ! Given the file specification fildef (e.g. *.*) this routine will
       ! return a file name that matches it, and keep on doing so until
       ! it runs out of names.  On reaching the last name in the series
       ! it will return lastf as true.
                                     
       ! UNIX version, written by Koji in 10 minutes.

       ! Changed by Timn to close the unit it writes the file names to
       ! when it's finished, and not to write out "No match." to the
       ! terminal if no files match.

       ! Purists version for f95.

       ! This statement is needed for the NAG f95 compiler.
       ! use f90_unix_env

       implicit none
                                             
       ! Passed variables.
       character(len=*) :: fildef, lstfil
       character(len=3) :: extdef
       logical :: lastf
                     
       character(len=128), save :: olddef, safe, nextfil, scratch
       character(len=128), external :: addon
       integer :: iflag, iostat
       integer, save :: lun

       character(len=17) :: filnam

       ! This statement is needed when not using the NAG f95 compiler.
       ! integer, external :: getpid
                  
       ! Trap a blank file specification.
       if (fildef.eq.' ' .and. (extdef(1:1).eq.' ' .or. &
       extdef(1:1).eq.char(0))) goto 110

       lastf=.false.
       if( addon(fildef, extdef) .ne. olddef ) then
         ! The input specification is new, save it for later comparison.
         olddef = addon(fildef, extdef)
         ! Now expand the file name.
         call unix_extend( fildef, safe, iflag )
         call find_free( lun )
         filnam='/tmp/wildf.'
         ! filnam='~/wildf.'
         write(filnam(12:17),'(i6.6)') getpid()
         ! You need this in case you go back to writing it in the home
         ! directory.
         call UNIX_EXTEND( filnam, scratch, iflag )
         call arkcom( '/bin/ls ' // &
         addon(safe, extdef) // '> ' // scratch )
         ! * It will crash if you don't have write permission to 
         ! your home directory, or if the disk is full.
         open( lun, file=scratch, status='old' )
         ! Now the file 'wildf.tmp' contains all the file names that
         ! correspond to the input specification, one per line.
         read( lun, '(a)', end=100) nextfil
       end if
       
       lstfil=nextfil
       ! more or less
       read(lun, '(a)', iostat=iostat) nextfil
       if (iostat < 0) then
       !   The current 'lstfil' is the last one.
         lastf = .true.
         close (lun)
         call arkcom( '/bin/rm ' // scratch )
         olddef = ' '
       end if
       goto 900

       ! Case where no files exist that satisfy the input
       ! specification.
100    continue
       close (lun)
       call arkcom( '/bin/rm ' // scratch )
110    lstfil = fildef
       lastf = .true.
       olddef = ' '

900    end
