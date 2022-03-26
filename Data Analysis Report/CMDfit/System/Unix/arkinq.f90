        subroutine arkinq(file, filex, exist, arktyp)

        implicit none
                                                     
!       Test to see if a file is there.  This is a SUN version.
                                                                    
!       Input values.
!       The file name, and the default extension to be used.
        character(len=*) :: file, filex
                                 
!       Returned values.
!       Logical for existance of file.
        logical :: exist
!       If possible, find out if the file is ARK, ASCII or ARK_FITS.
        character(len=*) :: arktyp
                            
        integer :: itest, iunit, iflag
        character(len=200) :: addon, work
!       The debugger.
        logical debug
        common /bug/ debug
                          
        arktyp=' '
        call unix_extend(file, work, iflag)
        if (debug) print*, '@ S/R ARKINQ thinks ', addon(work, filex)
        inquire(file=addon(work, filex), exist=exist)
!       Trap out what seems to be a bug in DEC FORTRAN inquires.
        if (addon(work, filex) .eq. ' ') exist=.false.
        if (exist) then
          call find_free(iunit)
          open(iunit, file=addon(work, filex), action='read', &
          form='unformatted', access='direct', recl=4, status='old')
!         Not using ARKOPN 'coz this is in system.a
          Read(iunit, rec=1, err=120) itest
!         An error means the file is formatted as opposed to unformatted.
!         ARK header is 2884 byte unformatted, sequential
!         UNIX appends a 4-byte number containing the record length.
          if (itest .eq. 2884) arktyp='ARK'
!         For DEC OSF/1
          if (itest .eq. 1347242323) arktyp='ARK_FITS'
!         For Solaris.
          if (itest .eq. 1397312848) arktyp='ARK_FITS'
120       if (debug) print*, '@ is of type ', arktyp
          if (debug) print*, '@ itest in s/r ARKOPN is ', itest
          close (iunit)
        else
          if (debug) print*, '@ Doesn''t exist.'
        end if
              
        end subroutine arkinq
