        subroutine arkopn(iunit, file, filex, stat, pos, form, acces, &
        irecl)

        implicit none
              
!       Opens a file with a machine independent call.            

!       For the access='direct' part of this the DEC OSF/1 version should
!       have recl=irecl/4, but for the Solaris and Linux version recl=irecl.

!       Returned values.
!       The number of the unit opened.
        integer :: iunit
                     
!       Input values.
!       The file name, and the default extension to be used.
        character(len=*) :: file, filex
!       The status to be used (OLD, NEW or UNKNOWN)
        character(len=*) :: stat
!       The position to start writing into the file (APPEND, OVERWRITE,
!       READONLY).
        character(len=*) :: pos
!       The form of the file (FORMATTED, UNFORMATTED or TAPE)
        character(len=*) :: form
!       The file access (DIRECT or SEQUENTIAL), need only be supplied if
!       form is UNFORMATTED.
        character(len=*) :: acces
!       The record length (the BLOCKSIZE is
!       set to the input record length).  Only used if access is 'DIRECT'.
        integer :: irecl
                     
!       Locals.
!       All the above strings translated to upper case.
        character(len=80) :: ufile, ustat, upos, uform, uacces
!       And some similar strings.
        character(len=80) :: uaction
!       The error flag from UNIX_EXTEND
        integer :: iflag
!       The VAX default file extension simulator.
        character(len=80) :: addon
!       An error flag.
        integer :: iostat
!       The debugger.
        logical :: debug
        common /bug/ debug
                          
!       Copy the strings into local variables.
        call unix_extend(file, ufile, iflag)
        ufile = addon(ufile, filex)

        ustat = stat    
        call upcase(ustat)

        uform = form
        call upcase(uform)

        uacces = 'SEQUENTIAL'
        if (uform .eq. 'UNFORMATTED') uacces = acces
        call upcase(uacces)

        uaction='READWRITE'

        upos = pos
        call upcase(upos)
!       We had better check that pos has a legal value, as the operating
!       system won't tell us!
        if (upos .eq. 'READONLY') then
          uaction='READ'
          upos='REWIND'
        else if (upos .eq. 'OVERWRITE') then
          upos='REWIND'     
!         This is a UNIX system which will crash if it attempts to open
!         an extant file with status='new', so cure the problem if the
!         user has specified 'OVERWRITE'.
          if (ustat .eq. 'NEW') ustat='UNKNOWN'
        else if (upos .ne. 'APPEND') then
          print*, 'ARK S/R ARKOPN called in illegal value for POS.'
          print*, 'POS was "'//upos(1:len(pos))//'"'
          goto 900
        end if
              
!       Get a free unit.
        call find_free(iunit)
                             
!       Write out an informational message.
        if (debug) then
          print*, '@ S/R ARKOPN about to open file ', ufile
          print*, '@ with argument list ', stat, ' ',  pos, ' ', &
          form, ' ', acces, ' ', irecl
          print*, '@ and unit number ', iunit
        end if
              
                                                                     
!       The hard work is done here.
        if (uacces .ne. 'DIRECT') then
          open (unit=iunit, status=ustat, form=uform, access=uacces, &
          file=ufile, action=uaction, position=upos, iostat=iostat)
        else
          open (unit=iunit, status=ustat, form=uform, access=uacces, &
          file=ufile, action=uaction, recl=irecl, iostat=iostat)
        end if

        if (iostat /= 0) then
          print*, 'Fortran open statement reports error ', iostat
          print*, 'for file name ', ufile
          stop
        end if
              
900     end subroutine arkopn
