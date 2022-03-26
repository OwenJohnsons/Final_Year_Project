c       The tapeio routines for reading tapes.  Based loosely on what 
c       UNIX should (but often doesn't) provide.  The idea is that you
c       call topen, and then tread until you reach the end of file.
c       You should then call topen before reading the next file (but
c       needn't call tclose).

        integer function topen(itlu, devnam, islabeled)

        integer itlu
        character*(*) devnam
        logical islabeled

c       This array holds the fortran unit numbers corresponding to
c       tape unit numbers 0 to 3.
        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk
        common /chatap/ device
        character*30 device
        integer iflag, tclose
        external tclose

        topen=0
        irec=0
        irw=0

        if (itlu.ge.0 .and. itlu.le.3) then
          if (iunit(itlu) .gt. 0) then
c           The tape unit is already open, let's close it.
            iflag=tclose(itlu)
            if (iflag .ne. 0) then
              print*, 'Error in ARK TAPEIO S/R TOPEN, which called'
              print*, 'S/R TCLOSE, which reported error ', iflag
              stop
            end if
          end if
          call find_free(iunit(itlu))
          dunit(itlu)=devnam
        else
          print*, 'Error in ARK routine topen, itlu = ', itlu
          print*, 'itlu must be between 0 and 3.'
          topen=-1
          goto 900
        end if

        device=devnam

900     end


        subroutine tblock(iblk)

        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk

        integer iblk

        nblk=iblk
                          
        end


        integer function tread(itlu, buffer)

        ! This is needed for the NAG compiler.
        !use f90_unix_proc

        integer itlu
        character*(*) buffer

        integer i

c       This array holds the fortran unit numbers corresponding to
c       tape unit numbers 0 to 3.
        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk
        common /chatap/ device
        character device*30, work*39
        logical there

c       Set tread to the end-of-file condition.
        tread=0

        if (itlu.lt.0 .or. itlu.gt.3) then
          print*, 'Error in ARK routine tread, itlu = ', itlu
          print*, 'itlu must be between 0 and 3.'
          tread = -1
          goto 900
        end if

        if (irw .eq. 2) then
          print*, 'Error in ARK routine tread.  These routines'
          print*, 'cannot support reading a writing to the same'
          print*, 'file.'
        else
          irw = 1
        end if

        if (iunit(itlu) .eq. -1) then
          print*, 'Error in ARK routine tread.  Tape unit ', itlu
          print*, 'has not been opened yet.'
          tread = -1
          goto 900
        end if

        if (irec .eq. 0) then
c         The tape unit has just been opened, and seeing as how a
c         read has been asked for, we'd better read the file off
c         tape.
c         So first see if the file name we are to dump to is being
c         used.
          inquire(file='arktio.tmp', exist=there)
          if (there) then
c           Get rid of it.
            call system('rm -f arktio.tmp')
          end if
c         For some reason you need a <cr> between fortran writes and
c         the output from dd.
          print*, ' '
          work='dd bs=????? of=arktio.tmp conv=sync if='
          write(work(7:11), '(i5.5)') nblk     
          call system(work//device)
          open(unit=iunit(itlu), file='arktio.tmp',
     &    form='unformatted', recl=nblk/4, iostat=i,
     &    access='direct', status='old')
!          open(unit=iunit(itlu), file='arktio.tmp',
!     &    form='unformatted', recl=nblk/4, blocksize=nblk, iostat=i,
!     &    access='direct', status='old')
          if (i .ne. 0) then
            print*, 'Error in ARK subroutine tread.'
            print*, 'Fortran open statement reports error number ', i
            tread=-1
            goto 900
          end if
        end if

        irec=irec+1
        read(iunit(itlu), iostat=i, rec=irec) buffer
c       Check for end-of-file.
        if (i .eq. 36) goto 900
        tread=len(buffer)

        if (i .ne. 0) then
          print*, 'Error in ARK subroutine tread.'
          print*, 'Fortran read statement reports error number ', i
          tread=-1
        end if

900     end


        integer function tclose(itlu)

        ! This is needed for the NAG compiler.
        !use f90_unix_proc

        integer itlu

c       This array holds the fortran unit numbers corresponding to
c       tape unit numbers 0 to 3.
        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk
        common /chatap/ device
        character device*30, work*39

        logical debug
        common /bug/ debug

        tclose=0

        if (itlu.lt.0 .or. itlu.gt.3) then
          print*, 'Error in ARK routine tclose, itlu = ', itlu
          print*, 'itlu must be between 0 and 3.'
          tclose = -1
          goto 900
        end if

        if (iunit(itlu) .ne. -1) close(iunit(itlu))

        if (irw .eq. 2) then
c         For some reason you need a <cr> between fortran writes and
c         the output from dd.
          print*, ' '
          work='dd bs=????? if=arktio.tmp conv=sync of='
          write(work(7:11), '(i5.5)') nblk     
          if (debug) print*, '@ ', work//device
          call system(work//device)
        end if

c       call system ('rm arktio.tmp')

        iunit(itlu)=-1

900     end

        integer function tskipf(itlu, nfiles, nrecs)

        integer itlu, nfiles, nrecs

c       This array holds the fortran unit numbers corresponding to
c       tape unit numbers 0 to 3.
        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk

        integer i
        character*70 spawn_line

        tskipf=0

        if (nfiles .gt. 0) then
          write(spawn_line, 120) dunit(itlu)(1:15), nfiles
120         format( 'mt -f ', A15, ' fsf ', I3.3 )
c         And now set the number of records read to zero.
          irec=0
        else if (nfiles .lt. 0) then
          write(spawn_line, 130) dunit(itlu)(1:15), abs(nfiles)
130         format( 'mt -f ', A15, ' bsf ', I3.3 )
          irec=0
        end if
       call arkcom(spawn_line)

        if (nrecs .gt. 0) then
          do 240 i=1, nrecs
            read(iunit(itlu))
240       continue
        else if (nrecs .lt. 0) then
          print*, 'Error in ARK routine tskipf.'
          print*, 'Attempt to skip a negative number of records.'
          tskipf=-1
        end if

        end

        integer function twrite(itlu, buffer)

        integer itlu
        character*(*) buffer

        integer i

c       This array holds the fortran unit numbers corresponding to
c       tape unit numbers 0 to 3.
        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk

        twrite=0

        if (itlu.lt.0 .or. itlu.gt.3) then
          print*, 'Error in ARK routine twrite, itlu = ', itlu
          print*, 'itlu must be between 0 and 3.'
          twrite = -1
          goto 900
        end if

        if (irw .eq. 1) then
          print*, 'Error in ARK routine twrite.  These routines'
          print*, 'cannot support reading a writing to the same'
          print*, 'file.'
        else
          irw = 2
        end if

        if (iunit(itlu) .eq. -1) then
          print*, 'Error in ARK routine twrite.  Tape unit ', itlu
          print*, 'has not been opened yet.'
          twrite = -1
          goto 900
        end if

        if (irec .eq. 0) then
          open(unit=iunit(itlu), file='arktio.tmp',
     &    form='unformatted', recl=nblk/4, iostat=i,
     &    access='direct', status='unknown')
!          open(unit=iunit(itlu), file='arktio.tmp',
!     &    form='unformatted', recl=nblk/4, blocksize=nblk, iostat=i,
!     &    access='direct', status='unknown')
          if (i .ne. 0) then
            print*, 'Error in ARK subroutine twrite.'
            print*, 'Fortran open statement reports error number ', i
            twrite=-1
            goto 900
          end if
        end if

        irec=irec+1
        write(iunit(itlu), rec=irec, iostat=i) buffer

        if (i .ne. 0) then
          print*, 'Error in ARK subroutine twrite.'
          print*, 'Fortran read statement reports error number ', i
          twrite=-1
        else
          twrite=len(buffer)
        end if

900     end

        block data
        integer iunit(0:3), irec, irw, nblk
        character*30 dunit(0:3)
        common /arktap/ iunit, dunit, irec, irw, nblk
        data iunit /4*-1/
        data nblk /2880/
        end
