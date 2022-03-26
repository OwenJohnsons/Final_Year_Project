module put_header

  ! Subroutines in this file add one line to an ARK_FITS header.
  ! see the subroutines in GETDROSS.FOR to recover the info.
  ! If header line of the same KEY value exists, you can over-write it
  ! (Override=1), add another line (Override=0) or stop writing it
  ! (Override=-1).

  use get_header

  implicit none

  integer, private, save :: ncolumns_header=0

contains

  Subroutine CLEAR_HEADER( )

    Integer K

    Do K = 1, size(dross)
      call blklin(dross(k))
    end Do
    dross(1)(1:8)='END     '

    ncolumns_header=0

  End subroutine clear_header

  subroutine put_header_column(type, form, comm)

    ! Writes the basic descriptors for a FITS column.

    character(len=*), intent(in) :: type, form, comm

    ncolumns_header=ncolumns_header+1
    call put_header_column_s('TTYPE', ncolumns_header, type) 
    call put_header_column_s('TFORM', ncolumns_header, form) 
    call put_header_column_s('TCOMM', ncolumns_header, comm) 

  end subroutine put_header_column

  subroutine put_header_column_s(key, column, value)

    ! Writes a keyword for the column-th column.

    character(len=*), intent(in) :: key, value
    integer, intent(in) :: column
    
    character(len=8) :: full_key, value8
    integer :: end

    full_key=key
    end=len_trim(full_key)
    if (ncolumns_header<0 .or. ncolumns_header>999) then
      print*, 'put_header_column_s cannot deal with ncolumns_header ', &
      ncolumns_header
    else if (ncolumns_header < 10) then
      write(full_key(end+1:end+1), '(i1)') ncolumns_header
    else if (ncolumns_header < 100) then
      write(full_key(end+1:end+2), '(i2)') ncolumns_header
    else 
      write(full_key(end+1:end+3), '(i3)') ncolumns_header
    end if

    if (len_trim(value) > 8) then
      call put_header_s(full_key, value, ' ', 0)
    else
      value8='        '
      value8(1:len_trim(value))=trim(value)
      call put_header_s(full_key, value8, ' ', 0)
    end if

  end subroutine put_header_column_s

  Subroutine put_header_c( Key, Value, Descriptor, Override )

    Character*( * ) Key, Descriptor
    Character Value
    Integer Override

    Character*8 out_Key
    Integer I, L

    Call MATCH_KEY( Key, out_Key, L, Override )
    If( L .le. 0 ) Return
    call blklin(Dross( L ))
    Dross( L )( 1: 8 ) = out_Key
    Dross( L )( 9: 9 ) = '='
    I = Ichar( Value )
    If( I .ge. 97 .and.  I .le. 122 ) Value = Char( I - 32 )
    Dross( L )( 30: 30 ) = Value
    If( Descriptor .eq. ' ' ) Return
    Dross( L )( 32: 32 ) = '/'
    I = Min( Len( Descriptor ) + 32, 80 )
    Dross( L )( 33: I ) = Descriptor( 1: I - 32 )

900 End subroutine put_header_c



    Subroutine Put_header_r( Key, Value, Descriptor, Override )

      Character*( * ) Key, Descriptor
      Real Value
      Integer Override

      Character*8 out_Key
      Integer I, L

      Call MATCH_KEY( Key, out_Key, L, Override )
      If( L .le. 0 ) Return
      call blklin(Dross( L ))
      Dross( L )( 1: 8 ) = out_Key
      Dross( L )( 9: 9 ) = '='
      If( Value .gt. 1.0E+04 .or. Value .lt. -1.0E+04 ) Then
        Write( Dross( L )( 15: 30 ), '(E16.10)' ) Value
      Else
        Write( Dross( L )( 15: 30 ), '(E16.10)' ) Value
      End If
      If( Descriptor .eq. ' ' ) Return
      Dross( L )( 32: 32 ) = '/'
      I = Min( Len( Descriptor ) + 32, 80 )
      Dross( L )( 33: I ) = Descriptor( 1: I - 32 )

    End subroutine put_header_r



    Subroutine Put_header_d( Key, Value, Descriptor, Override )

      Character*( * ) Key, Descriptor
      double precision Value
      Integer Override

      Character*8 out_Key
      Integer I, L

      Call MATCH_KEY( Key, out_Key, L, Override )
      If( L .le. 0 ) Return
      call blklin(Dross( L ))
      Dross( L )( 1: 8 ) = out_Key
      Dross( L )( 9: 9 ) = '='
      If( Value .gt. 1.0E+04 .or. Value .lt. -1.0E+04 ) Then
        Write( Dross( L )( 15: 30 ), '(E16.10)' ) Value
      Else
        Write( Dross( L )( 15: 30 ), '(E16.10)' ) Value
      End If
      If( Descriptor .eq. ' ' ) Return
      Dross( L )( 32: 32 ) = '/'
      I = Min( Len( Descriptor ) + 32, 80 )
      Dross( L )( 33: I ) = Descriptor( 1: I - 32 )

    End subroutine put_header_d



    Subroutine put_header_i( Key, Value, Descriptor, Override )

      Character*( * ) Key, Descriptor
      Integer Value
      Integer Override

      Character*8 out_Key
      Integer I, L

      Call MATCH_KEY( Key, out_Key, L, Override )
      If( L .le. 0 ) Return
      call blklin(Dross( L ))
      Dross( L )( 1: 8 ) = out_Key
      Dross( L )( 9: 9 ) = '='
      Write( Dross( L )( 20: 30 ), '(I11)' ) Value
      If( Descriptor .eq. ' ' ) Return
      Dross( L )( 32: 32 ) = '/'
      I = Min( Len( Descriptor ) + 32, 80 )
      Dross( L )( 33: I ) = Descriptor( 1: I - 32 )

    End subroutine put_header_i



    Subroutine put_header_s( Key, Value, Descriptor, Override )

      Character*( * ) Key, Descriptor
      Character*( * ) Value
      Integer Override

      Character*8 out_Key
      Integer I, K, L, N, Length, Pos

      Call MATCH_KEY( Key, out_Key, L, Override )
      If( L .le. 0 ) Return
      call blklin(Dross( L ))
      Dross( L )( 1: 8 ) = out_Key
      Dross( L )( 9: 11 ) = '= '''
      I = Len( Value )
      ! Make all blank strings have 1 character.
      Length = 1
      Do K = 1, I
        N = Ichar( Value( K: K ) )
        If( N .ge. 33 ) Then
          Length = K
          !           If( N .ge. 97 .and. N .le. 122) Value( K: K ) = Char( N - 32 )
        End If
      End Do
      Length = Min( Length, 68 )
      Pos = Length + 11
      Dross( L )( 12: Pos ) = Value( 1: Length )
      Pos = Pos + 1
      Dross( L )( Pos: Pos ) = ''''
      If( Pos .ge. 71 .or. Descriptor .eq. ' ' ) Return
      Pos = Max( 32, ( Pos + 9 ) / 8 * 8 )
      Dross( L )( Pos: Pos ) = '/'
      I = Min( Len( Descriptor ) + Pos, 80 )
      Dross( L )( Pos + 1: I ) = Descriptor( 1: I - Pos )

    End subroutine put_header_s



    Subroutine MATCH_KEY( Key, out_Key, L, Override )

      ! Tells the put_header_% routines where to write the line of header,
      ! or tells then -1 if they are not to.  This version emulates Koji's
      ! original using matkey, and also writes a 'END' to the header.

      Character*( * ) Key
      Character*8 out_Key
      Integer L, Override

      integer iend

      l=matkey(Key, out_Key)
      iend=lndrss()
      if (override .eq. 0) then
        if (out_key .eq. 'SIMPLE  ' &
        .or. out_key .eq. 'BITPIX  ' &
        .or. out_key .eq. 'REALTYPE' &
        .or. out_key .eq. 'NAXIS   ' &
        .or. out_key .eq. 'NAXIS1  ' &
        .or. out_key .eq. 'NAXIS2  ' &
        .or. out_key .eq. 'NAXIS3  ') then
          print*, 'Error in ARK S/R MATCH_KEY, attempt to add a'
          print*, 'second line of header with the keyword ', out_key
          print*, 'aborted.  Please call PUT_HEADER routine with'
          print*, 'override value other than 0.'
          goto 900
        end if
      end if

      !       Check for manditory header items.  If the item is manditory,
      !       put it in the right place and delete any other versions of it.
      if (     out_key.eq.'SIMPLE  ' .and. l.ne.1) then
        call sdown(1)
        l=1
        call rem_header('SIMPLE  ')
        iend=max(iend+1,2)
      else if (out_key.eq.'BITPIX  ' .and. l.ne.2) then
        call sdown(2)
        l=2
        call rem_header('BITPIX  ')
        iend=max(iend+1,3)
      else if (out_key.eq.'REALTYPE' .and. l.ne.2) then
        call sdown(2)
        l=2
        call rem_header('REALTYPE')
        iend=max(iend+1,3)
      else if (out_key.eq.'NAXIS   ' .and. l.ne.3) then
        call sdown(3)
        l=3
        call rem_header('NAXIS   ')
        iend=max(iend+1,4)
      else if (out_key.eq.'NAXIS1  ' .and. l.ne.4) then
        call sdown(4)
        l=4
        call rem_header('NAXIS1  ')
        iend=max(iend+1,5)
      else if (out_key.eq.'NAXIS2  ' .and. l.ne.5) then
        call sdown(5)
        l=5
        call rem_header('NAXIS2  ')
        iend=max(iend+1,6)
      else if (out_key.eq.'NAXIS3  ' .and. l.ne.6) then
        call sdown(6)
        l=6
        call rem_header('NAXIS3  ')
        iend=max(iend+1,7)
      else
        !         l is now equal to the line number of the header item,
        !         or negative if the item doesn't exist.              
        if (l.gt.0 .and. override.eq.-1) then
          !           The item exists and is not to be overwritten.
          l=-1
        end if
        if (l.lt.0 .or. override.eq.0) then
          !           If the line of header does not exist, or an extra
          !           line is to be written even if it does, then add
          !           a line of header.
          !           Check for a full header.
          if (iend .ge. size(dross)) then
            print*, 'Warning, header is full, overwriting item ', &
            dross(179)(1:8)
            l=179
          else
            l=iend
            iend=iend+1
          end if
        end if
      end if

      call blklin(dross(iend))
      dross(iend)(1:8)='END     '

900 end Subroutine MATCH_KEY

    subroutine sdown(itop)

      ! Moves the header down a line, to make room for a manditory item.

      integer itop

      integer i

      if (lndrss() .ge. size(dross)) then
        print*, 'Warning, header is full, overwriting item ', &
        dross(179)(1:8)
        call blklin(dross(179))
        dross(179)(1:3)='END'
      end if

      do 190 i=179, itop, -1
        dross(i+1)=dross(i)
190     continue
        call blklin(dross(itop))

      end subroutine sdown

      subroutine rem_header(key)

        ! Removes all the lines of header matching the key given.

        character key*(*), outkey*8
        integer i, iend

100     iend=matkey(key,outkey)
        if (iend .lt. 0) goto 900

        do 780 i=iend+1, size(dross)
          dross(i-1)=dross(i)
780       continue
          call blklin(dross(size(dross)))

          goto 100

900     end subroutine rem_header


        subroutine blklin(string)

          !       Makes sure every character in a string is a blank.

          character*(*) string

          integer i

          do 300 i=1, len(string)
            string(i:i)=' '
300         continue

          end subroutine blklin

        end module put_header
