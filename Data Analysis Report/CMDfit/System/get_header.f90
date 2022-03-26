      module get_header

!       The standard set of routines for getting headers from ARK files.
                                                                        
!       The routines are functions of the form;
!       integer get_header_%(key, value), where key is a character*(*) string
!       and value and % take the following meanings.
!       %     Value
!       -     -----
!       S     character*(*)
!       C     character*1
!       R     real
!       D     double precision
!       I     integer
!
!       The value of the function is the line number where the key was
!       found unless;
!       (1) If the key was not found -1 is returned.
!       (2) The key was found, but there was an error reading the value,
!           then -2 is returned.
!       If the input string is less than 8 characters, blanks are matched
!       for the remaining characters.  The match is performed in upper case,
!       translating the input string if need be.
                                                
!       History.
!       Original by Koji.
!       Get_header_s added by Tim.
!       Changes made to make it general for ARK and ARKFITS files:
!         1) Dross increased to 180 lines.
!         2) Keys made to be any length.
!         3) get_header_d added.
!       Made into Fortran 90 by Timn.

        implicit none

        integer, parameter, public :: mdross=360
        character(len=80), dimension(mdross), public :: dross
        data dross(1) /'END'/
        common /pasdrs/ dross
        
      contains
                                
       Integer Function GET_HEADER_S( Key, Value )

       Character*(*) Value, Key

       Integer I
        character*8 expkey

        i=matkey(key, expkey)
        if (i .gt. 0) then
         read(Dross( I )(11:80), *, err=100) value
          goto 200
!           If someone has missed the quotes of the string, read it this way.
100         read(Dross( I )(11:80), '(a)', err=800) value
200       continue
         If ( Value .gt. ' ' ) Then
           GET_HEADER_S = i
         Else
           GET_HEADER_S = -2
         End If
       else
         GET_HEADER_S = -1
        end if
        goto 900
                
800     get_header_s = -2
                         
900       End Function GET_HEADER_S



       Integer Function GET_HEADER_C(key, value)

       Character Value
        Character*(*) key

       Integer I, L
        character expkey*8

        i=matkey(key, expkey)
        if (i .gt. 0) then
         Value = Dross( I )( 30: 30 )
         If ( Value .gt. ' ' ) Then
           GET_HEADER_C = i
         Else
           GET_HEADER_C = -2
         End If
       else
         GET_HEADER_C = -1
        end if
              
900       End Function GET_HEADER_C


       Integer Function GET_HEADER_R(key, value)

       Character*(*) Key
       Real Value

       Integer I
       character expkey*8

       i = matkey(key, expkey)
       if (i .gt. 0) then
         Read( Dross( I )( 11: 30 ), *, Err = 100 ) Value
         GET_HEADER_R = i
       else
          GET_HEADER_R = -1
       end if
       goto 900
                
!      Condition if fail to read the header.
100    get_header_r = -2

900    End Function GET_HEADER_R



       Integer Function GET_HEADER_D( Key, Value )

       Character*(*) Key
       double precision Value

       Integer I, L
        character expkey*8

        i = matkey(key,expkey)
        if (i .gt. 0) then
         Read( Dross( I )( 11: 30 ), *, Err = 100 ) Value
         GET_HEADER_D = i
        else
          GET_HEADER_D = -1
        end if
        goto 900
                
!       Condition if fail to read the header.
100     get_header_D = -2

900       End Function GET_HEADER_D



       Integer Function GET_HEADER_I( Key, Value )

       Character*(*) Key
       Integer Value

       Integer I, L
        character expkey*8

        i = matkey(key,expkey)
        if (i .gt. 0) then
         Read( Dross( I )( 11: 30 ), *, Err = 100 ) Value
         GET_HEADER_i = i
        else
          GET_HEADER_i = -1
        end if
        goto 900
                
!       Condition if fail to read the header.
100     get_header_i = -i

900       End function GET_HEADER_I


        integer function matkey(key,expkey)
                                           
!       Returns the number of the first line whose key matches the
!       input key.  The search is not case sensitive for the input
!       key (but the key in the header must be in upper case).  If the
!       input key is shorter than 8 characters, blanks are matched as
!       rest of the key.  If the key is not found -1 is returned.
!       The upper case, 8 character version of the key is also returned.
                                                                        
        character*(*) key
        character expkey*8
        integer i, n, j
                                              
!       Convert to upper case.
        do i= 1, min(len(key),8)
          n = ichar( Key(i:i) )
          if (n.ge.97 .and. n.le.122 ) n = n - 32
!         And catch nulls.
          if (n .eq. 0) n=32
          expkey(i:i) = char(n)    
        end do
!       Pad with blanks.
        if (len(key) .lt. 8) then
          do 120 i=len(key)+1, 8
            expkey(i:i) = ' '
120       continue
        end if
              
!       Now if we are trying to match 'END' only examine the first 3
!       letters.
        j=8
        if (expkey(1:3) .eq. 'END') j=3
                                       
!       Look for a matching keyword.
        do 310 i=1, lndrss()
          if (dross(i)(1:j) .eq. expkey(1:j)) then
            matkey=i
            goto 900
          end if
310     continue
                
!       If we get to here we've failed to match the key.
        matkey=-1
                 
900     end function matkey
           
           
        Integer Function LNDRSS( )

!       Finds the next line of the header that should be written to.  Thus
!       if there is an 'END' statement, that line is returned, otherwise
!       it returns size of dross.

        integer :: k
 
        lndrss = size(dross)
        end_search: do k=1, size(dross)
          if (dross(k)(1:3) == 'END') then
            lndrss=k
            exit end_search
          end if 
        end do end_search
 
        End Function LNDRSS
      
        subroutine deblank
      
          integer :: k
                  
          do
            k=lndrss()
            if (dross(k-1)(1:8)  == '        ') then
              dross(k-1)='END     '
              dross(k)='        '
            else
              exit
            end if
          end do
          
        end subroutine deblank

      end module get_header
