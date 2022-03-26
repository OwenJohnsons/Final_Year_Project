        subroutine arkgcl(cbuf)

!       Returns the entire command line in cbuf.
!       Based on the XANADU subroutine rdforn.

!       This statement is needed for the NAG f95 compiler.
!       use f90_unix_env

        character(len=*), intent(out) :: cbuf

        integer narg, i, j, lbuf

        narg=iargc()
        cbuf=' '
        lbuf=0
        if (narg .gt. 0) then
          do 190 i=1, narg
            call getarg(i, cbuf(lbuf+1:len(cbuf)))
!           Add a space between arguments, and update lbuf.
            do 170 j=len(cbuf), 1, -1
              if(cbuf(j:j).ne.' ' .and. cbuf(j:j).ne.char(0)) then
                cbuf(j+1:j+1)=' '
                lbuf=j+1
                goto 190
              end if
170         continue
190       continue
        end if

        end subroutine arkgcl
