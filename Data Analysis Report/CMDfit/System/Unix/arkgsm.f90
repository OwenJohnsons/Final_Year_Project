        logical function arkgsm(symbol, value)
                                              
!       Gets a DCL global symbol or a UNIX enviroment variable.

        character(len=*), intent(in) :: symbol
        character(len=*), intent(out) :: value

        character(len=30) :: work

        work=symbol
        call upcase(work)
                                       
        call getenv(trim(work), value)
        if (value(1:1).ne.' ' .and. value(1:1).ne.char(0)) then
          arkgsm=.true.
        else
          arkgsm=.false.
        end if
              
        end function arkgsm
