      subroutine writen(prompt)

      character(len=*) :: prompt

      write(*,'(a)', advance='no') prompt(1:len_trim(prompt))//' '

      end subroutine writen
