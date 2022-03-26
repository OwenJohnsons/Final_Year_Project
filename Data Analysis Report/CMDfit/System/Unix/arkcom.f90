        subroutine arkcom(string)

!       Spawns an operating system command.

!       This statement is needed for the NAG f95 compiler.
!       use f90_unix_proc
                                 
        character(len=*) :: string
                            
        call system(string)
                              
        end subroutine arkcom
