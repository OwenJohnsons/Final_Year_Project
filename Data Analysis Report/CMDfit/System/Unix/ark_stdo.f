       subroutine ARK_STDO( )
C
C       This is an optional subroutine, to be called near the beginning of
C       a program, to "beautify" the terminal output for those program
C       that assumes Fortran Carriage Control.  Needs a blank subroutine
C       in vms_system.for; no harm if not called, other than the blank first
C       column.
C
C
C       This is the ULTRIX version.
c        open (unit=6, status='unknown', carriagecontrol='fortran')
       end
