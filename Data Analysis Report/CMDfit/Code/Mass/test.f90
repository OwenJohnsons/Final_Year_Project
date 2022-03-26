program test

  use mass_functions

  implicit none

  integer, parameter :: nstars=10000
  integer :: istars
  real :: mass

  open(unit=33, file='stars.dat')
  write(33,'(/,/)')

  do istars=1, nstars

    mass=power_law_mass(0.1, 10.0, -2.35)

    write(33,*) real(istars), mass

  end do

end program test

  
