module random

  implicit none

contains

  real function rnd_gauss(mean, sigma)

    ! This replaces Koji's old Gaussian routine.  It produces much
    ! less structure.



    real, intent(in) :: mean, sigma
    real :: x
    real, save :: y
    logical, save :: odd=.true.

    if (odd) then
      call grnf(x,y)
      rnd_gauss=mean+sigma*x
      odd=.false.
    else
      rnd_gauss=mean+sigma*y
      odd=.true.
    end if


  end function rnd_gauss

  SUBROUTINE GRNF (X,Y)

    ! Two Gaussian random numbers generated from two uniform random
    ! numbers. Copyright (c) Tao Pang 1997.
    ! Minor change by Timn to use random_number

    IMPLICIT NONE
    REAL, INTENT (OUT) :: X,Y
    REAL :: PI,R1,R2,harvest

    PI = 4.0*ATAN(1.0)
    call random_number(harvest)
    R1 = -ALOG(1.0-harvest)
    call random_number(harvest)
    R2 = 2.0*PI*harvest
    R1 = SQRT(2.0*R1)
    X  = R1*COS(R2)
    Y  = R1*SIN(R2)

  END SUBROUTINE GRNF

subroutine set_seed()

  ! A portable way of resetting the random number seed for Fortran 90.  
  ! Taken from http://web.ph.surrey.ac.uk/fortweb/index.html.

  integer :: i_seed
  integer, dimension(:), allocatable :: a_seed
  integer, dimension(1:8) :: dt_seed

  call random_seed(size=i_seed)
  allocate(a_seed(1:i_seed))
  call random_seed(get=a_seed)
  call date_and_time(values=dt_seed)
  a_seed(i_seed)=dt_seed(8)
  a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
  call random_seed(put=a_seed)
  deallocate(a_seed)

end subroutine set_seed

end module random
