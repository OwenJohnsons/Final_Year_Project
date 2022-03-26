real function fmake(mass)

  implicit none

  real, intent(in) :: mass

  ! Giles' version.
  fmake = (1.0/mass)*exp(-1.0*((log10(mass)-log10(0.079))**2.0)/(2.0*0.69))

  ! Dabringhausen, Hilker & Kroupa (2008) version.
  ! fmake = (1.0/mass)*exp(-1.0*((log10(mass)-log10(0.080))**2.0)/0.9522)

end function fmake
