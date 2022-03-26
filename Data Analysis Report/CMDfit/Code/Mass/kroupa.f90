module mass_functions

  implicit none

  contains

  subroutine calc_kroupa_consts(divide, alpha1, alpha2, low_mass, up_mass, const1, const2, divide_frac)

    ! This routine calculates the constants to ensure the Kroupa mass function integrates to one.

    real, intent(in) :: divide, alpha1, alpha2, low_mass, up_mass
    real, intent(inout) :: const1, const2
    ! The fraction of stars which lie below the power law break.
    real, intent(inout) :: divide_frac

    real, save :: old_divide=0.0, old_alpha1=0.0, old_alpha2=0.0, old_up_mass=0.0, old_low_mass=0.0
    logical, save :: first=.true.
    real :: norm

    if (first .or. &
      abs(old_divide-divide) > 2.0*tiny(divide) .or. &
      abs(old_alpha1-alpha1) > 2.0*tiny(alpha1) .or. &
      abs(old_alpha2-alpha2) > 2.0*tiny(alpha2) .or. &
      abs(old_up_mass-up_mass) > 2.0*tiny(up_mass) .or. &
      abs(old_low_mass-low_mass) > 2.0*tiny(low_mass)) then

      const1=1.0
      const2=const1*(divide**(alpha2-alpha1))

      ! Now calculate the normalisation.
      if (up_mass < divide) then
        norm = (const1/(1.0-alpha1))*(up_mass**(1.0-alpha1) - low_mass**(1.0-alpha1))
        divide_frac = 1.0
      else if (low_mass > divide) then
        norm = (const2/(1.0-alpha2))*(up_mass**(1.0-alpha2) - low_mass**(1.0-alpha2))
        divide_frac=0.0
      else
        divide_frac =        (const1/(1.0-alpha1))*(divide **(1.0-alpha1) - low_mass**(1.0-alpha1))
        norm = divide_frac + (const2/(1.0-alpha2))*(up_mass**(1.0-alpha2) - divide  **(1.0-alpha2))
        divide_frac=divide_frac/norm
      end if

      const1=const1/norm
      const2=const2/norm

      first=.false.
      old_divide=divide
      old_alpha1=alpha1
      old_alpha2=alpha2
      old_up_mass=up_mass
      old_low_mass=low_mass
      
    end if

  end subroutine calc_kroupa_consts

  real function kroupa_star_mass(low_mass, up_mass)

    ! This function returns a star mass drawn from the Kroupa function.
    
    real, intent(in) :: up_mass, low_mass

    real :: frac
    real, parameter :: divide=0.5, alpha1=1.3, alpha2=2.3
    real, save :: const1, const2, divide_frac

    call calc_kroupa_consts(divide, alpha1, alpha2, low_mass, up_mass, const1, const2, divide_frac)
  
    call random_number(frac)

    if (frac < divide_frac) then
      ! This lies below the power-law break.
      kroupa_star_mass = ( low_mass**(1.0-alpha1)  + (frac*(1.0-alpha1)/const1) )**(1.0/(1.0-alpha1))
    else
      frac=frac-divide_frac
      kroupa_star_mass = ( divide**(1.0-alpha2)  + (frac*(1.0-alpha2)/const2) )**(1.0/(1.0-alpha2))
    end if

  end function kroupa_star_mass


end module mass_functions

  real function fmake(mass)

    ! This function (for test purposes) returns the Kroupa mass function.

    use mass_functions
    
    implicit none

    real, intent(in) :: mass

    real, parameter :: divide=0.5, alpha1=1.3, alpha2=2.3
    real, parameter :: up_mass=10.0, low_mass=0.1
    real :: const1, const2, divide_frac

    call calc_kroupa_consts(divide, alpha1, alpha2, low_mass, up_mass, const1, const2, divide_frac)

    if (mass < divide) then
      fmake=const1/(mass**alpha1)
    else
      fmake=const2/(mass**alpha2)
    end if

  end function fmake
  
