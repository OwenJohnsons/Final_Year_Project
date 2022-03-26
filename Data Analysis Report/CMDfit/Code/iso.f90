program iso

  use cmdfit_subs
  use colteff_subs
  use define_star
  use red_ext_interp

  implicit none

  integer :: i, istar
  character(len=256) :: iso_file, bc_file, ext_file
  integer :: i_age, n_age, iostat
  real :: age, start_age, step_age, red, nom_ebmv
  character(len=50) :: ofname
  character(len=8) :: age_name
  integer :: i1, i2, ncomp, ipos
  real :: mass, low_mass, high_mass
  character(len=20) :: mass_function
  logical :: extrap, first, isgood
  real :: binfrac 

  real :: p_mass, p_teff, p_logg, p_lbol
  character(len=50) :: isoflag
  ! The number of stars in the simulation.
  integer, parameter :: nstar=10000
  character(len=20) , dimension(2) :: colnam

  ! The colour the redddening is in.
  character(len=10) :: red_col

  type (a_star) :: star

  character(len=50) :: instring

  print*, 'What range of masses do you want (low, high)?'
  read(*,*) low_mass, high_mass

  get_mass_function: do
    print*, 'What shall I do about stellar masses?'
    print*, '  iso - gives an isochrone equally spaced in mass.'
    print*, '  salpeter - a Salpeter mass function.'
    print*, '  kroupa - a Kroupa mass function.'
    read(*,*) mass_function
    if (trim(mass_function)=='kroupa' .or. trim(mass_function)=='salpeter') then
      get_binfrac: do
        print*, 'What binary fraction do you want?'
        read(*,*,iostat=iostat) binfrac
        if (iostat == 0) exit get_binfrac
        print*, 'Error in input.'
      end do get_binfrac
      exit get_mass_function
    else if (trim(mass_function) == 'iso') then
      binfrac=0.0
      exit get_mass_function
    end if
    print*, 'Error in input.'
  end do get_mass_function

  call choose_models(iso_file, bc_file, ext_file, colnam)

  do
    print*, 'What log10(age) do you want (or start log10(age), logarithmic step and number of isochrones)?'
    read(*,'(a50)',iostat=iostat) instring
    if (iostat /= 0) cycle
    read(instring,*,iostat=iostat) start_age, step_age, n_age
    if (iostat == 0) exit
    step_age=0
    n_age=1
    read(instring,*,iostat=iostat) start_age
    if (iostat == 0) exit
    print*, 'Error in input ', iostat
  end do

  step_age=real(nint(1000.0*step_age))/1000.0

  if (bell_there(ext_file)) then
    print*, 'What E(B-V) you want? (Measured in an energy intergrating system'
    print*, 'with the Bessell filters and the new ODF Kurucz atmospheres.)'
    do
      read(*,*,iostat=iostat) red
      if (iostat == 0) exit
      print*, 'Error in input.'
    end do
    call get_nom_ebmv(red, nom_ebmv, .false., extrap, ext_file)
    if (extrap) then
      print*, 'Fell outside range of E(B-V) tables.'
      stop
    end if
    print*, 'Using E(B-V) is ', nom_ebmv
    red_col='B-V'
  else
    print*, 'Give the required reddening in ', trim(colnam(2))//'.'
    do
      read(*,*,iostat=iostat) red
      if (iostat == 0) exit
      print*, 'Error in input.'
    end do
    red_col=trim(colnam(2))
  end if


  each_age: do i_age=1, n_age

    age=start_age+real(i_age-1)*step_age

    call zero_star(star)
    star%col(1:2)%err=0.01
    star%col(1:2)%flg='OO'
    first=.true.

    each_single_star: do istar=1, nstar

      if (trim(mass_function) == 'iso') then

        mass=low_mass + (high_mass-low_mass)*real(istar)/real(nstar)

        call make_star(age, iso_file, bc_file, ext_file, colnam, 'power_law', mass, mass, -0.1, binfrac, red, &
        star%col(1)%data, star%col(2)%data, p_mass, &
        ncomp, isoflag, p_teff, p_logg, p_lbol)
        mass=low_mass + (high_mass-low_mass)*real(istar)/real(nstar)

      else

        call make_star(age, iso_file, bc_file, ext_file, colnam, mass_function, low_mass, high_mass, -0.1, binfrac, red, &
        star%col(1)%data, star%col(2)%data, p_mass, &
        ncomp, isoflag, p_teff, p_logg, p_lbol)

      end if

      if (first) then
         ! Write out a simple isochrone file.age
        call iso_image_name(colnam, age, ofname)
        ipos=index(ofname, '.fit')
        print*, 'Writing output to ', ofname(1:ipos)//'iso'
        open (unit=12, file=ofname(1:ipos)//'iso')
        write(12,*) '# log10(age) = ', age
        write(12,*) '# Interiors from '//trim(iso_file)
        write(12,*) '# Bolometric corrections from '//trim(bc_file)
        write(12,*) '# Reddening is E('//trim(red_col)//')= ', red
        write(12,*) '# ', trim(colnam(2)), ' ', trim(colnam(1)),  &
        ' mass Teff logG lbol log(age/yr) flag isgood'
        ! And a cluster format catalogue.
        open (unit=21, file=ofname(1:ipos)//'cat')
        write(21,*) ' 2 colours'
        write(21,*) colnam(1), ' ', colnam(2)
        write(21,*)
        first=.false.
      end if

      isgood=.false.
      if (trim(isoflag) == "OK") isgood=.true.
      write(12,*) star%col(2)%data, star%col(1)%data, p_mass, p_teff, &
      p_logg, p_lbol, age, '"'//trim(isoflag)//'"', isgood
      if (trim(isoflag) == 'OK') then
        star%col(1:2)%flg='OO'
      else
        star%col(1:2)%flg='AA'
      end if
      call write_star(21, star, 2)
      
    end do each_single_star
    close(12)
    close(21)

  end do each_age

end program iso
