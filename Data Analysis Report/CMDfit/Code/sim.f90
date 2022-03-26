program sim

  ! Creates a simulated cluster CMD.

  use define_star
  use quad
  use cmdfit_subs
  use random

  implicit none

  type(a_star), dimension(:), allocatable :: star

  ! First set up the uncertainties to be used in each colour.  icorr is 
  ! between the errors.  Allowed values are:
  ! icorr =  1, for a CMD such as V, V-I;
  ! icorr = -1, for a CMD such as B, U-B;
  ! icorr =  0, gives you uncorrelated errors.
  ! The algebra is set up so that the uncertainties you supply are the total 
  ! uncertainty ignoring correlation (i.e. what most people quote).
  ! This means that for any icorr other than zero, the uncertainty in 
  ! magnitude must be less than that in colour.
  ! There are two ways to set up the uncertainties.  Choose which by
  ! setting this logical.
  logical, parameter :: fixed_uncer = .false.
  ! If this is set true, the uncertainties are equal for all data points, 
  ! and set by these two parameters.
  real, parameter :: fix_col_err=0.043, fix_mag_err=0.03
  ! Otherwise the magnitude uncertainty is taken to be of the form
  ! sqrt(mag_ind**2.0 (norm_10*10.0**(mag_norm*(true_mag-mag_shift)))**2.0),
  ! where mag_ind is a magnitude independent uncertainty which dominates at
  ! bright magnitudes, norm_10 and mag_norm are normalisation terms for
  ! what amounts to sky-limit term, and mag_shift indicates at what 
  ! magnitude the sky noise begins to dominate.
  ! real, parameter :: mag_ind=0.005, norm_10=0.03, mag_norm=0.3, mag_shift=11.5
  real, parameter :: mag_ind=0.02, norm_10=0.03, mag_norm=0.4, mag_shift=20
  ! Finally, the uncertainty in colour is taken to be a constant times this.
  real, parameter :: col_ratio=1.3
  integer :: icorr=0
  integer :: iband, iminus
  real :: age, dm, low_mass, high_mass

  integer :: istar
  real :: harvest
  character(len=50) :: ofname
  real :: v_err, i_err
  integer :: multiplicity
  real :: binary_fraction
  integer :: nstars
  real :: p_mass
  character(len=20) :: isoflag

  character(len=20), dimension(2) :: colnam
  character(len=200) :: iso_file, bc_file, ext_file
  real, parameter :: red=0.0

  call choose_models(iso_file, bc_file, ext_file, colnam)

  print*, '> Give binary fraction.'
  read(*,*) binary_fraction

  print*, '> Give output catalogue file name.'
  read(*,*) ofname

  print*, '> Give the required log10(age).'
  read(*,*) age

  print*, '> Give the required distance modulus.'
  read(*,*) dm

  ! The range used in the Naylor and Jeffries paper was 0.36 1.9Mo
  print*, '> Give range of masses (low, high).'
  read(*,*) low_mass, high_mass

  ! The number of stars used in Naylor and Jeffries was 100.
  print*, 'How many stars do you want?'
  read(*,*) nstars

  ! Now lets sort out whether you add or subtract the colour from the
  ! magnitude to get the other magnitude.
  iminus=index(colnam(2), '-')
  if (iminus == 0) then
    print*, 'Cannot find the minus sign in colour ', colnam(2)
    stop
  end if
  iband=index(colnam(2), trim(colnam(1)))
  if (iband == 0) then
    print*, 'Cannot find the magnitude ', trim(colnam(2)), &
    ' in the colour ', trim(colnam(2))
    print*, 'So assuming the colours are uncorrelated.'
    icorr=0
  else if (iband > iminus) then
    print*, 'Will add ', trim(colnam(1)), ' to ', trim(colnam(2)), &
         ' to create ', colnam(2)(1:iminus-1)
    icorr = -1
  else
    print*, 'Will subtract ', trim(colnam(2)), ' from ', &
         trim(colnam(1)), &
         ' to create ', trim(colnam(2)(iminus+1:len(colnam(2))))
    icorr = 1
  end if

  open(unit=21, file='sim.iso')
  write(21,*) 'Columns are color, mag, primary mass, flag and multiplicity.'
  write(21,*)
  write(21,*)


  allocate(star(nstars))
  call zero_star(star)
  call set_seed()

  do istar=1, nstars

10  call make_star(age, iso_file, bc_file, ext_file, colnam, 'kroupa', &
    low_mass, high_mass, 0.0, binary_fraction, red, &
    star(istar)%col(1)%data, star(istar)%col(2)%data, &
    p_mass, multiplicity, isoflag)

    write(21,*) star(istar)%col(2)%data, star(istar)%col(1)%data, p_mass, &
    trim(isoflag), multiplicity
    if (trim(isoflag) /= 'OK') then
      print*, 'Star ', istar, ' of mass ', p_mass, ' has flag ', trim(isoflag)
      print*, 'Trying again.'
      goto 10
    end if

    ! Deal with the distance modulus.
    star(istar)%col(1)%data=star(istar)%col(1)%data+dm

    ! Some old code which dealt with colour-colour diagrams.
    ! delta_b=rnd_gauss(0.0, ubv_uncer)
    ! delta_v=rnd_gauss(0.0, ubv_uncer)
    ! star(istar)%col(1)%data= v_for_ubv()+delta_v
    ! star(istar)%col(1)%err = ubv_uncer
    ! star(istar)%col(1)%flg = 'OO'
    ! star(istar)%col(2)%data= true_col+delta_b-delta_v
    ! star(istar)%col(2)%err = ubv_uncer
    ! star(istar)%col(2)%flg = 'OO'
    ! star(istar)%col(3)%data= true_mag+rnd_gauss(0.0, ubv_uncer)-delta_b
    ! star(istar)%col(3)%err = ubv_uncer
    ! star(istar)%col(3)%flg = 'OO'

    ! Now, assume the CMD is of the form V vs V-I.  Then we must 
    ! correlate the uncertainties in the following way.
    ! First get the error in each magnitude.
    if (fixed_uncer) then
      star(istar)%col(1)%err=fix_mag_err
      star(istar)%col(2)%err=fix_col_err
    else
      star(istar)%col(1)%err=sqrt(mag_ind**2.0 &
      + (norm_10*10.0**(mag_norm*(star(istar)%col(1)%data-mag_shift)))**2.0)
      star(istar)%col(2)%err=col_ratio*star(istar)%col(1)%err
    end if
    v_err=rnd_gauss(0.0, star(istar)%col(1)%err)
    i_err=rnd_gauss(0.0, sqrt(star(istar)%col(2)%err**2.0 - &
    (star(istar)%col(1)%err*real(icorr))**2.0))
    ! Now for the magnitude this is easy.
    star(istar)%col(1)%data=star(istar)%col(1)%data+v_err
    ! And for the colour.
    if (icorr == 0) then
      star(istar)%col(2)%data=star(istar)%col(2)%data+i_err
    else
      star(istar)%col(2)%data=star(istar)%col(2)%data+real(icorr)*(v_err-i_err)
    end if

    star(istar)%col(1)%flg ='OO'
    star(istar)%col(2)%flg='OO'

  end do
  call nxtcls_out(ofname)
  call write_cluster_file(2, colnam, star, nstars)
  close(21)

end program sim
     
