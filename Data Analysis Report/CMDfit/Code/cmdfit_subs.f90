module cmdfit_subs

  use colteff_subs
  use red_ext_interp
  use cmdfit_system
  use likelihood_mod
  use mass_functions

  implicit none

  real, private :: extra_mag, high_tau=-1.0

  contains

  real function high_calc_tau()

    high_calc_tau=high_tau/2.0

  end function high_calc_tau

  subroutine make_star(age, iso_file, bc_file, ext_file, colnam, &
  mass_function, low_mass, up_mass, &
  pwr_mass, bin_frac, red, mag, col, p_mass, &
  nstars, flag_out, p_teff, p_logg, p_lbol)

    ! The age required.
    real, intent(in) :: age
    ! The files to use for the interiors and bolometric corrections.
    character(len=*), intent(in) :: iso_file, bc_file, ext_file
    ! The names of the magnitude and/or colours to be created.
    character(len=20), dimension(2), intent(inout) :: colnam
    ! The form of the mass function to be used.
    character(len=*) :: mass_function
    ! If the form is "power_law" then you need to suppy maximum and
    ! minumum masses.  The function power_law_mass describes pwr_mass.
    ! If the form is "kroupa"  then pwr_mass is ignored.
    ! If the form is "salpeter" then -2.35 is used.
    real, intent(in) :: up_mass, low_mass, pwr_mass 
    ! The fraction of systems which are binaries.
    real, intent(in) :: bin_frac
    ! If the Bell extinction tables are there, the nominal E(B-V), if not the
    ! reddening. 
    real :: red
    ! The colour and magnitude of the stars created. 
    real, intent(out) :: mag, col
    ! The mass of the primary star.
    real, intent(out) :: p_mass
    ! The number of stars in the system.
    integer, intent(out) :: nstars
    ! A character flag, 'OK' if things are.
    character(len=*), intent(out) :: flag_out
    ! The effective temperature and gravity of the primary star.
    ! Calculating these costs time, so if the caller does not need them
    ! it faster to omit them.
    real, optional, intent(out) :: p_teff, p_logg, p_lbol

    ! Locals.
    logical, save :: first=.true.
    real, save :: col1, col2, mag1, mag2, ebmv, av
    real :: harvest, grad
    integer, parameter :: max_multiplicity=2
    real, dimension(max_multiplicity) :: mass, uflux, bflux, vflux, iflux
    real :: umb, bmv, vmi, v, logg, teff, test
    integer :: iminus
    logical :: v_vmi_like, v_bmv_like, umb_bmv_like
    character(len=20) :: flag, work_flag
    logical :: extrap

    type(chrone), save :: iso
    integer :: icol1_1, icol2_1, icol1_2, icol2_2
    logical :: ok
    character(len=50), dimension(:), allocatable :: workflag

    real :: work, red_in1, red_in2

    real, save :: last_age=0.0

    if (abs(age-last_age) > 2.0*tiny(age)) then
      ! Create an isochrone.       
      iso%iso_file=iso_file
      iso%bc_file=bc_file
      iso%age=age
      call iso_calc(iso)
      last_age=age
    end if

    ! First find out if the colours are of the form V vs V-I or V vs B-V.
    umb_bmv_like=.false.
    v_vmi_like=.false.
    v_bmv_like=.false.
    ! First is this colour/colour?
    iminus=index(colnam(1), '-') 
    if (iminus > 0) then
      ! O.K., is it in the correct order?
      if (trim(colnam(1)(iminus+1:len(colnam(1)))) /= &
      colnam(2)(1:index(colnam(2), '-')-1)) then
        print*, 'Colour ', colnam(1), ' not of form U-B vs B-V.'
        stop
      else
        umb_bmv_like=.true.
      end if
    else
      iminus=index(colnam(2), '-') 
      if (iminus == 0) then
        print*, 'Cannot find the minus in ', colnam(2)
        stop
      else
        if (trim(colnam(1)) == colnam(2)(1:iminus-1)) then
          v_vmi_like=.true.
        else if (trim(colnam(1)) == &
        trim(colnam(2)(iminus+1:))) then
          v_bmv_like=.true.
        else
          print*, 'Cannot find ', trim(colnam(1)), ' in ', &
          trim(colnam(2))
          stop
        end if
      end if
    end if

    ok=find_bc_combination(colnam(1), iso%colnam, iso%ncol, icol1_1, icol2_1)
    !if (first) then
    !  print*, ok, icol1_1, icol2_1, iso%colnam(icol1_1)
    !  if (icol2_1 > 0) print*, iso%colnam(icol2_1)
    !end if
    ok=find_bc_combination(colnam(2), iso%colnam, iso%ncol, icol1_2, icol2_2)
    !if (first) then
    !  print*, ok, icol1_2, icol2_2, iso%colnam(icol1_2), iso%colnam(icol2_2)
    !  print*, umb_bmv_like, v_vmi_like, v_bmv_like
    !  first=.false.
    !end if

    nstars=0
    each_component: do

      if (nstars == max_multiplicity) exit each_component

      if (nstars == 0) then 
        nstars=nstars+1
        if (trim(mass_function) == 'power_law') then
          mass(1)=power_law_mass(low_mass, up_mass, pwr_mass)
        else if (trim(mass_function) == 'salpeter') then
          mass(1)=power_law_mass(low_mass, up_mass, -2.35)
        else if (trim(mass_function) == 'kroupa') then
          mass(1)=kroupa_star_mass(low_mass, up_mass)
        else
          print*, 'I do not have mass function form for ', trim(mass_function)
          stop
        end if
        ! Work out the primary's temperature.
        if (present(p_teff)) p_teff = &
        10.0**linintflg(iso%mass, iso%teff, mass(nstars), flag, iso%massflg, &
        iso%teffflg)
        ! And gravity.
        if (present(p_logg)) p_logg = &
        linintflg(iso%mass, iso%logg, mass(nstars), flag, iso%massflg, &
        iso%loggflg)
        if (present(p_lbol)) p_lbol = &
        linintflg(iso%mass, iso%lbol, mass(nstars), flag, iso%massflg, &
        iso%lbolflg)
      else
        call random_number(harvest)
        if (mass(1) < 14.4) then 
          ! First does this star have a companion?  The orginal 
          ! condition (used in the work for the IAU general assembly
          ! paper) was teff < 30,000K, but this condition is better.
          if (harvest > bin_frac) exit each_component
          call random_number(harvest)
        else
          ! For O-stars the binary fraction is probably 67 percent, but
          ! let's be really pessimistic and say 75 percent, assuming
          ! a basic binary fraction of 0.5.
          if (harvest > 1.5*bin_frac) exit each_component
          ! Now the mass ratio.  There is an argument for a mass spike 
          ! that contains up to 20 percent of binaries at q=0.95 to 1.  
          ! Create this spike by losing all the low-q binaries.
          call random_number(harvest)
          if (harvest < 0.2) harvest = 1.0 - harvest/4.0
        end if
        nstars=nstars+1
        ! Now work out the mass of the companion.
        mass(nstars) = mass(nstars-1)*harvest 
        ! Hack to cut off at 0.25Mo for simulation of effect of
        ! binary wedge for original tau^2 paper.
        !if (mass(nstars) < 0.25) then
        !  flux1(nstars)=0.0
        !  flux2(nstars)=0.0
        !  cycle each_component
        !end if
      end if

      ! print*, 'Finding magnitude.'
      if (v_vmi_like .or. v_bmv_like) then
        v=linintflg(iso%mass, iso%pnt%col(icol1_1)%data, mass(nstars), flag, &
        iso%massflg, iso%pnt%col(icol1_1)%flag)
        if (trim(flag) == '1') flag='outside mass range'
      else if (umb_bmv_like) then
        ! In the colour-colour case get V from another isochrone.
        v=linintflg(iso%mass, iso%pnt%col(icol2_2)%data, mass(nstars), flag, &
        iso%massflg, iso%pnt%col(icol2_2)%flag)
        if (trim(flag) == '1') flag='outside mass range'
      end if
      ! Now, if the flag is not OK, continue, but note the fact.
      ! print*, 'Finding colour.'
      work_flag=flag
      if (v_bmv_like .or. umb_bmv_like) then
        if (icol2_2 == 0) then
          ! The files supplied B-V.
          bmv=linintflg(iso%mass, iso%pnt%col(icol1_2)%data, mass(nstars), flag, &
          iso%massflg, iso%pnt%col(icol1_2)%flag)
        else
          allocate(workflag(iso%npts))
          workflag=iso%pnt%col(icol1_2)%flag
          where(iso%pnt%col(icol2_2)%flag /= 'OK') workflag=iso%pnt%col(icol2_2)%flag
          bmv=linintflg(iso%mass, iso%pnt%col(icol1_2)%data-iso%pnt%col(icol2_2)%data, &
          mass(nstars), flag, iso%massflg, workflag)
          deallocate(workflag)
        end if
        if (trim(flag) == '1') flag='outside mass range'
        vmi=0.0
      else if (v_vmi_like) then
        if (icol2_2 == 0) then
          ! The files supplied V-I.
          vmi=linintflg(iso%mass, iso%pnt%col(icol1_2)%data, mass(nstars), flag, &
          iso%massflg, iso%pnt%col(icol1_2)%flag)
        else
          allocate(workflag(iso%npts))
          workflag=iso%pnt%col(icol1_2)%flag
          where(iso%pnt%col(icol2_2)%flag /= 'OK') workflag=iso%pnt%col(icol2_2)%flag
          vmi=linintflg(iso%mass, iso%pnt%col(icol1_2)%data-iso%pnt%col(icol2_2)%data, &
          mass(nstars), flag, iso%massflg, workflag)
          deallocate(workflag)
        end if
        if (trim(flag) == '1') flag='outside mass range'
        bmv=0.0
      end if
      if (flag /= 'OK') work_flag=flag
      if (umb_bmv_like) then
        ! In the colour-colour case U-B is stored as a magnitude.
        if (icol2_1 == 0) then
          ! The files supplied U-B.
          umb=linintflg(iso%mass, iso%pnt%col(icol1_1)%data, mass(nstars), flag, iso%massflg)
        else
          allocate(workflag(iso%npts))
          workflag=iso%pnt%col(icol1_1)%flag
          where(iso%pnt%col(icol2_1)%flag /= 'OK') workflag=iso%pnt%col(icol2_1)%flag
          umb=linintflg(iso%mass, iso%pnt%col(icol1_1)%data-iso%pnt%col(icol2_1)%data, &
          mass(nstars), flag, iso%massflg, workflag)
          deallocate(workflag)
        end if
        if (trim(flag) == '1') flag='outside mass range'
        if (flag /= 'OK') work_flag=flag
      else
        umb=0.0
      end if
      if (nstars == 1) then
        flag_out=work_flag
      else
        ! This is a secondary star.  
        if (trim(work_flag) /= 'OK') then
          ! Only possible error is that its dropped 
          ! off the bottom of the isochrone.  
          ! Set fluxes already to zero.
          uflux(nstars)=0.0
          bflux(nstars)=0.0
          vflux(nstars)=0.0
          iflux(nstars)=0.0
          cycle each_component
        end if
      end if
      ! And the gravity.
      ! print*, 'Finding gravity.'
      logg=linintflg(iso%mass, iso%logg, mass(nstars), flag, iso%massflg, iso%loggflg)
      !logg=isochr('mas', mass(nstars), 'log', age, iso_file, bc_file, colnam, flag)
      teff=10.0**linintflg(iso%mass, iso%teff, mass(nstars), flag, iso%massflg, iso%teffflg)
      !teff=10.0**isochr('mas', mass(nstars), 'tef', age, iso_file, bc_file, colnam, flag)
      ! For the reddening.
      if (abs(red) > tiny(red)) then
        !print*, 'About to do the reddening.', red, ext_file
        if (len_trim(ext_file) > 0) then
          if (v_bmv_like) then
            red_in1=bmv
            red_in2=v
          else if (v_vmi_like) then
            red_in1=vmi
            red_in2=v
          else if (umb_bmv_like) then
            red_in1=bmv
            red_in2=umb
          else
            print*, 'Hmm, that was unexpected.'
            stop
          end if
          if (bell_there(ext_file)) then
            call bell_ext(ext_file, colnam, red, logg, teff, red_in1, red_in2, extrap)
          else
            call reddening_point(colnam, ext_file, 0.0, red, red_in1, red_in2)
            extrap=.false.
          end if
          if (v_bmv_like) then
            bmv=red_in1
            v=red_in2
          else if (v_vmi_like) then
            vmi=red_in1
            v=red_in2
          else if (umb_bmv_like) then
            bmv=red_in1
            umb=red_in2
          end if
        else
          print*, 'No extinction file specified.'
          stop
        end if
        if (extrap) then
          if (trim(flag_out) == 'OK') flag_out='Extinction'
        end if
      end if
      uflux(nstars)=10.0**(-0.4*(v+bmv+umb))
      bflux(nstars)=10.0**(-0.4*(v+bmv))
      vflux(nstars)=10.0**(-0.4*v)
      iflux(nstars)=10.0**(-0.4*(v-vmi))

    end do each_component

    ! Add up all the components.
    if (v_vmi_like) then
      mag=-2.5*log10(sum(vflux(1:nstars)))
      col=-2.5*log10(sum(vflux(1:nstars))/sum(iflux(1:nstars)))    
    else if (v_bmv_like) then
      mag=-2.5*log10(sum(vflux(1:nstars)))
      col=-2.5*log10(sum(bflux(1:nstars))/sum(vflux(1:nstars)))
    else if (umb_bmv_like) then
      mag=-2.5*log10(sum(uflux(1:nstars))/sum(bflux(1:nstars)))
      col=-2.5*log10(sum(bflux(1:nstars))/sum(vflux(1:nstars)))
      extra_mag=-2.5*log10(sum(vflux(1:nstars)))
    end if
    p_mass=mass(1)


  end subroutine make_star

  real function v_for_ubv()

    v_for_ubv=extra_mag

  end function v_for_ubv

  subroutine iso_image_name(colnam, age, filnam)

    character(len=*), dimension(2), intent(in) :: colnam
    real :: age
    character(len=*), intent(out) :: filnam

    character(len=6) :: age_name

    write(age_name, 30) int(age+tiny(age)), nint(1000*(age-int(age+tiny(age))))
30  format(i2.2, '.', i3.3)

    filnam=trim(colnam(1)) //'_'// &
           trim(colnam(2)) //'_'// &
                age_name          //'.fit'

  end subroutine iso_image_name

  subroutine func2d(a_par, dat, colnam, correlated, &
    distrib, distrib_flag, posterior, rnpts)

    ! The calculates tau^2 with respect to a model.

    use ark_file_io
    use define_star

    ! The parameter array.
    real, intent(in), dimension(:) :: a_par
    ! The data points. The prior membership probabilities are stored as 
    ! colour 3 in dat. 
    type(a_star), dimension(:), intent(in) :: dat
    ! the posterior will be written out in colour 4.
    character(len=*), dimension(2), intent(inout) :: colnam
    ! The correlation flag for s/r likelihood.
    integer, intent(in) :: correlated
    real, intent(out), dimension(:) :: distrib
    character, intent(out), dimension(:) :: distrib_flag
    ! The posterior membership likelihoods.
    real, dimension(:), intent(out) :: posterior
    real, optional, intent(out) :: rnpts

    integer :: i
    real, save :: last_age=-1
    real, dimension(:,:), allocatable, save :: data, axdata
    integer, dimension(2), save :: naxis
    real, dimension(:), allocatable, save :: grad, axgrad
    integer, save :: n_grad
    character(len=256), save :: ext_file, ext_dir
    integer :: iflag

    integer :: ipos, icol, imag
    real :: mag, col
    character(len=50) :: ifname

    real :: sum_gauss

    ! The calculating the nomrmalisation for the model.
    integer :: npix
    real :: mag_min, mag_max
    ! And for the normalisation for the non-member 
    ! model we need the area of the CMD.
    real :: cmd_area

    if (abs(last_age-a_par(1)) > tiny(last_age)) then
      call iso_image_name(colnam, a_par(1), ifname)
      call nxtark_in(ifname)
      iflag=inpark(naxis, data, axdata)
      if (iflag < 0) then
        print*, 'Error reading 2D isochrone file ', ifname
        stop
      end if
      iflag=get_header_s('EXT_FILE', ext_file)
      if (iflag < 0) then
        print*, 'Error reading extinction file name from file ', ifname
        stop
      end if
      ! And add the directory name.
      if (get_header_s('EXT_DIR', ext_dir) >=0) &
      ext_file=trim(ext_dir)//ext_file
      last_age=a_par(1)
    end if

    ! Find the range of magnitudes in the data.
    find_range: do i=1, size(dat,1)
      ! Change to reddening-free absolute magnitude.
      mag=dat(i)%col(1)%data
      col=dat(i)%col(2)%data
      call reddening(colnam, ext_file, -a_par(2), -a_par(3), col, mag)
      if (i ==1) then
        mag_min=mag
        mag_max=mag
      else
        mag_min=min(mag_min, mag)
        mag_max=max(mag_max, mag)
      end if
    end do find_range
    call natural_norm(mag_max, mag_min, data, axdata)

    if (present(rnpts)) rnpts=0.0

    ! Find the area of the CMD covered by the datapoints.
    cmd_area = abs( (maxval(dat%col(1)%data,1)-minval(dat%col(1)%data,1)) & 
                   *(maxval(dat%col(2)%data,1)-minval(dat%col(2)%data,1)) )

    each_point: do i=1, size(dat,1)

      ! Change to reddening-free absolute magnitude.
      mag=dat(i)%col(1)%data
      col=dat(i)%col(2)%data
      call reddening(colnam, ext_file, -a_par(2), -a_par(3), col, mag)

      ! Find the nearest positions in X and Y.  (S/R likelihood will spot
      ! if we are close to the array edge).
      icol=locate_nearest(axdata(1:size(data,1),1), col)
      imag=locate_nearest(axdata(1:size(data,2),2), mag)

      if (present(rnpts)) then
        distrib(i)=likelihood(data, axdata, grad, icol, imag, &
        dat(i)%col(2)%err, dat(i)%col(1)%err, correlated, distrib_flag(i), &
        sum_gauss)
        rnpts=rnpts+sum_gauss
      else
        distrib(i)=likelihood(data, axdata, grad, icol, imag, &
        dat(i)%col(2)%err, dat(i)%col(1)%err, correlated, distrib_flag(i))
      end if

      if (distrib_flag(i) /= 'O') then
        ! Probably outside grid range, so set value to that of the non members.
        distrib(i) = 0.0
      end if

      ! Set the un-normalised posterior probability of membership before distrib
      ! gets overwritten to include the non-member distribution.
      posterior(i)=distrib(i)*dat(i)%col(3)%data

      ! Now adjust distrib to allow for the possibility it's not a member.
      distrib(i) = distrib(i)*dat(i)%col(3)%data + & 
      (1.0-dat(i)%col(3)%data)/cmd_area

      ! And normalise the posterior probability of membership.
      posterior(i)=posterior(i)/distrib(i)

    end do each_point

  end subroutine func2d


    subroutine grid2d(a_par, dat, colnam, correlated, &
      var_par, start, end, &
      data, axdata, best_value, best_par, best_distrib, rnpts)

      use define_star

      ! Performs a 2d grid search.

      ! The model parameters.
      real, dimension(:) :: a_par
      ! The data points.  The third colour should be the prior membership 
      ! probabilities, and the fourth will be filled with the posterior
      ! membership probabilities.
      type(a_star), dimension(:), intent(inout) :: dat
      ! The model and colours being used.
      character(len=*), dimension(2), intent(inout) :: colnam
      ! The correlation flag for s/r likelihood.
      integer, intent(in) :: correlated
      ! The numbers of the parameters to be searched.
      integer, dimension(2), intent(in) :: var_par
      ! The ranges of parameters to be searched.
      real, dimension(2), intent(in) :: start, end
      ! The output grid.
      real, dimension(:,:), intent(out) :: data, axdata
      ! And an optional fraction of the stars on the grid.
      real, dimension(:,:), intent(out), optional :: rnpts

      real, intent(out) :: best_value
      logical :: first_value
      real, intent(out), dimension(:) :: best_par
      real, intent(out), dimension(:) :: best_distrib

      ! Locals.
      character, allocatable, dimension(:) :: distrib_flag
      real, allocatable, dimension(:) :: distrib
      integer :: i1, i2, i3
      real, dimension(size(data,1), size(data,2), size(dat,1)) :: distrib_full
      real :: work, low_tau
      integer, dimension(2) :: loc
      ! The posterior probability that a given star is a member.
      real, dimension(size(data,1), size(data,2), size(dat,1)) :: posterior

      allocate(distrib_flag(size(dat,1)), distrib(size(dat,1)))

      ! Set up the values of best_par to be those of a_par.  We'll
      ! the correct the ones we fit.
      best_par=a_par

      first_value=.true.
      do i1=1, size(data,1)
        if (size(data,1) == 1) then
          a_par(var_par(1))=start(1)
        else
          a_par(var_par(1))=start(1)+(end(1)-start(1))*real(i1-1) &
          /real(size(data,1)-1)
        end if
        axdata(i1,1)=a_par(var_par(1))
        do i2=1, size(data,2)
          if (size(data,2) == 1) then 
            a_par(var_par(2))=start(2)
          else
            a_par(var_par(2))=start(2)+(end(2)-start(2))*real(i2-1) &
            /real(size(data,2)-1)
          end if
          axdata(i2,2)=a_par(var_par(2))
          if (present(rnpts)) then
            call func2d(a_par, dat, colnam, correlated, &
                 distrib_full(i1,i2,:), distrib_flag, posterior(i1, i2, :), &
                 rnpts(i1, i2))
          else 
            call func2d(a_par, dat, colnam, correlated, &
                 distrib_full(i1,i2,:), distrib_flag, posterior(i1, i2, :))
          end if
        end do
      end do

      data=0.0
      high_tau=-2.0*log(2.0*tiny(1.0))*size(dat,1)
      do i1=1, size(data,1)
        do i2=1, size(data,2)
          ! Construct the tau^2 for this grid point by adding
          ! the tau^2 of each point.
          each_point: do i3=1, size(dat,1)
            ! Protect against infinities.
            if (distrib_full(i1,i2,i3) > 2.0*tiny(1.0)) then
              data(i1,i2)=data(i1,i2)-2.0*log(distrib_full(i1,i2,i3))
            else
              ! At least one data point is returning a probablilty so low
              ! we can't cope, so set the entire grid point to a very
              ! high tau^2.
              data(i1,i2)=high_tau
              exit each_point
            end if
          end do each_point
        end do
      end do

      ! Find the best fit.
      loc=minloc(data)
      best_distrib=-2.0*log(distrib_full(loc(1),loc(2),:))
      best_par(var_par(1))=axdata(loc(1),1)
      best_par(var_par(2))=axdata(loc(2),2)
      best_value=sum(best_distrib)
      dat%col(4)%data=posterior(loc(1),loc(2),:)

      deallocate(distrib_flag, distrib)

    end subroutine grid2d

    end module cmdfit_subs
