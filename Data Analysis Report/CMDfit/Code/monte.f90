program monte

  use ark_file_io
  use colteff_subs
  use red_ext_interp
  use cmdfit_subs
  use cmdfit_system

  implicit none

  ! The number of stars in a standard isochrone model.  We use 10 times
  ! this if there is an age spread.
  integer, parameter :: mstar=1000000
  integer :: iflag
  integer :: iage, nage
  real :: age, age_low, age_high, age_step
  integer, dimension(2) :: naxis
  real, allocatable, dimension(:,:) :: data, axdata
  character(len=50) :: ofname, instring
  character(len=8) :: age_namemon
  integer :: i1, i2
  real :: col_start, col_end, mag_start, mag_end, low_mass, high_mass
  integer :: useful_top, useful_bottom

  character(len=200) :: iso_file, bc_file, ext_file, short_ext_file
  character(len=20) , dimension(2) :: colnam

  integer :: i, imag, icol, istar, multiple, icount, start_star, end_star
  real :: p_mass
  ! The number of stars in the simulation, and the usable number.
  integer :: nstar, ustar
  ! The number of stars in the simple isochrone.
  integer, parameter :: niso=10000
  ! The size of the pixels in magnitude.
  real, parameter :: pix_size=0.0025
  ! The binary fraction.
  real :: binary_fraction
  ! More edges.
  real :: mag_min, mag_max, col_min, col_max, col_step, mag_step
  ! A random number.
  real :: harvest
  ! Do we want and age spread?
  logical :: age_spread
  ! A position in a string.
  integer :: ipos

  real, dimension(:), allocatable :: star_mag, star_col
  character(len=50) :: isoflag, out_flag
  logical, dimension(:), allocatable :: isgood
  integer :: iostat
  logical :: extrap
  real :: red, nom_ebmv
  ! The colour the redddening is in.
  character(len=10) :: red_col
  ! Allow extraoplation in colour-Teff relationship?  You may want to do
  ! this for luminous evolved stars as the gravity grids run out.
  ! Its only allowed to extraoplate in gravity, and this currently only works
  ! for the Bessell atmospheres.
  logical, parameter :: extrap_col_teff=.true.
  !call setbug()

  ! This used to be hard coded 0.04 3.0
  print*, 'What range of masses do you want (low, high)?'
  read(*,*) low_mass, high_mass
  
  call choose_models(iso_file, bc_file, ext_file, colnam)

  if (bell_there(ext_file)) then
    print*, 'What E(B-V) you want? (Measured in an energy intergrating system'
    print*, 'with the Bessell filters and the new ODF Kurucz atmospheres.)'
    do
      read(*,*,iostat=iostat) red
      if (iostat == 0) exit
      print*, 'Error in input.'
    end do
    call get_nom_ebmv(red, nom_ebmv, .false., extrap, ext_file)
    red=nom_ebmv
    if (extrap) then
      print*, 'Fell outside range of E(B-V) tables.'
      stop
    end if
    print*, 'Nominal E(B-V) is ', nom_ebmv
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

  print*, 'What binary fraction do you want?'
  read(*,*) binary_fraction

  print*, 'What range of log10(ages) do you want?'
  read(*,*) age_low, age_high
  if (abs(age_low-age_high) > 2.0*tiny(age_low)) then    
    print*, 'What step in age do you want?'
    read(*,*) age_step
    if (age_step < 0.0) then
      age_spread=.true.
      age_step=abs(age_step)
      nage=1+nint((age_high-age_low)/age_step)
      nstar=10*mstar
      ! Round to an even number of stars per age bin.
      nstar=nage*(nstar/nage)
    else
      age_spread=.false.
      nage=1+nint((age_high-age_low)/age_step)
      if (abs(age_high-age_low-real(nage-1)*age_step)/age_step > 0.001) then
         print*, nage, age_high, age_low, age_step, abs(age_high-age_low-real(nage-1)*age_step)/age_step
        print*, 'Age step does not make sense.  You require ', &
        nint(age_high-age_low)/age_step, ' age steps.' 
        stop
      end if
      nstar=mstar
    end if
  else
    nage=1
    nstar=mstar
  end if

  allocate(star_mag(nstar), star_col(nstar), isgood(nstar))

  each_age: do iage=1, nage

    call random_seed()

    if (nage == 1) then
      age = age_low
    else
      age = age_low + real(iage-1)*(age_high-age_low)/real(nage-1)
    end if

    if (.not. age_spread) then

      print*, 'Doing log10(age)', age
      print*, 'With ', nstar, ' simulated objects.'

      ! Write out a simple isochrone file.
      call iso_image_name(colnam, age, ofname)
      ipos=index(ofname, '.fit', .true.)
      open (unit=12, file=ofname(1:ipos)//'iso')
      write(12,*) '# log10(age) = ', age
      write(12,*) '# Interiors from '//trim(iso_file)
      write(12,*) '# Bolometric corrections from '//trim(bc_file)
      write(12,*) '# Reddening is E('//trim(red_col)//')=', red
      write(12,*) '# ', trim(colnam(2)), ' ', trim(colnam(1)), ' primary_mass flag isgood'

    else

      call iso_image_name(colnam, 0.0, ofname)

    end if

    if (age_spread) then
      start_star = 1 + (iage-1)*nstar/nage
      end_star = start_star + nstar/nage - 1
    else
      start_star=1
      end_star=nstar
    end if

    create_each_star: do istar=start_star, end_star

      ! print*, istar, start_star, end_star

      if (istar - 100000*nint(real(istar)/100000) == 0) &
      print*, 'Created ', 100.0*real(istar)/real(nstar), &
      ' percent of stars.'!, age

      call make_star(age, iso_file, bc_file, ext_file, colnam, 'kroupa', &
      low_mass, high_mass, &
      0.0, binary_fraction, red, star_mag(istar), star_col(istar), &
      p_mass, multiple, isoflag)
      isgood(istar)=.false.
      if (trim(isoflag)=='OK') isgood(istar)=.true.
      if (trim(isoflag)=='Colour-Teff (grav)' .and. extrap_col_teff) &
      isgood(istar)=.true.
      ! if (isgood(istar) .and. istar < niso) write(45,*) p_mass, multiple
      out_flag=isoflag
      if (index(trim(out_flag), ' ') /= 0) out_flag='"'//trim(out_flag)//'"'
      ! Not sure why this block was ever in.  It has the effect of writing out the single stars, which
      ! when combined with the block below meant single stars were written out twice.
      !if (istar<=niso .and. multiple==1 .and. .not. age_spread) then
      !  write(12,*) star_col(istar), star_mag(istar), &
      !  p_mass, out_flag, isgood(istar)
      !end if
      if (istar<=niso .and. .not. age_spread) then
        write(12,*) star_col(istar), star_mag(istar), &
        p_mass, out_flag, isgood(istar)
      end if
      !print*, istar, star_mag(istar), star_col(istar), isgood(istar)

    end do create_each_star

    if (.not. age_spread) close(12)

    if (end_star==nstar .or. (.not. age_spread)) then

      ! Pack the arrays.
      icount=0
      do istar=1, nstar
        if (isgood(istar)) then
          icount=icount+1
          star_mag(icount)=star_mag(istar)
          star_col(icount)=star_col(istar)
        end if
      end do
      ustar=icount
      print*, 'There were ', ustar, ' usable stars.'

      ! Find the size of image needed.    
      mag_start=maxval(star_mag(1:ustar))
      mag_end  =minval(star_mag(1:ustar))
      ! Its often best to be a little wider in colour.
      col_start=minval(star_col(1:ustar))-0.2
      col_end  =maxval(star_col(1:ustar))+0.2
      ! Round these to "nice" values.  This makes the grids look nicer, since
      ! the aliasing effects are minimised.
      mag_start=real(int(mag_start/pix_size)  )*pix_size
      mag_end  =real(int(mag_end  /pix_size)+1)*pix_size
      col_start=real(int(col_start/pix_size)+1)*pix_size
      col_end  =real(int(col_end  /pix_size)  )*pix_size
      naxis(1)=nint((col_end-col_start)/pix_size)+1
      naxis(2)=nint((mag_start-mag_end)/pix_size)+1
      print*, 'Grid will be ', naxis
      allocate(data(naxis(1),naxis(2)), axdata(maxval(naxis),2))
      axdata=0.0
      data=0.0
      print*, 'Covering mags ', mag_start, ' to ', mag_end
      print*, 'and colours ', col_start, ' to ', col_end
      col_step=(col_end-col_start)/real(naxis(1)-1)
      do i=1, naxis(1)
        axdata(i,1)=col_start+col_step*real(i-1)
      end do
      mag_step=(mag_end-mag_start)/real(naxis(2)-1)
      do i=1, naxis(2)
        axdata(i,2)=mag_start+mag_step*real(i-1)
      end do

      place_star: do istar=1, ustar

        if (istar - 100000*nint(real(istar)/100000) == 0) &
        print*, 'Placed ',  100.0*real(istar)/real(ustar), &
        ' percent of stars in grid.'

        ! Find the pixel number.
        imag=int(0.5+(star_mag(istar)-mag_start)/mag_step)+1
        if (imag < 1) cycle place_star
        if (imag > naxis(2)) cycle place_star
        icol=int(0.5+(star_col(istar)-col_start)/col_step)+1
        if (icol < 1) cycle place_star
        if (icol > naxis(1)) cycle place_star

        ! Put the star into the image.
        data(icol,imag)=data(icol,imag)+1.0

      end do place_star

      ! Make the model integrate to one over all area.
      data=data/(pix_size*pix_size*sum(data))

      ! First write out the image file.
      ipos=index(ext_file, trim(data_dir()), .true.)
      if (ipos > 0) then
        short_ext_file=ext_file(len_trim(data_dir())+1:len_trim(ext_file))
        call put_header_s('EXT_DIR', trim(data_dir()), 'Extinction directory', 1)
      else
        short_ext_file=trim(ext_file)
      end if
      call put_header_s('EXT_FILE', trim(short_ext_file), 'Extinction file',1)
      if (red > tiny(red)) call put_header_r('RED', red, 'Reddening in E('//trim(red_col)//')',1)
      call put_header_r('BINFRAC', binary_fraction, 'Binary star fraction' , 1)
      call put_header_s('CTYPE1', trim(colnam(2)), ' ', 1)
      call put_header_s('CTYPE2', trim(colnam(1)), ' ', 1)
      if (age_spread) then
        call put_header_r('MIN_AGE', age_low, 'Youngest age used', 1)
        call put_header_r('MAX_AGE', age_high, 'Oldest age used', 1)
      end if
      call nxtark_out(ofname)
      iflag=makark(naxis, data, axdata)
      deallocate(data, axdata)

    end if

  end do each_age

end program monte
