program grid

  use define_star
  use cmdfit_subs
  use ark_file_io
  use colteff_subs

  implicit none

  integer, parameter :: n_par=3
  type(a_star), dimension(:), allocatable :: star
  real, dimension(n_par) :: a_par, best_par
  integer :: iflag, n_term, icol
  integer :: iostat, i
  integer :: nstars, mstars, istar
  integer, dimension(2) :: naxis
  integer, dimension(2) :: var_par
  real, dimension(2) :: start, end
  real, dimension(:,:), allocatable :: data, axdata, rnpts
  character(len=80), dimension(2) :: ctype
  real :: best_value, age_step
  type(a_star) :: test
  integer :: icol1, icol2
  integer :: ncol
  character(len=15), dimension(mcol) :: colstr
  character(len=50) :: ifname, fitsname
  integer :: iunit
  real, dimension(:), allocatable :: best_distrib

  type(a_star), dimension(:), allocatable :: dat
  character(len=10) :: colnam(2)
  character(len=256) :: ext_file, ext_dir
  real :: red_init

  real, allocatable, dimension(:) :: sorted_distrib
  real :: clipped_tau
  real :: prior_frac, prior_max
  integer :: iprior

  integer :: iminus, iband, correlated

  real :: mag_shift, colour_shift, global_mag_shift, pix_size, &
  global_colour_shift
  integer :: i_mag_shift
  real :: sys_mag, sys_col
  real :: low_tau

  ! Mode is 1 for searching in age and/or distance modulus.
  ! Anything else for searching in distance modulus and/or extinction.
  real, parameter :: mode=1

  call setbug()

  print*, 'First read the input catalogue.'
  iflag=read_cluster_file(mstars, ncol, colstr, star)

  print*, '* Will fit magnitude ', colstr(1)
  if (ncol == 2) then
    print*, '* Against colour ', colstr(2)
    icol1=1
    icol2=2
  else
    print*, '* Available magnitudes are '
    do icol=1, ncol
      print*, '* ', icol, colstr(icol) 
    end do
    print*, '> Give the number of the one you want.'
    read(*,*) icol1
    print*, '* Available colours are '
    do icol=1, ncol
      print*, '* ', icol, colstr(icol) 
    end do
    print*, '> Give the number of the one you want.'
    read(*,*) icol2
  end if
     
  nstars=mstars
  do istar=mstars, 1, -1
    if (star(istar)%col(icol1)%flg/='OO' .or. &
      star(istar)%col(icol2)%flg /= 'OO') then
      if (istar /= mstars) star(istar:mstars-1)=star(istar+1:mstars)
      nstars=nstars-1
    end if
  end do
  print*, 'Read ', nstars, ' of ', trim(colstr(icol1)), ' vs ', &
  trim(colstr(icol2))

  ! Now lets sort out whether you add or subtract the colour from the
  ! magnitude to get the other magnitude.
  if (index(colstr(icol1), '-') == 0) then
    ! This is truly a magnitude.
    iminus=index(colstr(icol2), '-')
    if (iminus == 0) then
      print*, 'Cannot find the minus sign in colour ', colstr(icol2)
      stop
    end if
    iband=index(colstr(icol2), trim(colstr(icol1)))
    if (iband == 0) then
      print*, 'Cannot find the magnitude ', trim(colstr(icol2)), &
      ' in the colour ', trim(colstr(icol2))
    else if (iband > iminus) then
      print*, 'Will add ', trim(colstr(icol1)), ' to ', trim(colstr(icol2)), &
      ' to create ', colstr(icol2)(1:iminus-1)
      correlated = 1
    else
      print*, 'Will subtract ', trim(colstr(icol2)), ' from ', &
      trim(colstr(icol1)), &
      ' to create ', trim(colstr(icol2)(iminus+1:len(colstr(icol2))))
      correlated = -1
    end if
  else
    ! Its colour-colour data.
    correlated=0
  end if

  print*, '> Give the extra uncertainty to be added to each magnitude.'
  read(*,*) sys_mag
  print*, '> Give the extra uncertainty to be added to each colour.'
  read(*,*) sys_col


  if (mode == 1) then
    ! Age.
    var_par(1)=1
    print*, '> Give range of log10(ages) to be searched.'
    read(*,*) start(1), end(1)
    if (abs(start(1)-end(1)) > 2.0*tiny(end(1))) then
      print*, '> And the log10(age) step to be used.'
      read(*,*) age_step
      naxis(1)=1+nint((end(1)-start(1))/age_step)
      if (abs(end(1)-start(1)-real(naxis(1)-1)*age_step)/age_step > 0.001) then
        print*, 'Age step does not make sense.  You require ', &
        nint(end(1)-start(1))/age_step, ' age steps.' 
        print*, end(1)-start(1), real(naxis(1)-1)*age_step
        stop
      end if
    else
      naxis(1)=1
    end if
  else 
    ! The extinction in the colour in question.
    var_par(1)=3
    print*, '> Give the range of reddening to search in ', &
    trim(colnam(2))
    read(*,*) start(1), end(1)
    print*, '> And the number of reddening points in the grid.'
    read(*,*) naxis(1)
  end if
  ! The distance modulus.
  print*, '> Give the range in distance modulus to be searched.'
  read(*,*) start(2), end(2)
  print*, '> And the number of distance modulus points in the grid.'
  read(*,*) naxis(2)
  var_par(2)=2

  colnam(1)=colstr(icol1)
  colnam(2)=colstr(icol2)
  iprior=-1
  do icol=3, ncol
    if (trim(colstr(icol)) == 'Prior') iprior = icol
  end do

  if (mode == 1) then

     ! read the first model to see if there is an extinction.
     call iso_image_name(colnam, start(1), fitsname)
     iflag=0
     iunit=0
     call fithed(fitsname, .true., iunit, iflag)
     close(iunit)
     if (iflag < 0) then
        print*, 'Failed to find file ', trim(ifname)
        stop
     end if
     iflag=get_header_r('RED', red_init)
     if (iflag > 0) then
       print*, '  The models have already been reddened by ', red_init, ' mags.'
       print*, '> Give extra reddening you wish to apply (or zero for no more).'
     else
        red_init=0.0
        ! The extinction in the colour in question.
        print*, '> Give the required reddening in ', trim(colnam(2))
     end if
     read(*,*) a_par(3)

  else
    ! Age
    print*, '> Give the log10(age).'
    read(*,*) a_par(1)
  end if

  if (iprior == -1) then
    print*, '> Give the fraction of stars which are likely to be members.'
    read(*,*) prior_frac
  else
    print*, '> Give the maximum probability for membership.'
    read(*,*) prior_max
    ! print*, 'Using individual prior probabilties of membership from the file.'
  end if

  allocate(data(naxis(1),naxis(2)), axdata(maxval(naxis),2), &
  rnpts(naxis(1),naxis(2)))

  n_term=2

  allocate(dat(nstars))

  do i=1, nstars
    ! Increase the uncertainties by the systematic.
    star(i)%col(icol1)%err  = sqrt(star(i)%col(icol1)%err**2.0 + sys_mag**2.0)
    ! Bug corrected here August 11 2006. 
    ! Used to be...
    ! star(i)%col(2)%err  = sqrt(star(i)%col(icol2)%err**2.0 + 2.0*(sys**2.0))
    ! The changed in the final stages of the NGC2169 paper (1st September 2006)
    ! to allow different uncertainties in mag and colour.
    star(i)%col(icol2)%err  = sqrt(star(i)%col(icol2)%err**2.0 + sys_col**2.0)
    ! Copy the star structure accross, which deals with colour 1.
    dat(i)=star(i)
    ! And now make sure colours are right.
    dat(i)%col(1)%data = star(i)%col(icol1)%data
    dat(i)%col(1)%err  = star(i)%col(icol1)%err
    colstr(1)=colstr(icol1)
    dat(i)%col(2)%data = star(i)%col(icol2)%data
    dat(i)%col(2)%err  = star(i)%col(icol2)%err
    colstr(2)=colstr(icol2)
    ! And set up the prior (colour 3) and posterior (colour 4) 
    ! probabilties of membership.
    if (iprior == -1) then
      dat(i)%col(3)%data = prior_frac
    else
      dat(i)%col(3)%data = star(i)%col(iprior)%data*prior_max
    end if
    dat(i)%col(3)%err  = 0.0
    dat(i)%col(3)%flg  = 'OO'
    colstr(3)='Prior'
    dat(i)%col(4)%data = 1.0
    dat(i)%col(4)%err  = 0.0
    dat(i)%col(4)%flg  = 'OO'
    colstr(4)='Post'
  end do

  if (sum(dat%col(1)%err) >= 0.99*sum(dat%col(2)%err)) then
    print*, 'The uncertanties in magnitude are larger than in colour so'
    print*, 'I will take them as uncorrelated.'
    correlated = 0
  end if

  allocate(best_distrib(nstars))
  call grid2d(a_par, dat, colnam, correlated, &
  var_par, start, end, &
  data, axdata, best_value, best_par, best_distrib, rnpts)

  low_tau=minval(best_distrib)

  ! Write out the grid.
  call nxtark_out('grid.fit')
  ! Force it to be a binary file as uncer wants the headers.
  call typark(1)
  call put_header_r('HI_TAU', high_calc_tau(), 'HIGHEST CALCULABLE TAU^2', 1)
  ctype='NOT SET'
  if (naxis(1) > 1) then
    if (mode == 1) then
      ctype(1)='LOG10(AGE)'
    else
      ctype(1)='EXTINCTION'
    end if
    call put_header_s('CTYPE1', ctype(1), ' ', 1)
    if (naxis(2) > 1) then
      ctype(2)='DISTANCE MODULUS'
      call put_header_s('CTYPE2', ctype(2), ' ', 1)
      i=makark(naxis, data, axdata)
      call nxtark_out('grid_npts.fit')
      i=makark(naxis, rnpts, axdata)
    else
      i=makark(naxis, data, axdata)
      ! And an ASCII file.
      open(unit=11, file='grid.asc')
      write(11,'(/,/)')
      do i=1, naxis(1)
        write(11,*) axdata(i,1), data(i,1)
      end do
      close(11)
    end if
  else
    ! This has to be called ctype2 for the paragraph of code below.
    ctype(2)= 'DISTANCE MODULUS'
    call put_header_s('CTYPE1', ctype(2), ' ', 1)
    i=makark(naxis(2), data(1,:), axdata(:,2))
    ! And an ASCII file.
    open(unit=11, file='grid.asc')
    write(11,'(/,/)')
    do i=1, naxis(2)
      write(11,*) axdata(i,2), data(1,i)
    end do
    close(11)
    call nxtark_out('grid_npts.fit')
    i=makark(naxis(2), rnpts(1,:), axdata(:,2))
    ! And an ASCII file.
    open(unit=11, file='grid_npts.asc')
    write(11,'(/,/)')
    do i=1, naxis(2)
      write(11,*) axdata(i,2), rnpts(1,i)
    end do
    close(11)
  end if


  print*, 'Lowest value ', best_value
  do i=1, 2
    if (ctype(i) == 'NOT SET') cycle
    print*, ' '
    if (ctype(i)(1:5) == 'LOG10') then
      print*, ' For parameter ', trim(ctype(i)(6:len(ctype)))//'/10^6'
      print*, ' Range searched ', 10.0**(axdata(1,i)-6.0), &
      10.0**(axdata(naxis(i),i)-6.0)
      print*, ' Best value ', 10.0**(best_par(var_par(i))-6.0)
    else
      print*, ' For parameter ', trim(ctype(i))
      print*, ' Range searched ', axdata(1,i), axdata(naxis(i),i)
      print*, ' Best value ', best_par(var_par(i))
    end if
  end do

  call iso_image_name(colnam, best_par(1), ifname)
  call nxtark_in(ifname)
  iflag=inpark(naxis, data, axdata)
  ! Get the name of the extinction file.
  iflag=get_header_s('EXT_FILE', ext_file)
  if (iflag < 0) then
    print*, 'Error reading extinction file name from file ', ifname
    stop
  end if
  ! And add the directory name.
  if (get_header_s('EXT_DIR', ext_dir) >=0) ext_file=trim(ext_dir)//ext_file

  ! Write a file of the tau^2s
  open (unit=66, file='tau_point.out')
  write(66,*) 'Columns are field, id, colour, magnitude and tau^2'
  write(66,'(/)')
  do i=1, nstars
    write(66,*) dat(i)%field, dat(i)%id, dat(i)%col(2)%data, &
    dat(i)%col(icol1)%data, best_distrib(i)
  end do
  close(66)

  allocate(sorted_distrib(nstars))
  sorted_distrib=best_distrib

  ! Now sort them.
  call tau_sort(sorted_distrib)

  ! Write out the points in apparent magnitude.
  call nxtcls_out('fitted.cat')
  call write_cluster_file(4, colstr, dat, nstars)

  ! Change to absolute magnitude and write out.
  do i=1, nstars
    call reddening(colnam, ext_file, -best_par(2), -best_par(3), &
    dat(i)%col(2)%data, dat(i)%col(1)%data)
  end do
  colstr(icol1)='('//trim(colstr(icol1))//')o'
  colstr(icol2)='('//trim(colstr(icol2))//')o'
  call nxtcls_out('fitted_abs.cat')
  call write_cluster_file(4, colstr, dat, nstars)
  close(32)

  ! Change to absolute magnitude and write out.
  ! This is complex as the magnitude shift is colour dependent.
  ! The following only does the job to the nearest half pixel.
  ! Find the magnitude shift at the middle pixel.
  global_colour_shift=(axdata(1,1)+axdata(naxis(1),1))/2.0
  global_mag_shift=0.0
  call reddening(colnam, ext_file, best_par(2), best_par(3), &
  global_colour_shift, global_mag_shift)
  ! Apply it.
  axdata(1:naxis(2),2)=axdata(1:naxis(2),2)+global_mag_shift
  do i=1, naxis(1)
    ! For each colour, work out the shift from nominal.
    colour_shift=axdata(i,1)
    mag_shift=0.0
    call reddening(colnam, ext_file, best_par(2), best_par(3), &
    colour_shift, mag_shift)
    mag_shift=mag_shift-global_mag_shift
    pix_size=(axdata(naxis(2),2)-axdata(1,2))/real(naxis(2)-1)
    i_mag_shift=nint(mag_shift/pix_size)
    if (i_mag_shift > 0) then
      data(i,1+i_mag_shift:naxis(2))=data(i,1:naxis(2)-i_mag_shift)
      data(i,1:i_mag_shift)=0.0
    else if (i_mag_shift < 0) then
      data(i,1:naxis(2)+i_mag_shift)=data(i,1-i_mag_shift:naxis(2))
      data(i,naxis(2)+i_mag_shift+1:naxis(2))=0.0
    end if
  end do
  ! Apply the colour shift..
  global_colour_shift=global_colour_shift-(axdata(1,1)+axdata(naxis(1),1))/2.0
  axdata(1:naxis(1),1)=axdata(1:naxis(1),1)+global_colour_shift

  call nxtark_out('best_model.fit')
  iflag=makark(naxis, data, axdata)

  open(unit=22, file='distrib.tau')
  write(22, '(/,/)')
  do i=1, nstars
    write(22,*) sorted_distrib(i), real(nstars-i+1)/real(nstars)
    write(22,*) sorted_distrib(i), real(nstars-i)/real(nstars)
  end do
  close(22)

  print*, ' '
  print*, 'Prior membership fraction ', sum(dat%col(3)%data)/real(nstars)
  print*, 'Posterior membership fraction is ', &
  sum(dat%col(4)%data)/real(nstars)
  if (iprior == -1) then
    print*, 'Difference and Prior_frac are ', &
    (sum(dat%col(3)%data)-sum(dat%col(4)%data))/real(nstars), prior_frac
  else
    print*, 'Difference and Prior_max are ', &
    (sum(dat%col(3)%data)-sum(dat%col(4)%data))/real(nstars), prior_max
  end if

end program grid
