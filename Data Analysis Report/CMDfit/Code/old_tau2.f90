module tau_subs

  implicit none

  contains

  subroutine prob_bin(chi2, prob, num, chi2_of_bin, prob_in_bin)

    ! Given a set of chi-squareds, this routine will work out a binned
    ! probability distribution of chi-squared.

    ! The input probability points.
    real, dimension(:), intent(inout) :: chi2, prob 
    ! The number of bins to be used.
    integer, intent(in) :: num
    real, dimension(num), intent(out) :: prob_in_bin, chi2_of_bin
    real, dimension(num) :: weight

    real :: chi2_min, chi2_max
    integer :: i, itest, num_of_chi2, i_val
    real, dimension(num) :: bound_low, bound_up
    real :: work

    ! Do the binning in a log space.
    work=minval(chi2)-1.0
    chi2=log(chi2-work)

    ! Set up some histogram boundaries.
    chi2_min=minval(chi2)
    chi2_max=maxval(chi2)
    ! print*, 'Corresponding chi-squared values are ', chi2_min, chi2_max 

    do i=1, num
      bound_low(i)=chi2_min+real(i-1)*(chi2_max-chi2_min)/real(num)
      bound_up (i)=chi2_min+real(i  )*(chi2_max-chi2_min)/real(num)
    end do

    chi2_of_bin=0.0
    prob_in_bin=0.0
    num_of_chi2=0
    weight=0.0
    itest=0
    each_value: do i_val=1, size(chi2)
      do i=1, num
        if (chi2(i_val)>bound_low(i) .and. chi2(i_val)<=bound_up(i)) then
          prob_in_bin(i)=prob_in_bin(i)+prob(i_val)
          chi2_of_bin(i)=chi2_of_bin(i)+chi2(i_val)*prob(i_val)
          itest=itest+1
          cycle each_value
        end if
      end do
      if (chi2(i_val) <= bound_low(1)) then
        prob_in_bin(1)=prob_in_bin(1)+prob(i_val)
        chi2_of_bin(1)=chi2_of_bin(1)+chi2(i_val)*prob(i_val)
        itest=itest+1
        cycle each_value
      end if
      if (chi2(i_val) > bound_up(num)) then
        prob_in_bin(num)=prob_in_bin(num)+prob(i_val)
        chi2_of_bin(num)=chi2_of_bin(num)+chi2(i_val)*prob(i_val)
        itest=itest+1
        cycle each_value
      end if
      print*, 'Problem in prob_bin', bound_low(1), chi2(i_val), bound_up(num), i_val
      print*, i_val
      stop
    end do each_value

    where(chi2_of_bin > 0.0) 
      chi2_of_bin=chi2_of_bin/prob_in_bin
    elsewhere
      chi2_of_bin=(bound_up+bound_low)/2.0
    endwhere

    ! Back into the linear space.
    chi2_of_bin=exp(chi2_of_bin) + work

  end subroutine prob_bin

end module tau_subs


program tau

  use ark_file_io
  use likelihood_mod
  use quad
  use define_star
  use tau_subs
  use cmdfit_system

  implicit none

  integer :: ipos, i1, i2, i, j
  real, dimension(:,:), allocatable, save :: data, axdata, output
  integer, dimension(2) :: naxis
  real, dimension(:), allocatable, save :: grad, axgrad
  integer, save :: n_grad
  integer :: iflag
  character :: flag
  character(len=10) :: sflag
  character(len=50) :: ifname

  ! For the cluster file.
  integer :: nstars, ncol, istar, jstar
  character(len=10), dimension(iso_mcol) :: colstr
  type(a_star), dimension(:), allocatable :: star
  character(len=30) :: clus_fil

  ! Integers for colour and magnitude in pixel space.
  integer :: icol, imag, i_mag_sig, i2_low, i2_high

  ! To correct for the number of degrees of freedom.
  integer :: nparams
  
  ! The number of bins.
  integer, parameter :: num=1000
  real, dimension(num) :: tot_chi2_bin, tot_prob_bin
  real, dimension(num) :: star_prob_bin, star_chi2_bin
  integer :: num_tot, tot_start, tot_end
  real, dimension(num*num) :: prod_prob, prod_chi2
  real, allocatable, dimension(:) :: prob, chi2

  integer :: n_val, i_val
  character(len=50) :: str_range
  real :: low, high
  integer :: jcol1, jcol2
  integer :: iminus, iband, correlated

  ! Distribution of tau for one realisation.
  integer, parameter :: n_one=10000
  real, dimension(n_one) :: tau_one, frac_one
  integer :: near

  real :: expect, sigma
  real, allocatable, dimension(:) :: integ_x, integ_y

  call setbug()

  print*, 'First read in the best fitting model FITS file.'
  iflag=inpark(naxis, data, axdata)
  if (iflag < 0) then
    print*, 'Error reading 2D isochrone file ', ifname
    stop
  end if

  allocate(output(naxis(1), naxis(2)))

  ! Now read the gradient file.
  call lstark_in(ifname)
  ipos=index(ifname, '.fit')
  if (ipos > 0) then
    ifname(ipos:ipos+3)='.grd'
    call nxtark_in(ifname)
    iflag=inpark(n_grad, grad, axgrad)
    if (iflag == -2) then
      print*, 'No gradient file found, using natural normalisation.'
    else if (iflag < 0) then
      print*, 'Error reading gradient file file ', ifname
      stop
    end if
  else
    print*, 'Not attempting to find gradient file, using natural normalisation.'
  end if

  print*, '> Give the file of data points (or type end for chunk of image).'
  read(*,'(a)') clus_fil

  if (clus_fil == 'end') then

    print*, '> Give the number of datapoints and free parameters.'
    read(*,*) nstars, nparams

    allocate(star(nstars))
    jcol1=1
    jcol2=2

    print*, '> Give uncertainty in magnitude.'
    read(*,*) star(1)%col(jcol1)%err
    star%col(jcol1)%err=star(1)%col(jcol1)%err

    star%col(jcol2)%err=sqrt(2.0)*star%col(jcol1)%err
    print*, 'Assuming uncertainty in colour of ', star(1)%col(jcol2)%err

    print*, 'And that the uncertainties are uncorrelated.'
    correlated = 0

  else

    call nxtcls_in(clus_fil)
    iflag=read_cluster_file(nstars, ncol, colstr, star)

    if (ncol > 2) then
      print*, 'The available magnitudes are'
      do jcol1=1, ncol
        print*, jcol1, colstr(jcol1)
      end do
      print*, 'Give the number of the magnitude you want.'
      read(*,*) jcol1
      print*, 'The available colours are'
      do jcol2=1, ncol
        print*, jcol2, colstr(jcol2)
      end do
      print*, 'Give the number of the colour you want.'
      read(*,*) jcol2
    else
      jcol1=1
      jcol2=2
    end if

    ! Check for flagged stars.
    istar=1
    do 
      if (istar > nstars) then
        exit
      else if (istar == nstars) then
        if (star(istar)%col(jcol1)%flg/='OO' &
        .or. star(istar)%col(jcol2)%flg/='OO') nstars=nstars-1
        exit
      else
        if (star(istar)%col(jcol1)%flg/='OO' &
          .or. star(istar)%col(jcol2)%flg/='OO') then
          print*, 'Removing ', star(istar)%id, star(istar)%col(jcol1)%flg, &
          star(istar)%col(jcol2)%flg
          star(istar:nstars-1)=star(istar+1:nstars)
          nstars=nstars-1
        else
          istar=istar+1
        end if
      end if
    end do

    !do istar=1, nstars
    !  print*, star(istar)%id, star(istar)%col(jcol1)%data, &
    !  star(istar)%col(jcol2)%data
    !end do
    !print*, nstars

    ! Now lets sort out whether you add or subtract the colour from the
    ! magnitude to get the other magnitude.
    iminus=index(colstr(jcol2), '-')
    if (iminus == 0) then
      print*, 'Cannot find the minus sign in colour ', colstr(jcol2)
      stop
    end if
    iband=index(colstr(jcol2), trim(colstr(jcol1)))
    if (iband == 0) then
      print*, 'Cannot find the magnitude ', trim(colstr(jcol2)), &
      ' in the colour ', trim(colstr(jcol2))
      print*, 'Assuming colours are uncorrelated.'
      correlated=0
    else if (iband > iminus) then
      print*, 'Will add ', trim(colstr(jcol1)), ' to ', trim(colstr(jcol2)), &
      ' to create ', colstr(jcol2)(1:iminus-1)
      correlated = 1
    else
      print*, 'Will subtract ', trim(colstr(jcol2)), ' from ', &
      trim(colstr(jcol1)), &
      ' to create ', trim(colstr(jcol2)(iminus+1:len(colstr(jcol2))))
      correlated = -1
    end if
    if (sum(star%col(jcol1)%err) >= 0.99*sum(star%col(jcol2)%err)) then
      print*, 'The uncertanties in magnitude are larger than in colour so'
      print*, 'I will take them as uncorrelated.'
      correlated = 0
    end if

    print*, '> Give the number of free parameters.'
    read(*,*) nparams

    ! Normalise the data.
    call natural_norm(maxval(star(1:nstars)%col(jcol1)%data,1), &
    minval(star(1:nstars)%col(jcol1)%data,1), data, axdata)    

  end if

  ! Set up the single realisation tau (x) axis.
  do i=1, n_one
    tau_one(i)=100.0*(real(i-n_one/2)-0.5)/real(n_one)
  end do
  ! And empty the fraction axis.
  frac_one=0.0

  ! Now, for each star
  each_star: do istar=1, nstars

    if (clus_fil == 'end') then

      print*, '> Give range of magnitudes to be used (<cr>=all).'
      read(*,'(a)') str_range

      if (str_range == ' ') then
        i2_low=1
        i2_high=naxis(2)
      else

        read(str_range,*) low, high
        i2_low =minloc(abs(axdata(:,2)-low ),1)
        i2_high=minloc(abs(axdata(:,2)-high),1)

        if (i2_low > i2_high) then
          iflag=i2_low
          i2_low=i2_high
          i2_high=iflag
        end if

      end if

      print*, 'Range of pixels is ', i2_low, i2_high

    else

      print*, 'Doing star ', star(istar)%id

      ! Find the nearest pixel in colour-magnitude space.
      imag=minloc(abs(axdata(1:naxis(2),2)-star(istar)%col(jcol1)%data),1)
      icol=minloc(abs(axdata(1:naxis(1),1)-star(istar)%col(jcol2)%data),1)

      ! Calculate how many pixels the magnitude uncertainty represents.
      i_mag_sig=nint(star(istar)%col(jcol1)%err/&
      abs(axdata(imag,2)-axdata(imag+1,2)))

      i2_low =max(imag-3*i_mag_sig,1)
      i2_high=min(imag+3*i_mag_sig,naxis(2))
      if (i2_low<1 .or. i2_high>naxis(2)) then
        print*, 'Warning, this star falls outside the image'
        i2_low=max(1, i2_low)
        i2_high=min(i2_high, naxis(2))
      end if
    end if

    if (star(istar)%col(jcol1)%err > 0.0) then
      output=0.0
      do i2=i2_low, i2_high
        do i1=1, naxis(1) 
          output(i1,i2)=likelihood(data, axdata, grad, i1, i2, &
          star(istar)%col(jcol2)%err, star(istar)%col(jcol1)%err, correlated, &
          flag)
        end do
      end do
    else
      ! The error function is a delta function.
      output=data
    end if
    ! Make a file of the output data to look at it.
    !iflag=makark(naxis, output, axdata)
    !stop

    n_val=0
    do i1=1, naxis(1)
      do i2=i2_low, i2_high
        if (output(i1,i2) > tiny(output(i1,i2))) n_val=n_val+1
      end do
    end do

    print*, 'Found ', n_val, ' pixels above zero.'

    ! Put the 2D array into a 1D array, missing out the zero probabilty
    ! datapoints.
    allocate(prob(n_val), chi2(n_val))
    i_val=0
    prob=0.0
    chi2=0.0
    do i1=1, naxis(1)
      do i2=i2_low, i2_high
        if (output(i1,i2) > tiny(output(i1,i2))) then
          i_val=i_val+1
          prob(i_val)=output(i1,i2)
          chi2(i_val)=-2.0*log(output(i1,i2))
        end if
      end do
    end do

    if (i_val /= n_val) then
      print*, 'Problem in tau.'
      stop
    end if

    prob=prob/sum(prob)
    do j=1, n_val
      near=minloc(abs(chi2(j)-tau_one),1)
      frac_one(near)=frac_one(near)+prob(j)
    end do

    call prob_bin(chi2, prob, num, star_chi2_bin, &
    star_prob_bin)

    deallocate(prob, chi2)

    if (istar == 1) then
      tot_prob_bin=star_prob_bin
      tot_chi2_bin=star_chi2_bin

    else

      ! Normalise, to stop numerical overflow.
      !tot_prob_bin=1.0e22*tot_prob_bin/maxval(tot_prob_bin)
      tot_prob_bin=tot_prob_bin/maxval(tot_prob_bin)

      ! Find the number of bins with zero in at the end of the array.
      low_empty_bins: do j=1, num
        if (tot_prob_bin(j) >= tiny(tot_prob_bin(j))) then
          tot_start=j
          exit low_empty_bins
        end if
      end do low_empty_bins

      ! Find the number of bins with zero in at the end of the array.
      num_tot=num
      find_empty_bins: do j=num, 1, -1
        if (tot_prob_bin(j) >= tiny(tot_prob_bin(j))) then
          tot_end=j
          exit find_empty_bins
        end if
      end do find_empty_bins

      num_tot=tot_end-tot_start+1

      ! Create a product array.
      ! rewind(66)
      do i=1, num
        do j=tot_start, tot_end
          prod_prob((i-1)*num_tot + j-tot_start+1)= &
          tot_prob_bin(j)*star_prob_bin(i)
          ! write(66,*) (i-1)*num_tot + j-tot_start+1, j, tot_chi2_bin(j), i, star_chi2_bin(i)
          prod_chi2((i-1)*num_tot + j-tot_start+1)= &
          tot_chi2_bin(j)+star_chi2_bin(i)
        end do
      end do

      !call sort(prod_chi2, prod_prob, num*num)

      call prob_bin(prod_chi2(1:num*num_tot), prod_prob(1:num*num_tot), num, &
      tot_chi2_bin, tot_prob_bin)

    end if

    print*, 'Modal tau is ', tot_chi2_bin(maxloc(tot_prob_bin)), &
    ' for star ', istar, ' out of ', nstars

    ! if (clus_fil == 'end') iflag=makark(naxis, output, axdata)

  end do each_star

  allocate (integ_x(1:num-1), integ_y(1:num-1))
  do i=1, num-1
    integ_x(i)=(tot_chi2_bin(i)+tot_chi2_bin(i+1))/2.0
    integ_y(i)=sum(tot_prob_bin(i:num))/sum(tot_prob_bin(1:num))
  end do

  ! Lets measure the variance of the distribution, by measuring the
  ! 68.2 percent confidence interval.
  sigma = (integ_x(minloc(abs(0.159-integ_y),1)) - &
           integ_x(minloc(abs(0.841-integ_y),1)))/2.0
  print*, 'The sigma of the distribution is ', sigma, 'which should'
  print*, 'be approximately the square root of twice the number of free '
  print*, 'parameters ', nstars-nparams, 'i.e. ', sqrt(real(2*(nstars-nparams)))

  ! Note this is the exectation value were there no free parameters.
  expect=sum(tot_chi2_bin(1:num)*tot_prob_bin(1:num))/sum(tot_prob_bin(1:num))

  ! Not sure what one is supposed to do the the x-axis to allow for
  ! free parameters for the ditribution of individual points.
  open(unit=12, file='one.tau')
  write(12,'(/,/)')
  do j=1, n_one
    write(12,*) tau_one(j), sum(frac_one(j:n_one))/real(nstars)
  end do
  close(12)

  ! Now apply the scaling to allow for free parameters.
  integ_x = (integ_x-expect)
  integ_x = integ_x*real(nstars-nparams)/real(nstars)
  integ_x=integ_x+expect-real(nparams)

  open(unit=12, file='integ.tau')
  write(12,'(/,/)')
  do i=1, num-1
    write(12,*) integ_x(i), integ_y(i)
  end do
  close(12)

  open(unit=12, file='tau.diff')
  write(12,'(/,/)')
  do i=1, num
    write(12,*) tot_chi2_bin(i), tot_prob_bin(i)
  end do
  close(12)

  print*, ' '
  print*, 'Reading best-fit tau^2 value from grid.fit'
  call nxtark_in('grid.fit')
  iflag=inpark(naxis, data, axdata)
  print*, 'Pr(', minval(data), ') is ', &
  quadint(integ_x, integ_y, minval(data), sflag), &
  ' for ', nparams, ' free parameters.'  
  print*, 'Reduced tau^2 is ', 1.0 + &
  ((minval(data)-expect)/real(nstars-nparams))

end program tau
