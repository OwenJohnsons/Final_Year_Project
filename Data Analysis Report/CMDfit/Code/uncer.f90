module uncer_subs

  use likelihood_mod
  use quad

  implicit none

contains

  real function sum_tau2(tau2_grid, axdata, axisnam, min1, max1, min2, max2)

    ! This sums the probability from a given set of tau2 pixels.
    ! The only clever bit is that if it finds one of the axes is a
    ! log axis it will weight the pixels so that the sum is done
    ! as though it were a linear axis.

    real, dimension(:,:), intent(in) :: tau2_grid, axdata
    character(len=*), dimension(:), intent(in) :: axisnam
    ! The region to be summed.  Clumsy to have to put it in, but otherwise
    ! you have to do clumsy things with subsetting axdata.
    integer, intent(in) :: min1, max1, min2, max2

    ! Locals.
    integer :: iaxis1, iaxis2
    real, parameter :: big=2.0*log(huge(tau2_grid))
    ! Set these up to be the size of the pixels.
    real, dimension(size(tau2_grid,1)) :: delta1
    real, dimension(size(tau2_grid,2)) :: delta2

    delta1=1.0
    if (axisnam(1)(1:3) == 'LOG') delta1(1:size(delta1))=10.0**axdata(1:size(tau2_grid,1),1)
    delta2=1.0
    if (axisnam(2)(1:3) == 'LOG') delta2(1:size(delta2))=10.0**axdata(1:size(tau2_grid,2),2)

    sum_tau2=0.0
    do iaxis2=min2, max2
      do iaxis1=min1, max1
        if (tau2_grid(iaxis1, iaxis2) < big/2.0) &
        sum_tau2=sum_tau2 + &
        (delta1(iaxis1)*delta2(iaxis2)*exp(-0.5*(tau2_grid(iaxis1,iaxis2))))
      end do
    end do

  end function sum_tau2

  subroutine create_cum_prob(tau2_grid, axgrid, axnam, tau2_cum, &
  cum_prob, one_sig_tau2)

    ! Creates the plot of tau^2 against cumulative probability.

    ! The input tau2 grid.
    real, dimension(:,:), intent(in) :: tau2_grid, axgrid
    character(len=*), dimension(:) :: axnam
    ! The output plot.
    real, allocatable, dimension(:), intent(out) :: tau2_cum, cum_prob
    ! And the tau2 corresponding to one sigma.
    real :: one_sig_tau2

    ! Local variables.
    integer, dimension(2) :: naxis
    integer :: icount, i1, i2, i
    real :: tau2_min, tot_prob
    character(len=10) :: flag
    real, allocatable, dimension(:) :: pix_prob

    naxis(1)=size(tau2_grid,1)
    naxis(2)=size(tau2_grid,2)
    allocate(tau2_cum(naxis(1)*naxis(2)), &
    pix_prob(naxis(1)*naxis(2)), cum_prob(naxis(1)*naxis(2)))

    ! The numerics can fail if tau**2 is too high.
    tau2_min=minval(tau2_grid)

    icount=0
    do i2=1, naxis(2)
      do i1=1, naxis(1)
        icount=icount+1
        tau2_cum(icount)=tau2_grid(i1,i2)
        ! Now mulitply by the area of the pixel.
        pix_prob(icount)=sum_tau2(tau2_grid-tau2_min, &
        axgrid, axnam, i1, i1, i2, i2)
        !print*, pix_prob(icount)
      end do
    end do
    call tau_sort(tau2_cum, pix_prob)

    tot_prob=sum_tau2(tau2_grid-tau2_min, axgrid, axnam, 1, &
    naxis(1), 1, naxis(2))

    do icount=1, naxis(1)*naxis(2)
      cum_prob(icount)=&
      100.0*sum(pix_prob(1:icount))/tot_prob
    end do

    ! Find the 1 sigma contour by quadratic interpolation.
    one_sig_tau2=quadint(cum_prob, tau2_cum, 68.26, flag)

  end subroutine create_cum_prob



  subroutine conf2(tau2_grid, axgrid, axnam, long_axis)

    ! Print out the positions above and below which (1-0.68)/2 
    ! of the probability lies.

    ! The tau^2 grid and associated axis information.
    real, intent(in), dimension(:,:) :: tau2_grid, axgrid
    character(len=*), dimension(2) :: axnam
    ! Which axis is NOT to be summed over.
    integer, intent(in) :: long_axis

    real, parameter :: conf_level=(1.0-0.6826)/2.0
    !real, parameter :: conf_level=0.5
    integer :: i, ilow, ihigh
    character(len=8) :: flag
    real :: low, high, best, tot_prob, sum_prob
    integer, dimension(2) :: iwhere, naxis

    real, allocatable, dimension(:) :: cum

    ! You can't use size in an argument with the 64-bit machines, or it takes
    ! it to be i8.  So this is Joseph Stead's fix.
    naxis(1)=size(tau2_grid,1)
    naxis(2)=size(tau2_grid,2)

    sum_prob=sum_tau2(tau2_grid, axgrid, axnam, 1, naxis(1), 1, &
    naxis(2))
    tot_prob=0.0
    allocate(cum(size(tau2_grid,long_axis)))
    from_bottom: do i=1, size(tau2_grid,long_axis)
      if (long_axis == 1) then
        tot_prob=tot_prob+sum_tau2(tau2_grid, axgrid, axnam, &
        i, i, 1, naxis(2))
      else
        tot_prob=tot_prob+sum_tau2(tau2_grid, axgrid, axnam, &
        1, naxis(1), i, i)
      end if
      cum(i)=tot_prob/sum_prob
    end do from_bottom
    ! Find the confidence limits by quadratic interpolation.
    low=quadint(cum, axgrid(1:size(cum),long_axis), conf_level, flag)
    ilow=minloc(abs(low-axgrid(1:size(cum),long_axis)),1)
    high=quadint(cum, axgrid(1:size(cum),long_axis), 1.0-conf_level, flag)
    ihigh=minloc(abs(high-axgrid(1:size(cum),long_axis)),1)

    iwhere=minloc(tau2_grid)
    best=axgrid(iwhere(long_axis), long_axis)
    print*, ' '
    if (trim(axnam(long_axis)(1:5)) == 'LOG10') then
      print*, 'For ', &
      trim(axnam(long_axis)(6:len(axnam(long_axis))))//'/10^6', &
      ' one parameter of interest.'
      best=10.0**(best-6.0)
      high=10.0**(high-6.0)
      low= 10.0**(low -6.0)
    else
      print*, 'For ', trim(axnam(long_axis)), ' one parameter of interest.'
    end if
    if (best < low) then
      ! Ah, this can only be an upper limit.
      high=quadint(cum, axgrid(1:size(cum),long_axis), 1.0-2.0*conf_level, flag)
      print*, 'Since the best value is ', best, ' and the lower limit is '
      print*, low, ' we can only derive an upper limit.'
      print*, 100.0*(1.0-2.0*conf_level), ' percent limit is <', high 
    else
      print*, best, '+', high-best, low-best
      print*, 'i.e. ', low, ' -> ', high, '(=', ihigh-ilow, ' grid pixels.)'
      best=(low+high)/2.0
      print*, best, '+/-', (high-low)/2.0 ,'(symmetrised)'
    end if
  
  end subroutine conf2

end module uncer_subs

program uncer

  ! Calculates the uncertainties using a new experimental method.

  use ark_file_io
  use uncer_subs
  use quad

  implicit none

  integer :: iflag
  real, allocatable, dimension(:,:) :: axgrid, tau2_grid
  real, allocatable, dimension(:) :: cum_prob, tau2_cum
  real, allocatable, dimension(:,:) :: oneDgrid
  integer, dimension(2) :: naxis

  integer :: i, i1, i2, iaxis
  real :: tot_prob
  ! The minimum tau**2, removed to make numerics work.
  real :: tau2_min
  real :: work

  real :: one_sig_tau2, one_sig_1d_tau2
  character(len=80), dimension(2) :: axnam

  call nxtark_in('grid.fit')
  iflag=inpark(naxis, tau2_grid, axgrid)
  do iaxis=1, 2
    if (naxis(iaxis) == 1) cycle
    if (iaxis == 1) then 
      iflag=get_header_s('CTYPE1', axnam(1))
      if (iflag < 0) axnam(1)='AXIS 1'
    else if (iaxis == 2) then
      iflag=get_header_s('CTYPE2', axnam(2))
      if (iflag < 0) axnam(2)='AXIS 2'
    end if
  end do

  call create_cum_prob(tau2_grid, axgrid, axnam, tau2_cum, cum_prob, &
  one_sig_tau2)

  print*, 'The 68 percent confidence contour is at tau^2 =', one_sig_tau2

  call nxtark_out('uncer.out')
  call typark(2)
  iflag=makark(naxis(1)*naxis(2), cum_prob, tau2_cum)

  !deallocate(cum_prob, tau2_cum)

  tau2_min=minval(tau2_grid)
  tau2_grid=tau2_grid-tau2_min
  do iaxis=1, 2
    if (naxis(iaxis) == 1) cycle
    call conf2(tau2_grid, axgrid, axnam, iaxis)
  end do

end program uncer
