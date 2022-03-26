module quad

implicit none

contains

  real function quadintflg(x, y, xval, flag, xflg, yflg)

    ! Performs quadratic interpolation.
    ! flag returned as 'OK' for success, '1' if outside x-range,
    ! and hence extrapolated

    real, dimension(:), intent(in) :: x, y
    character(len=*), dimension(:), intent(in), optional :: xflg, yflg
    real, intent(in) :: xval
    character(len=*), intent(out) :: flag

    integer :: i, j, ifirst, ilast, extrap
    real :: a, b, c, work, work1

    ! Start as we hope to go on.
    flag='OK'

    if (size(x) /= size(y)) then
      print*, 'Error in arrays provided to quadintflag, they are of different sizes.'
      stop
    end if

    ! Find first and last good values.
    ifirst=1
    ilast=size(x,1)
    if (present(xflg)) then
      xflglo: do i=ifirst, ilast
        ifirst=i
        if (xflg(i) == 'OK') exit xflglo
      end do xflglo
      xflghi: do i=ilast, ifirst, -1
        ilast=i
        if (xflg(i) == 'OK') exit xflghi
      end do xflghi
    end if
    if (present(yflg)) then
      yflglo: do i=ifirst, ilast
        ifirst=i
        if (yflg(i) == 'OK') exit yflglo
      end do yflglo
      yflghi: do i=ilast, ifirst, -1
        ilast=i
        if (yflg(i) == 'OK') exit yflghi
      end do yflghi
    end if

    if (ilast-ifirst <= 0) then
      ! Cannot interpolate with just one datapoint!
      flag='1'
      quadintflg=y(ifirst)
    end if      

    if (ilast-ifirst == 2) then
      ! Only enough data for linear interpolation.
      quadintflg=y(ifirst) + &
      (y(ilast)-y(ifirst))*(x(ifirst)-xval)/(x(ifirst)-x(ilast))
    else
      !i=minloc(abs(x-xval),1) 
      i=locate_nearest(x(ifirst:ilast), xval, extrap)

      ! The line below has the effect of
      ! always interpolating between the two points either side and the
      ! point whose index is one more.  This had the advantage of compatibility with 
      ! the numerical recipes routine we used to use, and to ensures 
      ! we do not change the interpolating polynomial between points, only on a point
      ! (where the value is forced) so there are no gliches.
      ! The down side is that it restricts your range of available values.
      if ((xval-x(i))/(x(size(x))-x(1)) > 0.0) i=i+1

      ! Don't run off the top.
      if (i >= ilast) i=ilast-1
      ! Or the bottom.
      if (i <= ifirst) i=ifirst+1
      work1=((x(i-1)*x(i-1)-x(i)*x(i))/(x(i-1)-x(i)) - &
      (x(i)*x(i)-x(i+1)*x(i+1))/(x(i)-x(i+1)))
      if (abs(work1) > 10.0*tiny(work1)) then
        a = ((y(i-1)-y(i))/(x(i-1)-x(i)) - (y(i)-y(i+1))/(x(i)-x(i+1)))/work1
      else
        ! The quadratic term is small
        a=0.0
      end if
      b = (y(i-1)-y(i)-a*(x(i-1)*x(i-1)-x(i)*x(i)))/(x(i-1)-x(i))
      c = y(i-1) - a*x(i-1)*x(i-1) - b*x(i-1)
      quadintflg = a*xval*xval + b*xval + c
      if (present(xflg)) then
        do j=i-1, i+1
          if (xflg(j) /= 'OK') flag=xflg(j)
        end do
      end if
      if (present(yflg)) then
        do j=i-1, i+1
          if (yflg(j) /= 'OK') flag=yflg(j)
        end do
      end if
    end if

    ! Did we extrapolate?
    if (xval < min(x(ifirst),x(ilast))) then
      if (x(ifirst) < x(ilast)) then
        quadintflg=y(ifirst)
      else
        quadintflg=y(ilast)
      end if
      flag='1'
    else if (xval > max(x(ifirst),x(ilast))) then
      if (x(ifirst) > x(ilast)) then
        quadintflg=y(ifirst)
      else
        quadintflg=y(ilast)
      end if
      flag='1'
    end if

  end function quadintflg

  real function quadint(x, y, xval, flag)

    real, dimension(:), intent(in) :: x, y
    real, intent(in) :: xval
    character(len=*), intent(out) :: flag

    character(len=10), allocatable, dimension(:) :: xflg, yflg

    allocate(xflg(size(x)), yflg(size(y)))

    xflg='OK'
    yflg='OK'

    quadint=quadintflg(x, y, xval, flag, xflg, yflg)

    deallocate(xflg, yflg)

  end function quadint

  real function linint(x, y, xval, flag)

    real, dimension(:), intent(in) :: x, y
    real, intent(in) :: xval
    character(len=*), intent(out) :: flag

    character(len=10), allocatable, dimension(:) :: xflg, yflg

    allocate(xflg(size(x)), yflg(size(y)))

    xflg='OK'
    yflg='OK'

    linint=linintflg(x, y, xval, flag, xflg, yflg)

    deallocate(xflg, yflg)

  end function linint

  real function linintflg(x, y, xval, flag, xflg, yflg)

    real, dimension(:), intent(in) :: x, y
    ! Note these are interchangeable.  So if you only have one set of
    ! flags you can call it in just one set.
    character(len=*), dimension(:), intent(in), optional :: xflg, yflg
    real, intent(in) :: xval
    character(len=*), intent(out) :: flag
    
    integer :: i, extrap, j

    if (size(x) /= size(y)) then
      print*, 'Error in arrays provided to quadintflag, they are of different sizes.'
      stop
    end if

    flag='OK'

    !i=minloc(abs(x-xval),1) 
    i=locate_nearest(x, xval, extrap)
    ! Make i the data point above xval.
    if ((xval-x(i))/(x(size(x))-x(1)) > 0.0) i=i+1
    ! Don't run off the top.
    if (i >= size(x)) i=size(x)
    ! Or the bottom.
    if (i <= 1) i=2
    linintflg=y(i-1) + (y(i)-y(i-1))*(x(i-1)-xval)/(x(i-1)-x(i))

    ! Was either datapoint flagged?
    if (present(xflg)) then
      if (xflg(i-1) /= 'OK') flag=xflg(i-1)
      if (xflg(i) /= 'OK') flag=xflg(i)
    end if
    if (present(yflg)) then
      if (yflg(i-1) /= 'OK') flag=yflg(i-1)
      if (yflg(i) /= 'OK') flag=yflg(i)
    end if

    if (extrap /=0) flag='1'

  end function linintflg

  integer function locate_nearest(xdata, xval, extrap)

    ! Given an array which is monotonically increasing or decreasing, this function
    ! finds the closest point to xval.

    real, intent(in), dimension(:) :: xdata
    real, intent(in) :: xval
    integer, optional, intent(out) :: extrap

    integer :: npts, ilow, ihigh, itest
    real :: change

    npts=size(xdata,1)

    change=1.0
    if (xdata(npts) < xdata(1)) change=-1.0

    if (present(extrap)) extrap=0
    if (change*(xval-xdata(npts)) > 0.0) then
      locate_nearest=npts
      if (present(extrap)) extrap=1
      return
    else if (change*(xval-xdata(1)) < 0.0) then
      locate_nearest=1
      if (present(extrap)) extrap=-1
      return
    end if

    ilow=1
    ihigh=npts
    do 
      itest=(ilow+ihigh)/2
      !if (change*(xval-xdata(itest)) >= 0.0) then
      !  ilow=itest
      !else if (change*(xval-xdata(itest)) <= 0.0) then
      !  ihigh=itest
      !end if
      if (change*(xval-xdata(itest)) >= 0.0) then
        ilow=itest
      else
        ihigh=itest
      end if
      if (ihigh - ilow < 3) then
        locate_nearest=minloc(abs(xdata(ilow:ihigh)-xval),1)+ilow-1
        return
      end if
    end do

  end function locate_nearest

end module quad

