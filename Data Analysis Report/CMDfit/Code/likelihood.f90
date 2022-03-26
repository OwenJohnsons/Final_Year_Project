module likelihood_mod

  implicit none

  contains

    subroutine tau_sort(tau_data, pix_prob)

      real, intent(inout), dimension(:)  :: tau_data
      real, optional, intent(inout), dimension(:)  :: pix_prob

      integer :: k, l, m
      real :: safe, safe_pix

      outer: do k=2, size(tau_data)
        safe=tau_data(k)
        if (present(pix_prob)) safe_pix=pix_prob(k)
        do l=1,k-1
          if (safe < tau_data(l)) then
            do m=k, l+1, -1
              tau_data(m)=tau_data(m-1)
              if (present(pix_prob)) pix_prob(m)=pix_prob(m-1)
            end do
            tau_data(l)=safe
            if (present(pix_prob)) pix_prob(l)=safe_pix
            cycle outer
          endif
        end do
      end do outer

    end subroutine tau_sort


  subroutine natural_norm(mag_max, mag_min, data, axdata)

    ! This subroutine normalises the model, so that it integrates to
    ! one between the faintest and brightest stars.

    ! Bug corrected in 2012 where it dod not multiply by the pixel area.
    
    real, intent(in) :: mag_max, mag_min
    real, intent(inout), dimension(:,:) :: data, axdata

    real :: pix_area
    integer :: imag_low, imag_high
    integer, dimension(2) :: naxis

    ! Find the area of a pixel.
    naxis(1)=size(data,1)
    naxis(2)=size(data,2)
    pix_area = (abs(axdata(naxis(1),1)-axdata(1,1))/real(naxis(1)-1)) &
             * (abs(axdata(naxis(2),2)-axdata(1,2))/real(naxis(2)-1))
        
    ! Find the nearest pixel position in mag (avoiding problems if
    ! magnitude axes run in the wrong direction) for the brightest
    ! and faintest stars.
    imag_low =min( minloc(abs( mag_min-axdata(1:naxis(2),2) ),1), &
                   minloc(abs( mag_max-axdata(1:naxis(2),2) ),1))
    imag_high=max( minloc(abs( mag_min-axdata(1:naxis(2),2) ),1), &
                   minloc(abs( mag_max-axdata(1:naxis(2),2) ),1))
    data=data/(pix_area*sum(data(:,imag_low:imag_high)))

  end subroutine natural_norm

  real function likelihood(data, axdata, grad, icol, imag, sig_col, &
    sig_mag, correlated, flag, sum_gauss)
    
    ! The unsmoothed image and axes.
    real, dimension(:,:), intent(in) :: data, axdata
    real, dimension(:), allocatable :: grad
    ! The co-ordinates of the position being considered.
    integer, intent(in) :: icol, imag
    real, intent(in) :: sig_col, sig_mag
    ! If correlated is zero the colour and magnitude uncertainties are taken
    ! to be uncorrelated.  If its 1, they are correlated as would
    ! be expected in I vs R-I; if its -1 as for V vs V-I (i.e. its whether 
    ! you add or subtract the colour to make the other magnitude).  
    ! If the uncertainties are correlated, make sure sig_col > sig_mag.
    integer, intent(in) :: correlated
    character, intent(out) :: flag

    real, optional, intent(out) :: sum_gauss

    ! We cut off the Gaussian smoothing the 2D distribution at some sigma.
    real, parameter :: sig_cut=4.0
    real, parameter :: pi=3.141592654

    ! Gradient in mag-mag space, and uncertainty in second magnitude.
    real :: mag_grad, mag_sig
    integer :: iicol, iimag, i1, i2
    real :: col, mag
    real :: work
    real, save :: area, pix1, pix2
    real :: d_col, d_mag, gauss, rho, gauss_norm
    logical, save :: first=.true.
    
    if (first) then
      ! Find the area of a pixel.
      pix1=(axdata(size(data,1),1)-axdata(1,1))/real(size(data,1)-1)
      pix2=(axdata(size(data,2),2)-axdata(1,2))/real(size(data,2)-1)
      area=abs(pix1*pix2)
      first=.false.
    end if

    ! Find the nearest positions in X and Y.
    col=axdata(icol,1)
    mag=axdata(imag,2)

    ! How far out do we need to go? 
    iicol=nint(sig_cut*sig_col/abs(pix1))
    iimag=nint(sig_cut*sig_mag/abs(pix2))

    ! For correlated uncertainties calculate mag_sig.
    if (correlated == 0) then
      mag_sig=sig_col
    else
      if (sig_mag > sig_col) then
        mag_sig=sig_mag/2.0
      else
        mag_sig=max(sqrt(sig_col**2.0 - sig_mag**2.0), sig_mag/2.0)
      end if
    end if
 
    ! Now we normalise the model so it is one within the magnitude
    ! range.


    if (present(sum_gauss)) then
      sum_gauss=0.0
      
      fudge: do i1=icol-iicol, icol+iicol

        ! This sums up the volume under the uncertainty ellipse which
        ! lies on the area of the model.  If it falls outside the colour
        ! range of the image we still add up the volume.  This is because
        ! monte trims in colour to save space.  In magnitude its a
        ! different story.  The model may well be non-zero outside the
        ! magnitude range, so we need to calculate what fraction of this 
        ! error ellipse falls within the monte grid.

        ! This works for all values of i1, even if they lie outside the grid.
        work=axdata(1,1) + (pix1*real(i1-1)) - col

        do i2=max(imag-iimag, 1), min(imag+iimag, size(data,2))
          
          d_mag=(axdata(i2,2)-mag)/sig_mag
          if (correlated == 0) then
            d_col=work/sig_col
          else if (correlated == 1) then
            d_col=(axdata(i2,2)-mag)+work
            d_col=d_col/mag_sig
          else if (correlated == -1) then
            d_col=(axdata(i2,2)-mag)-work
            d_col=d_col/mag_sig
          else 
            print*, 'Unacceptable value of ', correlated
            stop
          end if
          gauss=exp(-0.5*(d_mag**2.0 + d_col**2.0))
          sum_gauss=sum_gauss+(area*gauss)
          
        end do
      end do fudge

      sum_gauss=sum_gauss/(2.0*pi*mag_sig*sig_mag)

    end if


    ! Removed this if.  It used to check you didn't skip off the edge of the
    ! image.  Now the assmumption is that if you do, the values of the image
    ! pixels should be zero.
    !if (icol-iicol>0       .and. imag-iimag>0 .and. &
    !  icol+iicol<=size(data,1) .and. imag+iimag<=size(data,2)) then

      work=0.0

    
      outer: do i1=max(icol-iicol, 1), min(icol+iicol, size(data,1))

        if (allocated(grad)) then
          if (grad(i1) < tiny(grad(i1))) cycle outer
          ! I.e. if there are any non-zero data points at this colour.
          ! Calculate the normalisation.
          if (correlated == 0) then
            rho=(1.0/sig_mag)**2.0 + 1.0/(sig_col*grad(i1))**2.0
          else 
            mag_grad=1.0/(1.0-(1.0/grad(i1)))
            ! Avoid a numberical error if the two uncertainties are equal.
            mag_sig=max(sqrt(sig_col**2.0 - sig_mag**2.0), sig_mag/2.0)
            ! write(22,*) 'Got ', grad(i1), mag_grad, i1, mag_sig
            rho=(1.0/sig_mag)**2.0 + 1.0/(mag_sig*mag_grad)**2.0
            ! write(22,*) 'Which gave ', rho
          end if
          ! rho=sqrt(rho/pi)
          ! Correct the factor two error.
          rho=sqrt(0.5*rho/pi)
          gauss_norm=1.0
        else
          rho=1.0
          gauss_norm=1.0/(2.0*pi*mag_sig*sig_mag)
        end if



        do i2=max(imag-iimag, 1), min(imag+iimag, size(data,2))
          if (data(i1,i2) >= tiny(data(i1,i2))) then
            ! Find the colour difference from this point to the
            ! data point, and put it in terms of sigma.
            if (correlated == 0) then
              d_col=(axdata(i1,1)-col)/sig_col
            else if (correlated == 1) then
              ! But now comes the swindle to account for correlated
              ! uncertainties.  The colour distance we calculate, is 
              ! actually the distance in the other magnitude.
              d_col=(axdata(i2,2)-mag)+(axdata(i1,1)-col)
              ! Thus we have to divide by the uncertainty in the other 
              ! magnitude, reconstructed from the given uncertainties.
              d_col=d_col/mag_sig
            else if (correlated == -1) then
              d_col=(axdata(i2,2)-mag)-(axdata(i1,1)-col)
              d_col=d_col/mag_sig
            else 
              print*, 'Unacceptable value of ', correlated
              stop
            end if
            ! The same for the magnitude.  
            d_mag=(axdata(i2,2)-mag)/sig_mag
            ! gauss=exp(-1.0*(d_mag**2.0 + d_col**2.0))
            ! Correct the factor two error.
            gauss=exp(-0.5*(d_mag**2.0 + d_col**2.0))
            ! Now multiply the Gaussian at this point by the value of the
            ! model, the area of the pixel, and apply the normalisation.
            ! write(22,*) work, rho, area, gauss, data(i1,i2), i1, i2
            work=work+rho*area*gauss_norm*gauss*data(i1,i2)
          end if
        end do
      end do outer

      likelihood=work
      flag = 'O'

    !else

    !  likelihood=0.0
    !  flag = 'A'

    !end if

  end function likelihood

  real function new_likelihood(data, axdata, naxis, col, mag, sig_col, &
  sig_mag, correlated, sum_gauss)

    ! Consider a model which is in the form of an image (data). This routine
    ! finds the probability of there being a data point at co-ordinates 
    ! (col, mag), given the uncertainties (sig_col, sig_mag).

    ! The unsmoothed image and axes.
    real, dimension(:,:), intent(in) :: data, axdata
    integer, dimension(2) :: naxis
    ! The co-ordinates of the position being considered.
    real, intent(in) :: col, mag
    real, intent(in) :: sig_col, sig_mag
    ! If correlated is zero the colour and magnitude uncertainties are taken
    ! to be uncorrelated.  If its 1, they are correlated as would
    ! be expected in I vs R-I; if its -1 as for V vs V-I (i.e. its whether 
    ! you add or subtract the colour to make the other magnitude).  
    ! If the uncertainties are correlated, make sure sig_col > sig_mag.
    integer, intent(in) :: correlated

    real, optional, intent(out) :: sum_gauss

    ! We cut off the Gaussian smoothing the 2D distribution at some sigma.
    real, parameter :: sig_cut=4.0
    real, parameter :: pi=3.141592654

    ! Gradient in mag-mag space, and uncertainty in second magnitude.
    real :: mag_grad, mag_sig
    integer :: iicol, iimag, i1, i2
    integer :: icol, imag
    real :: work
    real, save :: area, pix1, pix2
    real :: d_col, d_mag, gauss, rho, gauss_norm

    ! Find the area of a pixel.
    pix1=(axdata(naxis(1),1)-axdata(1,1))/real(naxis(1)-1)
    pix2=(axdata(naxis(2),2)-axdata(1,2))/real(naxis(2)-1)
    area=abs(pix1*pix2)

    ! Find out which pixel numbers this new star corresponds to.
    ! Note that field stars could lie outside the model area.
    work=col-axdata(1,1)
    icol = 1 + nint(real(naxis(1)-1)*work/(axdata(naxis(1),1)-axdata(1,1)))
    work=mag-axdata(1,2)
    imag = 1 + nint(real(naxis(2)-1)*work/(axdata(naxis(2),2)-axdata(1,2)))
    !icol=minloc(abs(simcol-axdata(1:naxis(1),1)),1)
    !imag=minloc(abs(simmag-axdata(1:naxis(2),2)),1)

    ! How far out do we need to go? 
    iicol=nint(sig_cut*sig_col/abs(pix1))
    iimag=nint(sig_cut*sig_mag/abs(pix2))

    ! For correlated uncertainties calculate mag_sig.
    if (correlated == 0) then
      mag_sig=sig_col
    else
      if (sig_mag > sig_col) then
        mag_sig=sig_mag/2.0
      else
        mag_sig=max(sqrt(sig_col**2.0 - sig_mag**2.0), sig_mag/2.0)
      end if
    end if

    ! Now we normalise the model so it is one within the magnitude
    ! range.


    if (present(sum_gauss)) then
      sum_gauss=0.0

      fudge: do i1=icol-iicol, icol+iicol

        ! This sums up the volume under the uncertainty ellipse which
        ! lies on the area of the model.  If it falls outside the colour
        ! range of the image we still add up the volume.  This is because
        ! monte trims in colour to save space.  In magnitude its a
        ! different story.  The model may well be non-zero outside the
        ! magnitude range, so we need to calculate what fraction of this 
        ! error ellipse falls within the monte grid.

        ! This works for all values of i1, even if they lie outside the grid.
        work=axdata(1,1) + (pix1*real(i1-1)) - col

        do i2=max(imag-iimag, 1), min(imag+iimag, size(data,2))

          d_mag=(axdata(i2,2)-mag)/sig_mag
          if (correlated == 0) then
            d_col=work/sig_col
          else if (correlated == 1) then
            d_col=(axdata(i2,2)-mag)+work
            d_col=d_col/mag_sig
          else if (correlated == -1) then
            d_col=(axdata(i2,2)-mag)-work
            d_col=d_col/mag_sig
          else 
            print*, 'Unacceptable value of ', correlated
            stop
          end if
          gauss=exp(-0.5*(d_mag**2.0 + d_col**2.0))
          sum_gauss=sum_gauss+(area*gauss)

        end do
      end do fudge

      sum_gauss=sum_gauss/(2.0*pi*mag_sig*sig_mag)

    end if


    work=0.0

    outer: do i1=max(icol-iicol, 1), min(icol+iicol, size(data,1))

      rho=1.0
      gauss_norm=1.0/(2.0*pi*mag_sig*sig_mag)

      do i2=max(imag-iimag, 1), min(imag+iimag, size(data,2))
        if (data(i1,i2) >= tiny(data(i1,i2))) then
          ! Find the colour difference from this point to the
          ! data point, and put it in terms of sigma.
          if (correlated == 0) then
            d_col=(axdata(i1,1)-col)/sig_col
          else if (correlated == 1) then
            ! But now comes the swindle to account for correlated
            ! uncertainties.  The colour distance we calculate, is 
            ! actually the distance in the other magnitude.
            d_col=(axdata(i2,2)-mag)+(axdata(i1,1)-col)
            ! Thus we have to divide by the uncertainty in the other 
            ! magnitude, reconstructed from the given uncertainties.
            d_col=d_col/mag_sig
          else if (correlated == -1) then
            d_col=(axdata(i2,2)-mag)-(axdata(i1,1)-col)
            d_col=d_col/mag_sig
          else 
            print*, 'Unacceptable value of ', correlated
            stop
          end if
          ! The same for the magnitude.  
          d_mag=(axdata(i2,2)-mag)/sig_mag
          ! gauss=exp(-1.0*(d_mag**2.0 + d_col**2.0))
          ! Correct the factor two error.
          gauss=exp(-0.5*(d_mag**2.0 + d_col**2.0))
          ! Now multiply the Gaussian at this point by the value of the
          ! model, the area of the pixel, and apply the normalisation.
          ! write(22,*) work, rho, area, gauss, data(i1,i2), i1, i2
          work=work+rho*area*gauss_norm*gauss*data(i1,i2)
        end if
      end do
    end do outer

    new_likelihood=work

  end function new_likelihood

end module likelihood_mod
