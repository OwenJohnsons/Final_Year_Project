program tau2_monte

  ! This is the code which calculates the probability that the data are a good
  ! fit to the model by exmaining the value of tau^2.  

  ! This is the Monte Carlo version described in Bell, Mamajek & Naylor 
  ! (in prep).

  use define_star
  use ark_file_io
  use likelihood_mod
  use random

  implicit none

  integer :: iflag

  real, dimension(:,:), allocatable :: data, axdata
  integer, dimension(2) :: naxis
  integer :: iaxis1, iaxis2

  integer :: nstars, ncol
  character(len=5), dimension(mcol) :: colstr
  type(a_star), dimension(:), allocatable :: star

  ! The colours we are using.
  integer :: jcol1, jcol2

  ! The magnitude boundaries.  The ith star lies between
  ! boundaries i and i+1.
  real, dimension(:), allocatable :: boundaries
  integer :: n_magbins

  integer :: istar, jstar, jpts, icluster
  integer, parameter :: ncluster=1000

  ! The summed tau2 for a cluster.
  real, dimension(ncluster) :: tau2_cluster
  ! The tau2 for an individual star.
  real, allocatable, dimension(:) :: tau2_star
  integer :: ntau2_star
  logical :: lost_cluster

  ! And individual pixel.
  type a_pixel
     ! It has a colour, magnitude and a cumulative fraction.
     real :: col, mag, frac
  end type a_pixel

  ! A bin in magnitude.
  type a_magbin
     ! It has a number of pixels.
     integer :: npts, ipts
     type(a_pixel), allocatable, dimension(:) :: pixel
  end type a_magbin

  ! And we have an array of bins in magnitude.
  type(a_magbin), dimension(:), allocatable :: magbins

  real :: harvest
  ! A simulate star.
  real :: simcol, simmag

  integer :: ishift

  character :: flag

  real :: prob
  real :: tau2_min

  integer :: iminus, iband, correlated
  real :: delta_mag1, delta_mag2
  real :: min_mag, max_mag

  integer :: ilost, nparams

  double precision :: tau2_mean

  ! The area of the CMD and range of colour and magnitude
  real :: cmd_area, col_min, col_max, mag_min, mag_max
  ! Temporary storage.
  real :: work
  ! Has the correct pixel been found?
  logical :: found_pixel

  print*, &
  'First read in the best fitting model FITS file (normally best_model.fit).'
  iflag=inpark(naxis, data, axdata)
  if (iflag < 0) then
    print*, 'Error reading 2D isochrone file.'
    stop
  end if

  print*, 'And now the file of data points (normally fitted.cat).'
  iflag=read_cluster_file(nstars, ncol, colstr, star)
  print*, 'Read ', nstars, ' data points.'

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

  ! Sort by magnitude.
  call sort_star(star, nstars, 'increasing_col(?)%data', jcol1)

  ! Check the datapoints lie on the model.
  ! Find the maximum and minium magnitude values of the model.
  max_mag=maxval(axdata(1:naxis(2),2))
  min_mag=minval(axdata(1:naxis(2),2))
  
  if (star(1)%col(1)%data < min_mag) then
    print*, 'The minimum magnitude in the model is ', min_mag
    print*, 'but in the data is ', star(1)%col(1)%data
    print*, 'i.e. the data fall off the top of the model, which means the'
    print*, 'fitting is unreliable.'
    stop
  else if (star(nstars)%col(1)%data > max_mag) then
    print*, 'The maximum magnitude in the model is ', max_mag
    print*, 'but in the data is ', star(nstars)%col(1)%data
    print*, 'i.e. the data fall off the bottom of the model, which means the'
    print*, 'fitting is unreliable.'
    stop
  end if

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
    ! E.g. I vs R-I since R=(R-I)+I. 
    ! So delta (R-I) = delta R - delta I = delta_mag2 - delta_mag1
    correlated = 1
  else
    print*, 'Will subtract ', trim(colstr(jcol2)), ' from ', &
    trim(colstr(jcol1)), &
    ' to create ', trim(colstr(jcol2)(iminus+1:len(colstr(jcol2))))
    ! E.g. V vs V-I since V-(V-I)=I.  
    ! So delta (V-I) = delta V - delta I = delta_mag2 - delta_mag1
    correlated = -1
  end if
  if (sum(star%col(jcol1)%err) >= 0.99*sum(star%col(jcol2)%err)) then
    print*, 'The uncertanties in magnitude are larger than in colour so'
    print*, 'I will take them as uncorrelated.'
    correlated = 0
  end if

  allocate(boundaries(nstars+1))
  ! Set the boundaries to lie between each observed star.
  do istar=2, nstars
    boundaries(istar) = &
    (star(istar-1)%col(jcol1)%data + star(istar)%col(jcol1)%data)/2.0
  end do
  ! Set the top and bottom bin boundaries using the positions of the
  ! second brightest and faintest stars.
  boundaries(1) = star(1)%col(jcol1)%data - & 
  (star(2)%col(jcol1)%data-star(1)%col(jcol1)%data)/2.0
  boundaries(nstars+1) =  star(nstars)%col(jcol1)%data + & 
  (star(nstars)%col(jcol1)%data-star(nstars-1)%col(jcol1)%data)/2.0

  allocate(magbins(nstars))
  magbins%npts=0

  ! Count how many non-zero pixels there are in each magnitude range.
  do iaxis1=1, naxis(1)
    do iaxis2=1, naxis(2)
      if (data(iaxis1,iaxis2) > 2.0*tiny(data(iaxis1,iaxis2))) then
        do istar=1, nstars
          if (axdata(iaxis2,2)>boundaries(istar) .and. &
          axdata(iaxis2,2)<boundaries(istar+1)) then
            magbins(istar)%npts=magbins(istar)%npts+1
          end if
        end do
      end if
    end do
  end do

  !do istar=1, nstars
  !  print*, '**', boundaries(istar)
  !  print*, star(istar)%col(jcol1)%data, magbins(istar)%npts
  !end do
  !print*, '**', boundaries(nstars+1)    

  ! Allocate arrays which have enough elements for each pixel.
  do istar=1, nstars
    allocate(magbins(istar)%pixel(magbins(istar)%npts))
  end do

  ! Now set the colour and magnitude for each pixel, and set a cumulative 
  ! probability for that pixel.
  magbins%ipts=0
  do iaxis1=1, naxis(1)
    do iaxis2=1, naxis(2)
      if (data(iaxis1,iaxis2) > 2.0*tiny(data(iaxis1,iaxis2))) then
        do istar=1, nstars
          if (axdata(iaxis2,2)>boundaries(istar) .and. &
          axdata(iaxis2,2)<boundaries(istar+1)) then
            magbins(istar)%ipts=magbins(istar)%ipts+1
            magbins(istar)%pixel(magbins(istar)%ipts)%col=axdata(iaxis1,1)
            magbins(istar)%pixel(magbins(istar)%ipts)%mag=axdata(iaxis2,2)
            if (magbins(istar)%ipts == 1) then
              magbins(istar)%pixel(magbins(istar)%ipts)%frac = &
              data(iaxis1,iaxis2)
            else
              magbins(istar)%pixel(magbins(istar)%ipts)%frac = &
              magbins(istar)%pixel(magbins(istar)%ipts-1)%frac + &
              data(iaxis1,iaxis2)
            end if
          end if
        end do
      end if
    end do
  end do

  ! Now normalise the fractions.
  do istar=1, nstars
    if (magbins(istar)%npts > 0) then
      magbins(istar)%pixel(:)%frac = &
      magbins(istar)%pixel(:)%frac / &
      magbins(istar)%pixel(magbins(istar)%npts)%frac
    end if
  end do

  ! Normalise the model.
  call natural_norm(maxval(star(1:nstars)%col(jcol1)%data,1), &
  minval(star(1:nstars)%col(jcol1)%data,1), data, axdata)

  ! Find the area of the CMD covered by the datapoints.
  col_min=minval(star%col(1)%data,1)
  col_max=maxval(star%col(1)%data,1)
  mag_min=minval(star%col(2)%data,1)
  mag_max=maxval(star%col(2)%data,1)
  cmd_area = (col_max-col_min)*(mag_max-mag_min)

  ! Simulate a few clusters.
  tau2_cluster=0.0
  ilost=0
  allocate(tau2_star(ncluster*nstars))
  ntau2_star=0
  call set_seed()
  each_cluster: do icluster=1, ncluster
    ! print*, 'Doing cluster ', icluster, ' of ', ncluster
    lost_cluster=.false.
    each_star: do istar=1, nstars

      if (magbins(istar)%npts /= 0) then
        jstar=istar
      else
        ! Sometimes a star can end up with no pixels in its bin.
        find_shift: do ishift=1, nstars
          ! So look one star on either side.
          if (magbins(max(istar-ishift,1))%npts /= 0) then
            ! Perhaps the star a tad brighter is OK.
            jstar=istar-ishift
            exit find_shift
          else if (magbins(min(istar+ishift,nstars))%npts /= 0) then
            ! Or perhaps the one below.
            jstar=istar+ishift
            exit find_shift
          end if
          if (ishift == nstars) then
            print*, 'Two stars at around a magnitude of ', &
            star(istar)%col(1)%data
            print*, 'have the same magnitude in a way I cannot cope with.'
            stop
          end if
        end do find_shift
        print*, 'For the ', istar, 'th star out of ', nstars
        print*, 'magnitude ', star(istar)%col(1)%data
        print*, 'Using star ', jstar, ' of magnitude ', star(jstar)%col(1)%data
        Print*, 'instead, because there is no probability in the magnitude ' 
        print*, 'range of the other star.'
        
      end if
        
      ! So is it a field star or not?
      call random_number(harvest)
      if (harvest > star(istar)%col(3)%data) then
        ! Its a field star.  Place it randomly in colour.
        call random_number(harvest)
        simcol=col_min + harvest*(col_max-col_min)
        ! But somewhere in the bin magnitude range.
        call random_number(harvest)
        simmag=boundaries(istar) + &
        harvest*(boundaries(istar+1)-boundaries(istar))
        !print*, '-',  boundaries(istar)
        !print*, star(istar)%col(1)%data, simmag
        !print*, '+', boundaries(istar+1)
        found_pixel=.true.
      else
        found_pixel=.false.
        call random_number(harvest)
        ! Now go through each pixel in the magnitude bin.
        find_pixel: do jpts=1, magbins(jstar)%npts
          ! harvest is always less than 1, but is it greater than the cumulative
          ! fraction?
          if (magbins(jstar)%pixel(jpts)%frac > harvest) then
            simmag = magbins(jstar)%pixel(jpts)%mag
            simcol = magbins(jstar)%pixel(jpts)%col
            found_pixel=.true.
            exit find_pixel
          end if
        end do find_pixel
      end if

      ! Now deal with the uncertainties.
      ! The magnitude is always easy to deal with.
      delta_mag1=rnd_gauss(0.0, star(istar)%col(jcol1)%err)
      ! The colour is harder.
      if (correlated == 0) then
        ! As is the colour if the uncertanties are uncorrelated.
        delta_mag2=rnd_gauss(0.0, star(istar)%col(jcol2)%err)
      else
        ! First find the uncertainty in the magnitude which only 
        ! appears in the colour, and sample from it.
        delta_mag2=rnd_gauss(0.0,sqrt(star(istar)%col(jcol2)%err**2.0 - &
        star(istar)%col(jcol1)%err**2.0))
        if (correlated == 1) then
          delta_mag2 = delta_mag2 - delta_mag1
        else if (correlated == -1) then
          delta_mag2 = delta_mag1 - delta_mag2
        end if
      end if

      simmag = simmag + delta_mag1
      simcol = simcol + delta_mag2

      if (.not. found_pixel) then
        print*, 'Error finding pixel'
        stop
      end if

      ! Now find the tau2 for this simulated data point.
      prob=new_likelihood(data, axdata, naxis, simcol, simmag, &
      star(istar)%col(jcol2)%err, star(istar)%col(jcol1)%err, &
      correlated)

      ! Allow for the possibility it's not a member.
      prob = prob*star(istar)%col(3)%data + &
      (1.0-star(istar)%col(3)%data)/cmd_area

      ntau2_star=ntau2_star+1
      if (prob < 2.0*tiny(prob)) then
        print*, 'Lost ', simcol, star(istar)%col(jcol2)%err, simmag, &
        star(istar)%col(jcol1)%err
        tau2_star(ntau2_star)=huge(tau2_star(ntau2_star))
        if (.not. lost_cluster) then 
          ilost=ilost+1
          lost_cluster=.true.
          tau2_cluster(icluster)=huge(tau2_cluster(icluster))
        end if
      else
        tau2_star(ntau2_star)=-2.0*log(prob)
      end if
      if (.not. lost_cluster) tau2_cluster(icluster)=tau2_cluster(icluster)+tau2_star(ntau2_star)
      !write(34,*) istar, icluster, tau2_star(ntau2_star)

    end do each_star

    !write(33,*) icluster, tau2_cluster(icluster)

  end do each_cluster

  ! print*, sum(tau2_cluster)/real(ncluster)

  print*, 'Number of very high tau^2 fits ', ilost
  print*, ' '

  ! Now apply the scaling to allow for free parameters.
  ! For this we first need to find the mean.
  tau2_mean=0.0d0
  do jpts=1, ncluster
    if (tau2_cluster(jpts) < 0.5*huge(tau2_cluster(jpts))) &
    tau2_mean = tau2_mean + dble(tau2_cluster(jpts))
  end do
  tau2_mean=tau2_mean/dble(ncluster-ilost)
  ! Then ask for the number of free parameters.
  print*, '> Give the number of free parameters.'
  read(*,*) nparams
  ! And finally apply the scaling.
  tau2_cluster = tau2_cluster-tau2_mean
  tau2_cluster = tau2_cluster*sqrt(real(nstars-nparams)/real(nstars))
  tau2_cluster=tau2_cluster+tau2_mean-real(nparams)

  ! Write out the distribution for the fitted points.
  open(unit=23, file='model_distrib.tau')
  write(23,*) '# This is the expected distribution of individual tau^2s'
  write(23,*) '# for each star.'
  write(23,*) '# tau2  fraction'
  call tau_sort(tau2_star(1:ntau2_star))
  do jpts=1, ntau2_star
    write(23,*) tau2_star(jpts), 1.0-(real(jpts)/real(ntau2_star))
  end do
  close(23)

  ! Write out the full distribution.
  open(unit=23, file='pr_tau2.dat')
  write(23,*) '# This is the expected distribution of total tau^2'
  write(23,*) '# for the cluster.'
  write(23,*) '# tau2  fraction'
  call tau_sort(tau2_cluster)
  do jpts=1, ncluster
    if (tau2_cluster(jpts) < 0.5*huge(tau2_cluster(jpts))) &
    write(23,*) tau2_cluster(jpts), 1.0-(real(jpts)/real(ncluster))
  end do
  close(23)

  ! And find where the measured tau2 lies.
  print*, 'Reading best-fit tau^2 value from grid.fit'
  call nxtark_in('grid.fit')
  iflag=inpark(naxis, data, axdata)
  if (iflag < 0) then
    print*, 'Failed to read file grid.fit, inpark returned flag ', iflag
    stop
  end if
  tau2_min=minval(data)
  icluster=count(tau2_cluster > tau2_min)

  print*, 'Pr(', minval(data), ') is ', real(icluster)/real(ncluster-ilost)

end program tau2_monte
