program for_gaia

  ! Given a FITS file and a set of data points (in cluster format), this
  ! programme creates a file which can be plotted on the image using the 
  ! regions menu in GAIA.

  use ark_file_io
  use define_star

  implicit none

  integer, dimension(2) :: naxis
  real, dimension(:,:), allocatable :: data, axdata

  integer :: nstars, ncol
  character(len=3), dimension(mcol) :: colstr
  type(a_star), dimension(:), allocatable :: star

  integer :: iflag, istar
  integer :: i1, i2, ilow, ihigh
  real :: start1, start2, step1, step2, rlow, rhigh, r1, r2

  integer :: ipos, icol, icol1, icol2
  character(len=50) :: ofname

  iflag=read_cluster_file(nstars, ncol, colstr, star)
  call lstcls_in(ofname)
  ipos=index(ofname, '.')
  if (ipos < 1) ipos=len_trim(ofname)+1
  ofname=ofname(1:ipos-1)//'.ard'

  if (ncol > 2) then
    print*, 'Which number colour do you want as y-axis?'
    do icol=1, ncol
      print*, icol, colstr(icol)
    end do
    read(*,*) icol1
    print*, 'Which number colour do you want as x-axis?'
    do icol=1, ncol
      print*, icol, colstr(icol)
    end do
    read(*,*) icol2
  else
    icol1=1
    icol2=2
  end if
    
  print*, '> Give the FITS file for the best fit to the data.'
  iflag=inpark(naxis, data, axdata)

  open(unit=35, file=ofname)

  step1=(axdata(naxis(1),1)-axdata(1,1))/real(naxis(1)-1)
  start1=axdata(1,1)-step1
  step2=(axdata(naxis(2),2)-axdata(1,2))/real(naxis(2)-1)
  start2=axdata(1,2)-step2

  do istar=1, nstars

    if (star(istar)%col(icol1)%flg/='OO') cycle
    if (star(istar)%col(icol2)%flg/='OO') cycle

    ! And write a file of positions in X-Y space for Gaia.
    !i1=minloc(abs(star(istar)%col(icol2)%data-axdata(1:naxis(1),1)),1)
    !i2=minloc(abs(star(istar)%col(1    )%data-axdata(1:naxis(2),2)),1)
    r1=(star(istar)%col(icol2)%data-start1)/step1
    r2=(star(istar)%col(icol1)%data-start2)/step2
    ! Circle at position.
    ! For VBV diagrams set the circle radius to 2.0
    ! for the EBV diagrams use 4.0
    !write(34,*) 'CIRCLE(', i1, ',', i2, ',', 5.0, ')'
    write(35,*) 'CIRCLE(', r1, ',', r2, ',', 5.0, ')'

    ! And error bars.
    !ilow =minloc(abs(star(istar)%col(icol2)%data-star(istar)%col(icol2)%err-axdata(1:naxis(1),1)),1)
    !ihigh=minloc(abs(star(istar)%col(icol2)%data+star(istar)%col(icol2)%err-axdata(1:naxis(1),1)),1)
    rlow =(star(istar)%col(icol2)%data-star(istar)%col(icol2)%err-start1)/step1
    rhigh=(star(istar)%col(icol2)%data+star(istar)%col(icol2)%err-start1)/step1
    !if (ilow /= ihigh) write(34,*) 'LINE(', ilow, ',', i2, ',', ihigh, ',', i2, ')'
    if (abs(rlow-rhigh) > tiny(rlow)) write(35,*) 'LINE(', rlow, ',', r2, ',', rhigh, ',', r2, ')'

    !ilow =minloc(abs(star(istar)%col(    1)%data-star(istar)%col(    1)%err-axdata(1:naxis(2),2)),1)
    !ihigh=minloc(abs(star(istar)%col(    1)%data+star(istar)%col(    1)%err-axdata(1:naxis(2),2)),1)
    rlow =(star(istar)%col(icol1)%data-star(istar)%col(icol1)%err-start2)/step2
    rhigh=(star(istar)%col(icol1)%data+star(istar)%col(icol1)%err-start2)/step2
    if (abs(rlow-rhigh) > tiny(rlow)) write(35,*) 'LINE(', r1, ',', rlow, ',', r1, ',', rhigh, ')'    
    !if (ilow /= ihigh) write(34,*) 'LINE(', i1, ',', ilow, ',', i1, ',', ihigh, ')'    

  end do

  close(35)

end program for_gaia
