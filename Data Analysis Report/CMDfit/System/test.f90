program test

  use ark_file_io

  implicit none

  real, dimension(3,3,3) :: data
  real, dimension(3,3) :: axdata
  integer, dimension(3) :: naxis
  integer :: iflag

  naxis=3

  data=5.0
  axdata=0.0

  call nxtark_out('test.fits')
  iflag=makark(naxis, data, axdata)

end program test
