program ascii2cluster

  ! This is a code which will read a file which is a series of columns of IDs, colours
  ! and uncertainties and write them into a cluster format file.  You will be prompted
  ! for which columns are which colours.

  use define_star

  implicit none

  character(len=10), dimension(mcol) :: colstr
  type(a_star), dimension(:), allocatable :: star
  integer :: iostat, ncol, iunit, istar, nstars, icol, ncols
  integer :: column_id
  character(len=15), allocatable, dimension(:) :: cols
  integer, allocatable, dimension(:) :: column_col, column_uncer
  character(len=30) :: filnam

  print*, '> How many colours do you want?'
  read(*,*) ncol

  allocate(column_col(ncol), column_uncer(ncol))

  do icol=1, ncol
    print*, '> Give the name for colour ', icol
    read(*,*) colstr(icol)
    print*, '> Give the column numbers for ', trim(colstr(icol)), ' and its uncertainty.'
    read(*,*) column_col(icol), column_uncer(icol)
  end do

  print*, '> Which column should be used for the ID number (or zero if none).'
  read(*,*) column_id

  ncols=max(max(maxval(column_col), maxval(column_uncer)), column_id)
  allocate(cols(ncols))

  print*, '> Give input file name.'
  read(*,*) filnam
  open(file=filnam, status='old', unit=1)
  nstars=0
  do 
    read(1,*,iostat=iostat)
    if (iostat < 0) exit
    nstars=nstars+1
  end do
  allocate(star(nstars))

  rewind(1)
  do istar=1, nstars
    read(1,*) cols
    call zero_star(star(istar))
    star(istar)%col(1:ncol)%flg='OO'
    do icol=1, ncol
      read(cols(column_col(icol)),*) star(istar)%col(1)%data
      read(cols(column_uncer(icol)),*) star(istar)%col(1)%err
      if (column_id /= 0) read(cols(column_id),*) star(istar)%id
    end do
  end do

  call write_cluster_file(ncol, colstr, star, nstars)

end program ascii2cluster





