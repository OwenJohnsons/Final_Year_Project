module bolcor_readfiles

  use quad
  use cmdfit_system

  implicit none

  ! The maximum number of gravities and temperatures.
  integer, parameter, public :: m_grav=20, m_teff=500
  ! Define a structure for the bolometric corrections which makes sure the indices clear.
  type a_teff
    real, dimension(m_grav) :: logg
  end type a_teff
  type a_corr
    type(a_teff), dimension(m_teff) :: teff
  end type a_corr

  ! And we need an array of flags to go with it.
  type flag_teff
    character(len=2), dimension(m_grav) :: logg
  end type flag_teff
  type flag_corr
    type(flag_teff), dimension(m_teff) :: teff
  end type flag_corr

contains

  subroutine sort_grav(tau_data)
    
    real, intent(inout), dimension(:)  :: tau_data

    integer :: k, l, m
    real :: safe

    outer: do k=2, size(tau_data)
      safe=tau_data(k)
      do l=1,k-1
        if (safe < tau_data(l)) then
          do m=k, l+1, -1
            tau_data(m)=tau_data(m-1)
          end do
          tau_data(l)=safe
          cycle outer
        endif
      end do
    end do outer

  end subroutine sort_grav

  subroutine readfiles(file, n_bc, iteff, ilogg, bc, bc_flag, temp, grav, &
  n_teff, n_grav)


    ! Allows table of bolometric corrections to be read in.
    ! It assumes the first column is effective temperature, then gravity (if
    ! available), then the BCs.

    character(len=*), intent(in) :: file
    integer, intent(in) :: n_bc
    integer, intent(in) :: iteff, ilogg
    type(a_corr), intent(inout), dimension(:), allocatable :: bc
    type(flag_corr), intent(inout), dimension(:), allocatable :: bc_flag
    real, dimension(m_teff), intent(inout) :: temp
    real, dimension(m_grav), intent(inout) :: grav
    integer, intent(out) :: n_teff, n_grav

    character(len=80), save :: lastfile=' '
    
    integer :: i, j, int_grav, int_temp, iostat, j_grav
    real, dimension(n_bc)  :: bc_in
    real :: grav_in, temp_in, temp_safe

    real, dimension(m_grav) :: col_grav, work_grav
    integer :: k_grav
    character(len=20) :: flag, iflg
    integer :: ncol
    character(len=10), dimension(:), allocatable :: colnames

    !if (trim(lastfile) == trim(file)) return

    lastfile=file

    ! Find the number of different temperatures and gravities
    call read_head(file, 10, ncol, colnames)
    n_teff=1
    n_grav=0
    if (ilogg > 0) n_grav=1
    read(10,*) temp_safe, grav(1)
    do
      read(10,*,iostat=iostat) temp_in, grav_in
      if (iostat < 0) exit

      if (nint(temp_in) /= nint(temp_safe) ) then
        n_teff=n_teff+1
        temp_safe=temp_in
      end if

      if (ilogg > 0) then
        if (minval(abs(grav(1:n_grav)-grav_in)) > 0.1) then
          n_grav=n_grav+1
          if (n_grav > m_grav) then
            print*, 'm_grav needs expanding in readfiles.'
            stop
          end if
          grav(n_grav)=grav_in
        end if
      end if
      
    end do

    if (n_teff > m_teff) then
      print*, 'In module bolcor_readfiles m_teff needs enlarging to ', n_teff
      stop
    end if

    !print*, 'Found ', n_teff, ' effective temperatures.'

    if (ilogg > 0) then
      ! Sort the gravities into the (traditional) decreasing numeric order.
      grav=-1.0*grav
      call sort_grav(grav(1:n_grav))
      grav=-1.0*grav
      !print*, 'Found ', n_grav, 'gravaties.' 
    end if

    close(10)
    call read_head(file, 10, ncol, colnames)
    read(10,*) temp(1)
    i=1
    do
      read(10,*) temp_in
      if (nint(temp_in) /= nint(temp(i))) then
        i=i+1
        temp(i)=temp_in
      end if
      if (i == n_teff) exit
    end do

    ! Sort the temperatures into the (traditional) decreasing numeric order.
    call sort_grav(temp(1:n_teff))

    temp=log10(temp)

    allocate(bc(n_bc), bc_flag(n_bc))

    ! Set the flags as unknown.
    do k_grav=1, m_grav
      do int_temp=1, m_teff
        bc_flag(1:n_bc)%teff(int_temp)%logg(k_grav)='AA'
      end do
    end do

    close(10)
    call read_head(file, 10, ncol, colnames)
    do 
      if (ilogg > 0) then
        read(10,*,iostat=iostat) temp_in, grav_in, bc_in
      else
        read(10,*,iostat=iostat) temp_in, bc_in
      end if
      if (iostat < 0) exit
      ! We work in log temperature.
      temp_in=log10(temp_in)
      ! Find the array indices most appropriate for the read temperature 
      ! and gravity.
      if (ilogg > 0) then
        !int_grav=minloc(abs(grav-grav_in),1)
        int_grav=locate_nearest(grav(1:n_grav), grav_in)
        if (abs(grav(int_grav)-grav_in) > 0.01) then
          print*, 'Cannot find logg of ', grav_in
          print*, 'In the array ', grav
          stop
        end if
      else
        int_grav=1
      end if
      !int_temp=minloc(abs(temp-temp_in),1)
      int_temp=locate_nearest(temp(1:n_teff), temp_in)
      if (abs(temp(int_temp)-temp_in) > 0.001) then
        print*, 'Cannot find Teff of ', temp_in
        print*, 'In the array ', temp(1:n_teff)
        stop
      end if
      ! And fill those elements of the flag arrary which are OK.
      bc(1:n_bc)%teff(int_temp)%logg(int_grav)=bc_in(1:n_bc)
      where(bc_in(1:n_bc) > -98.0) bc_flag(1:n_bc)%teff(int_temp)%logg(int_grav)='OK'
    end do

    close(10)

    !open(unit=22, file='junk.dat')
    !do int_temp=1, n_teff
    !  do int_grav=1, n_grav
    !    write(22,*) temp(int_temp), grav(int_grav), bc(3)%teff(int_temp)%logg(int_grav), &
    !    bc_flag%teff(int_temp)%logg(int_grav)
    !  end do
    !end do
    !stop

  end subroutine readfiles

  type(a_corr) function subtract_corr(bc1, bc2, n_teff, n_grav)

    type(a_corr), intent(in) :: bc1, bc2
    integer, intent(in) :: n_teff, n_grav

    integer :: i_teff, i_grav

    subtract_corr = bc1
    do i_teff=1, n_teff
      do i_grav=1, n_grav
        subtract_corr%teff(i_teff)%logg(i_grav) = &
        bc1%teff(i_teff)%logg(i_grav) - bc2%teff(i_teff)%logg(i_grav)
      end do
    end do

  end function subtract_corr

  real function interp(grav_in, temp_in, col, bc_flag, grav, temp, n_teff, n_grav, flg)

    real, intent(in) :: grav_in, temp_in
    type(a_corr), intent(in) :: col
    type(flag_corr), intent(in) :: bc_flag
    real, dimension(:), intent(in) :: grav, temp
    integer, intent(in) :: n_teff, n_grav
    character(len=*), intent(out) :: flg

    real, dimension(m_grav) :: col_grav, work_grav
    real, dimension(n_teff) :: col_temp, work_temp
    character(len=20) :: flag
    character(len=2), dimension(n_grav) :: grav_flag

    integer :: j_grav, k_temp, j_temp

    flg='OK'

    ! In the below code I have tried using linear (as opposed to quadratic) interpolation 
    ! in both gravity and effective temperature
    ! with the Kurucz atmospheres and the Bessell & Murphy system responses.  These are
    ! probably the poorest sampled in terms of temperature and so are a good test.
    ! I find the difference can be as much as 3 hundreths, so its best to stick with
    ! quadratic interpolation.
    k_temp=0
    grav_flag='OK'
    do j_temp=1, n_teff
      col_temp(k_temp+1)=quadintflg( &
      !col_temp(k_temp+1)=linintflg( &
      grav(1:n_grav), col%teff(j_temp)%logg(1:n_grav), grav_in, flag, &
      grav_flag(1:n_grav), bc_flag%teff(j_temp)%logg(1:n_grav))
      if (flag == 'OK') then
        k_temp=k_temp+1
        work_temp(k_temp)=temp(j_temp)
      end if
    end do
    if (k_temp > 1) then
      interp=quadint(work_temp(1:k_temp), col_temp(1:k_temp), temp_in, flag)
      !interp=linint(work_temp(1:k_temp), col_temp(1:k_temp), temp_in, flag)
      if (flag == '1') flg='(grav)'
    else
      flg='1'
      interp=0.0
    end if

  end function interp

end module bolcor_readfiles


module bolcor_subs

  use bolcor_readfiles
  use cmdfit_system

  implicit none

contains

  integer function name_match(list, name)

    character(len=*), intent(in) :: name
    character(len=*), dimension(:), intent(in) :: list

    integer :: i

    name_match=0

    do i=1, size(list,1)
      if (trim(name) == trim(list(i))) then
        name_match=i
        exit
      end if
    end do

  end function name_match


  logical function find_bc_combination(colour, names, n_names, icol1, icol2)

    character(len=*), intent(in) :: colour
    character(len=*), dimension(:) :: names
    integer, intent(in) :: n_names
    integer, intent(out), optional :: icol1, icol2

    integer :: iminus, iname1, iname2

    find_bc_combination=.false.
    if (present(icol1)) icol1=0
    if (present(icol2)) icol2=0

    iname1 = name_match(names, trim(colour))
    ! Is this a simple magnitude, if so we are done.
    if (iname1>0 .and. iname1<=n_names) then
      iname2=0
    else
      ! If not maybe we can do it by subtracting two BCs.

      iminus=index(colour, '-')

      if (iminus == 0) then
        print*, 'Colour ', trim(colour), ' is not available.'
        return
      end if

      iname1 = name_match(names, trim(colour(1:iminus-1)))
      if (iname1<1 .or. iname1>n_names) then
        print*, 'Colour ', trim(colour), ' is not available.'
        return
      end if

      iname2 = name_match(names, trim(colour(iminus+1:len(colour))))
      if (iname2<1 .or. iname2>n_names) then
        print*, 'Colour ', trim(colour), ' is not available.'
        return
      end if

    end if

    find_bc_combination=.true.
    if (present(icol1)) icol1=iname1
    if (present(icol2)) icol2=iname2

  end function find_bc_combination

  subroutine bolcor(iso, m_bol_sun)

    ! Takes the isochrone and fills in the colours at each point 
    ! i.e. iso%pnt(1:iso%npts)%col(1:iso%ncol).

    type(chrone), intent(inout) :: iso
    real, intent(out) :: m_bol_sun

    character(len=10), allocatable, dimension(:) :: s_column
    integer :: iminus, n_column, i_column, icol, ipt

    integer :: n_grav, n_teff
    type(a_corr), dimension(:), allocatable :: bc
    type(flag_corr), dimension(:), allocatable :: bc_flag
    real, dimension(m_teff) :: temp
    real, dimension(m_grav) :: grav
    character(len=2), dimension(:), allocatable :: work_flag
    type(flag_corr) :: work_bc_flag
    type(cols) :: colnums
    logical :: first=.true.
    real :: work
    integer :: iostat

    if (first) then
      print*, 'Using bolometric correction file ', trim(iso%bc_file)
      first=.false.
    end if

    call read_head(iso%bc_file, 2, n_column, s_column, colnums)
    if (colnums%iteff == 0) print*, 'Could not find effective temperature in BC file.'
    if (colnums%ilogg == 0) print*, 'Could not find logG in BC file.'

    iostat=read_header_real(2, 'm_bol_sun', m_bol_sun)
    if (iostat < 0) then
      m_bol_sun=4.74
      print*, ' '
      print*, 'Failed to find bolometric magnitude of the Sun in the head of file'
      print*, trim(iso%bc_file)
      print*, 'Using value of ', m_bol_sun
      print*, ' '
    end if

    close(2)

    iso%ncol=0
    do i_column=1, n_column
      if (s_column(i_column)(1:3) == 'BC_') then
        iminus=index(s_column(i_column), '-')
        iso%ncol=iso%ncol+1
        if (iminus == 0) then
          ! Simply remove the BC at the beginning.
          iso%colnam(iso%ncol)=s_column(i_column)(4:len(s_column(i_column)))
        else
          ! Recall that colours are negative of BCs.
          iso%colnam(iso%ncol)=s_column(i_column)&
          (iminus+1:len_trim(s_column(i_column)))//'-'//s_column(i_column)(4:iminus-1)
        end if
      end if
    end do
    
    call readfiles(iso%bc_file, iso%ncol, colnums%iteff, &
    colnums%ilogg, bc, bc_flag, temp, grav, n_teff, n_grav)
    
    ! Now for each colour find the bolometric correction.
    do ipt=1, iso%npts
      do icol=1, iso%ncol
        if (n_grav > 0) then
          ! Easy.
          iso%pnt(ipt)%col(icol)%data=interp(iso%logg(ipt), iso%teff(ipt), bc(icol), &
          bc_flag(icol), grav, temp, n_teff, n_grav, iso%pnt(ipt)%col(icol)%flag)
        else
          ! Even easier.
          iso%pnt(ipt)%col(icol)%data=quadintflg(temp(1:n_teff), &
          bc(icol)%teff(1:n_teff)%logg(1), iso%teff(ipt), & 
          iso%pnt(ipt)%col(icol)%flag, yflg=bc_flag(icol)%teff(1:n_teff)%logg(1))
        end if
      end do
    end do
    

  end subroutine bolcor


end module bolcor_subs
