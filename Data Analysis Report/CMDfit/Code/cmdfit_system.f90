module cmdfit_system

  ! Use this module for stuff which calls nothing else, but may be
  ! used by many other modules.

  implicit none

  ! The numbers for the standard columns which may appear.
  type cols
    integer :: imass, ilogg, ilbol, ilogte, iteff, iage
  end type cols

  ! Maximum number of colours in an isochrone.
  integer, parameter :: iso_mcol=200

  ! We are heading towards a structure which is as follows.
  ! iso%iso_file
  ! iso%pnt(ipt)%mass
  ! iso%pnt(ipt)%col(icol)%data
  ! iso%pnt(ipt)%col(icol)%flag
  type a_col
    real :: data
    character(len=50) :: flag
  end type a_col
  
  type a_pnt
    type(a_col), dimension(iso_mcol) :: col 
  end type a_pnt
  
  type chrone 
    ! Now the specification of the iscohrone.  Its age, colours and model.
    real :: age
    character(len=10), dimension(2) :: colreq
    character(len=10), dimension(iso_mcol) :: colnam
    character(len=200) :: iso_file, bc_file
    integer :: npts, ncol=2
    type(a_pnt), dimension(:), allocatable :: pnt
    real, dimension(:), allocatable :: mass, teff, logg, lbol
    character(len=50), dimension(:), allocatable :: &
    massflg, teffflg, loggflg, lbolflg
  end type chrone

  ! The minimun difference in log age allowable.  This is needed because
  ! taking the log of the age can result in considerable rounding error.
  ! Its set by the (curently 3) places after the decimal in the age in the
  ! file names.
  real, parameter :: age_diff_min=0.001

contains

  character(len=100) function data_dir()
    
    ! A little function to return where the datafiles are stored.
    character(len=100) :: here
    integer :: length
    
    call get_environment_variable('CMDDATA', here, length)

    if (length > len(here)) then
      print*, 'The environment variable CMDDATA is of length ', length
      print*, 'But the string to save it is of length ', len(here)
      stop
    end if

    data_dir=trim(here)//'/'
    
  end function data_dir

  logical function bell_there(ext_file)

    character(len=*) :: ext_file

    bell_there=.false.
    if (len_trim(ext_file) > 2) then
      if (ext_file(len_trim(ext_file)-2:len_trim(ext_file)) == 'ext') bell_there=.true.
    end if

  end function bell_there
  
  integer function count_header(iunit)

    integer, intent(in) :: iunit
    character(len=50) :: instring
    
    integer :: ihead

    ihead=0
    count_head: do 
      read(iunit,*) instring
      if (instring(1:1) == '#') then
        ihead=ihead+1
      else
        exit count_head 
      end if
    end do count_head

    count_header=ihead

  end function count_header
  
  integer function read_header_real(iunit, name, value)

    integer, intent(in) :: iunit
    character(len=*), intent(in) :: name
    real, intent(out) :: value

    character(len=500) :: instring
    character(len=50) :: name_in
    integer :: iend, iostat

    read_header_real=-1
    rewind(iunit)
    do
      read(iunit, '(a500)', iostat=iostat) instring
      if (iostat < 0) exit
      instring=adjustl(instring)
      if (instring(1:1) == '#') then
        read(instring(2:len(instring)),*,iostat=iostat) name_in
        if (iostat == 0) then
          if (name_in == name) then
            iend=index(instring, trim(name_in))
            iend=iend+len_trim(name_in)
            read(instring(iend:len(instring)),*) value
            read_header_real=0
          end if
        end if
      end if
    end do

  end function read_header_real
    


  subroutine read_head(ifname, iunit, ncol, colnames, colnums)

    ! Open a file, read the column headers, and find a group of standard names.

    character(len=*), intent(in) :: ifname
    integer, intent(in) :: iunit
    integer, intent(out) :: ncol
    character(len=*), dimension(:), allocatable, intent(out) :: colnames
    type(cols), intent(out), optional :: colnums

    character(len=500) :: instring, wkstring
    integer :: ihead, i, ipos, ihed

    inquire(file=ifname, number=i)

    open(iunit, file=ifname, action='read') 

    ihead=count_header(iunit)

    if (ihead > 0) then
      rewind(iunit)
      do i=1, ihead-1
        read(iunit,*)
      end do
      read(iunit,'(a500)') instring
      instring=adjustl(instring)
      ! Remove the tabs.
      do i=1, len(instring)
        if (instring(i:i) == char(9)) instring(i:i)=' '
      end do

      ncol=0
      wkstring=instring(2:len(instring))
      find_ncol: do
        ipos=index(trim(wkstring), ' ', .true.)
        if (ipos > 0) then
          ncol=ncol+1
          wkstring=wkstring(1:ipos-1)
        else
          exit find_ncol
        end if
      end do find_ncol
      
      allocate(colnames(ncol))

      ihed=ncol+1
      wkstring=instring(2:len(instring))
      get_names: do
        ipos=index(trim(wkstring), ' ', .true.)
        if (ipos > 0) then
          ihed=ihed-1
          colnames(ihed)=adjustl(trim(wkstring(ipos+1:len(wkstring))))
          wkstring=wkstring(1:ipos-1)
        else
          exit get_names
        end if
      end do get_names
        
      if (present(colnums)) then
        ! Now get the names of the column heads.
        colnums%imass=0 
        colnums%ilogg=0
        colnums%ilbol=0 
        colnums%ilogte=0
        colnums%iteff=0
        colnums%iage=0
        each_header: do ihed=1, ncol
          if (trim(colnames(ihed)) == 'Mini') colnums%imass=ihed
          if (trim(colnames(ihed)) == 'LogG') colnums%ilogg=ihed
          if (trim(colnames(ihed)) == 'logG') colnums%ilogg=ihed
          if (trim(colnames(ihed)) == 'logg') colnums%ilogg=ihed
          if (trim(colnames(ihed)) == 'logL/Lo') colnums%ilbol=ihed
          if (trim(colnames(ihed)) == 'logTe') colnums%ilogte=ihed
          if (trim(colnames(ihed)) == 'Teff') colnums%iteff=ihed
          if (trim(colnames(ihed)) == 'log(age/yr)') colnums%iage=ihed
        end do each_header

      end if

    end if

  end subroutine read_head

  integer function read_line(iunit, ncol, data)
    
    integer, intent(in) :: iunit, ncol
    real, dimension(:), intent(out) :: data

    integer :: iostat
    character(len=500) :: test_str
    
    each_line: do
      read(iunit,'(a500)',iostat=iostat) test_str
      read_line=iostat
      test_str=adjustl(test_str)
      if (iostat == 0) then
        if (test_str(1:1) == '#') cycle each_line
        read(test_str,*) data(1:ncol)
      end if
      exit each_line
    end do each_line
    
  end function read_line

  subroutine sort(xdata, ydata, ndata)

    implicit none

    real, dimension(:), intent(inout) :: xdata, ydata
    integer, intent(in) :: ndata

    integer :: i, j, l, ir
    real :: xswap, yswap

    if (ndata < 1) then
      print*, ' Error in s/r sort.  Number of data points is ', ndata
      stop
    else if (ndata > 1) then
      ! If ndata==1, then the data are already sorted!
      ir=ndata
      l=ir/2+1
      do
        if (l > 1)then
          l=l-1
          xswap=xdata(l)
          yswap=ydata(l)
        else
          xswap=xdata(ir)
          yswap=ydata(ir)
          xdata(ir)=xdata(1)
          ydata(ir)=ydata(1)
          IR=IR-1
          if (ir == 1) then
            xdata(1)=xswap
            ydata(1)=yswap
            exit
          end if
        end if
        i=l
        j=l+l
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(xdata(j) < xdata(j+1)) j=j+1
          ENDIF
          IF(xswap < xdata(j))THEN
            xdata(i)=xdata(j)
            ydata(i)=ydata(j)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
          GO TO 20
        ENDIF
        xdata(i)=xswap
        ydata(i)=yswap
      end do

    end if

  end subroutine sort

end module cmdfit_system
