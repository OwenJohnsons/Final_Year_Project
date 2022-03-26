module define_star

  use put_header

  implicit none

  ! The maximum number of colours.  If you change this you also need to 
  ! change the format statement with a line number of 10.
  integer, parameter :: mcol=12

  type sexag_pos
    integer :: ra_h  = 0
    integer :: ra_m  = 0
    real    :: ra_s  = 0.0
    integer :: dc_d  = 0 
    integer :: dc_m  = 0
    real    :: dc_s  = 0.0
    character(len=1) :: dc_sign = '+'
  end type sexag_pos

  type a_colour
    real :: data = 0.0
    real :: err = 0.0
    logical :: neg_flux = .false.
    character(len=2) :: flg = 'AA'
  end type a_colour

  type a_star
    integer(kind=8) :: field = 0
    integer(kind=8) :: id    = 0
    integer :: ccd   = 0
    double precision :: ra = 0.0d0
    double precision :: dec = 0.0d0
    real :: x = 0.0
    real :: y = 0.0
    type(a_colour), dimension(mcol) :: col
  end type a_star

  double precision, parameter :: twopi=8.0d0*datan(1.0d0)
  integer, private, save, dimension(99) :: bytes_written, lines_written
  ! The format the file is being written.  A=ASCII, V=VOTable, F=FITS, N=Native.
  character(len=1), save :: file_format=' '

  ! For the match routines.
  type(a_star), dimension(:), allocatable, save, private :: stars
  integer, save, private :: nstars
  ! The range in RA and Dec of the catalogue to be matched.
  double precision, save, private :: min_alp, max_alp, min_del, max_del
  integer, dimension(:), allocatable, save, private :: ipoint

  logical, private :: pmt_nxtcls_out=.true., pmt_nxtcls_in=.true.
  character(len=100) :: nam_nxtcls_out, nam_nxtcls_in, nam_lstcls_in

  ! An array which contains the formats for each FITS field.
  character(len=8), private, allocatable, dimension(:) :: tform
  ! And an array of any offset.
  real, private, allocatable, dimension(:) :: tzero

  character(len=80), dimension(3), public, save :: cluster_header=' '
  ! Normally set to ' ', but if set, its copied into the header.
  character(len=80), save :: creation_date=' '

  interface zero_star

    module procedure zero_star_one
    module procedure zero_star_array

  end interface

contains

  subroutine write_star(iunit, star, ncol)

    ! Writes out a star. 

    integer, intent(in) :: iunit
    type(a_star), intent(in) :: star
    integer, optional :: ncol

    integer :: icol, jcol
    character(len=15) :: form
    type(a_star) :: st
    real :: flux, ferr
    double precision :: alpha, delta
    type(sexag_pos) :: sexag

    jcol=mcol
    if (present(ncol)) jcol=ncol

    inquire(unit=iunit, form=form)
    if (form == 'UNFORMATTED') then

      if (file_format == 'N') then

        ! This is the old unformatted format.
        call rad2sexag(star%ra, star%dec, sexag)
        write(iunit) star%field, star%ccd, star%id, sexag%ra_h, &
        sexag%ra_m, sexag%ra_s, sexag%dc_sign, sexag%dc_d, sexag%dc_m, &
        sexag%dc_s, star%x, star%y, (star%col(icol)%data, star%col(icol)%err, &
        star%col(icol)%flg,icol=1,jcol)

      else

        ! We write this out in an odd order to ensure that the first two 
        ! columns stand the best chance of being X and Y.
        if (ncol == 1) then
          write(iunit) clean_star_data(star%col(1))
        else 
          write(iunit) clean_star_data(star%col(2))
          write(iunit) clean_star_data(star%col(1))
          if (ncol > 2) then
            do icol=3, ncol
              write(iunit) clean_star_data(star%col(icol))
            end do
          end if
        end if

        write(iunit) swapped_i8(star%field)
        write(iunit) swapped_i4(star%ccd)
        write(iunit) swapped_i8(star%id)
        write(iunit) swapped_r8(360.0d0*star%ra/twopi)
        write(iunit) swapped_r8(360.0d0*star%dec/twopi)
        write(iunit) swapped_r4(star%x)
        write(iunit) swapped_r4(star%y)
      
        do icol=1, ncol
          write(iunit) swapped_r4(star%col(icol)%data)
          write(iunit) swapped_r4(star%col(icol)%err)
          write(iunit) star%col(icol)%flg
          write(iunit) boolean(star%col(icol)%neg_flux)
        end do

        ! 4 bytes each for ccd, xpos and ypos, plus 8 each for field, starid, RA
        ! and dec = 3*4 + 4*8 = 44.
        ! Each colour 3 reals, two characters and a logical = 3*4 + 2 + 1 = 15

        bytes_written(iunit)=bytes_written(iunit)+44+ncol*15
        lines_written(iunit)=lines_written(iunit)+1

      end if

    else

      st=star

      if (st%x<-9999.99 .or. st%x>99999.99 &
      .or. st%y<-9999.99 .or. st%y>99999.99) then
        st%x=0.0
        st%y=0.0
      end if

      ! Now negative fluxes, noted in the code by neg_flux.
      do icol=1, jcol
        if (st%col(icol)%neg_flux .and. st%col(icol)%flg(2:2)=='O') then
          ! These are flagged in the output file by M.
          st%col(icol)%flg(2:2)='M'
          ! First get the flux from the magnitude.
          flux=10.0**(-0.4*st%col(icol)%data)
          ! And now uncertainty in flux.
          ferr=10.0**((st%col(icol)%err-st%col(icol)%data)/2.5) - flux
          ! Now the horrible way these are encoded in cluster files.
          ! The magnitude is placed in the uncertainty column with the sign 
          ! changed.
          st%col(icol)%err=-st%col(icol)%data
          ! And the magnitude column is the flux uncertainty converted to a 
          ! magnitude.
          st%col(icol)%data=-2.5*log10(ferr)
        end if
        if (abs(st%col(icol)%data) > 999.999 .or. isnan(st%col(icol)%data)) then
          if (st%col(icol)%flg(1:1)/='A' .and. st%col(icol)%flg(2:2)/='A') then
            print*, 'Problem writing data for star ', st%field, st%ccd, st%id
            print*, 'Colour ', icol, ' flagged ', st%col(icol)%flg 
          end if
          st%col(icol)%data=0.0
        end if
        if (abs(st%col(icol)%err) > 999.999 .or. isnan(st%col(icol)%err)) then
          if (st%col(icol)%flg(1:1)/='A' .and. st%col(icol)%flg(2:2)/='A') then
            print*, 'Problem writing uncertainty for star ', &
            st%field, st%ccd, st%id
            print*, 'Colour ', icol, ' flagged ', st%col(icol)%flg 
          end if
          st%col(icol)%err=0.0
        end if
      end do
      
      ! Fit the star id, field number and CCD in the format statement.
      if (st%id    > 999999) st%id   =mod(st%id,    1000000)
      if (st%id    < -99999) st%id   =mod(st%id,    100000 )
      if (st%field >    999) st%field=mod(st%field, 1000   )
      if (st%ccd   >    100) st%ccd  =mod(st%ccd,   100    )

      call rad2sexag(star%ra, star%dec, sexag)

      write(iunit,10) st%field+real(st%ccd)/100.0, st%id, sexag%ra_h, &
      sexag%ra_m, sexag%ra_s, sexag%dc_sign, sexag%dc_d, sexag%dc_m, &
      sexag%dc_s, st%x, st%y, (st%col(icol)%data, st%col(icol)%err, &
      st%col(icol)%flg, icol=1,jcol)

    end if

10  format(1x,f6.2,2x,i6,2x,i2.2,1x,i2.2,1x,f6.3,1x, &
    a1,i2.2,1x,i2.2,1x,f5.2,2x, &
    2(f9.3,2x),12(f9.3,2x,f9.3,2x,a2))

  end subroutine write_star

  integer function read_star(iunit, star, ncol)

    integer, intent(in) :: iunit
    type(a_star), intent(inout) :: star
    integer, optional :: ncol

    integer :: icol, jcol, iostat, ifield
    character(len=4) :: dc_d
    real :: field_ccd
    character(len=15) :: form, access
    real :: ferr, flux
    integer :: i4
    integer*8 :: i8
    character(len=1) :: wrkbyte
    type(sexag_pos) :: sexag

    jcol=mcol
    if (present(ncol)) jcol=ncol

    star%col(1:ncol)%neg_flux=.false.

    inquire(unit=iunit, form=form, access=access)
    if (trim(form) == 'FORMATTED') then

      read(iunit,*,iostat=iostat) field_ccd, star%id, sexag%ra_h, &
      sexag%ra_m, sexag%ra_s, dc_d, sexag%dc_m, sexag%dc_s, star%x, &
      star%y, (star%col(icol)%data, star%col(icol)%err, &
      star%col(icol)%flg,icol=1,jcol)

      if (iostat == 0) then

        ! Sort out small negative declinations and convert into radians.
        read(dc_d,*) sexag%dc_d      
        sexag%dc_sign='+'
        if (dc_d(1:1) == '-') sexag%dc_sign='-'
        sexag%dc_d=abs(sexag%dc_d)
        call sexag2rad(sexag, star%ra, star%dec)

        ! Sort out the field and ccd numbers.
        star%field=int(field_ccd)
        if (100*star%field - nint(100.0*field_ccd) == 0) then
           star%ccd=0
        else
          star%ccd=nint(100.0*(field_ccd-real(star%field)))
        end if

        ! Convert any old-style flags.
        do icol=1, ncol
          if (star%col(icol)%flg(2:2) == ' ') then
            star%col(icol)%flg(2:2)=star%col(icol)%flg(1:1)
            star%col(icol)%flg(1:1)='O'
          end if
          call flagconv(star%col(icol)%flg)
        end do

      end if

    else if (trim(form) == 'UNFORMATTED') then

      if (trim(access) =='STREAM') then
        ifield=1
        do icol=1, ncol
          iostat=get_real(iunit, ifield, flux)
          if (iostat /= 0) goto 500
        end do
        iostat=get_long_integer(iunit, ifield, star%field)
        iostat=get_integer(iunit, ifield, star%ccd)
        iostat=get_long_integer(iunit, ifield, star%id)
        read(iunit, iostat=iostat) i8
        star%ra=unswap_r8(i8)
        star%ra  = twopi*star%ra/360.0d0
        read(iunit, iostat=iostat) i8
        star%dec=unswap_r8(i8)
        star%dec = twopi*star%dec/360.0d0
        ifield=ifield+2
        iostat=get_real(iunit, ifield, star%x)
        iostat=get_real(iunit, ifield, star%y)
        do icol=1, ncol
          iostat=get_real(iunit, ifield, star%col(icol)%data)
          iostat=get_real(iunit, ifield, star%col(icol)%err)
          read(iunit, iostat=iostat) star%col(icol)%flg
          read(iunit, iostat=iostat) wrkbyte
          if (wrkbyte=='T' .or. wrkbyte=='t' .or. wrkbyte=='1') then
            star%col(icol)%neg_flux=.true.
          else if (wrkbyte=='F' .or. wrkbyte=='f') then
            star%col(icol)%neg_flux=.false.
          else if (ichar(wrkbyte) == 0) then
            ! Technically this means undefined, but in TOPCAT it also seems to be false.
            star%col(icol)%neg_flux=.false.            
          else
            print*, 'Cannot intepret ', wrkbyte, ' as a Boolean.', ichar(wrkbyte)
            print*, 'Assuming false.'
            star%col(icol)%neg_flux=.false.                        
          end if
          ifield=ifield+2
        end do

      else

        read(iunit,iostat=iostat) star%field, star%ccd, star%id, &
        sexag%ra_h, sexag%ra_m, sexag%ra_s, sexag%dc_sign, sexag%dc_d, sexag%dc_m, &
        sexag%dc_s, star%x, star%y, (star%col(icol)%data, star%col(icol)%err, &
        star%col(icol)%flg,icol=1,jcol)
        call sexag2rad(sexag, star%ra, star%dec)

      end if

    else
      print*, 'Attempting to read file of form ', form
      stop
    end if

500 if (iostat == 0) then

      if (jcol < mcol) then
        do icol=jcol+1, mcol
          star%col(icol)%data=0.0
          star%col(icol)%err=0.0
          star%col(icol)%neg_flux=.false.
          star%col(icol)%flg='AA'
        end do
      end if

      do icol=1, ncol
        if (star%col(icol)%flg(2:2) == 'M') then
          ! This is where we sort out the ghastly way negative fluxes are
          ! are stored in cluster files.
          ! The data column stores the uncertainty as a magnitude, get
          ! it as a flux.
          ferr=10.0**(-0.4*star%col(icol)%data)
          ! The uncertainty stores the flux as a magnitude, sometimes with the
          ! sign negative.  So first get the flux from it.
          flux=10.0**(0.4*star%col(icol)%err)
          ! And now set the magnitude.
          star%col(icol)%data=-star%col(icol)%err
          ! And is a negative flux.
          star%col(icol)%neg_flux=.true.
          ! Which is unflagged.
          star%col(icol)%flg(2:2)='O'
          ! Finally create an uncertainty in magnitude space.
          !star%col(icol)%err = safe_M -2.5*log10((10.0**star%col(icol)%err)-1.0)
          star%col(icol)%err= 2.5*log10((ferr+flux)/flux)
        end if
      end do

    end if

    read_star=iostat

  end function read_star

  subroutine zero_star_array(star)

    type(a_star), intent(inout), dimension(:) :: star

    integer :: icol

    star%field=0
    star%ccd=0
    star%id=0
    star%ra=0.0d0
    star%dec=0.0d0
    star%x=0.0
    star%y=0.0
    do icol=1, mcol
      star%col(icol)%data=0.0
      star%col(icol)%err=0.0
      star%col(icol)%neg_flux=.false.
      star%col(icol)%flg='AA'
    end do

  end subroutine zero_star_array


  subroutine zero_star_one(star)

    type(a_star), intent(inout) :: star

    integer :: icol

    star%field=0
    star%ccd=0
    star%id=0
    star%ra=0.0d0
    star%dec=0.0d0
    star%x=0.0
    star%y=0.0
    do icol=1, mcol
      star%col(icol)%data=0.0
      star%col(icol)%err=0.0
      star%col(icol)%neg_flux=.false.
      star%col(icol)%flg='AA'
    end do

  end subroutine zero_star_one

  subroutine flux_to_mag(col)

    ! Changes the components of a colour from flux to magnitude.

    type(a_colour), intent(inout) :: col

    real :: flux, ferr, tinyflux

    ! The smallest absolute value of the flux that we can represent in 
    ! in magnitude space.  We never want a magnitude as faint as 
    ! 100th (since it will overrun our format 
    ! statements), which implies a flux greater than 10**-100/2.5. Nor do 
    ! we want a flux so near zero that when we take the log of it we get a 
    ! floating point error.
    tinyflux=max(100.0*tiny(1.0), 1.0e-40)

    flux=col%data
    ferr=col%err

    ! If this is a flux, then col%neg_flux should not be set.
    if (col%neg_flux) then
      print*, 'Found col%neg_flux set for a flux.'
      stop
    end if

    if (flux < 0.0) then
      col%neg_flux = .true.
    else
      col%neg_flux = .false.
    end if

    flux=abs(flux)
    flux=max(tinyflux, flux)
    ferr=max(tinyflux, ferr)

    col%data=-2.5*log10(flux)
    col%err = 2.5*log10((ferr+flux)/flux)

  end subroutine flux_to_mag


  subroutine mag_to_flux(col)

    ! Convert a magnitude into a flux, using the flag to allow
    ! for negative fluxes, and avoiding overflows for small fluxes.

    ! Note that col%neg_flux is always set false for a flux.

    type(a_colour), intent(inout) :: col

    real :: mag, err

    if (col%flg(2:2) =='M') then
      print*, 'Subroutine mag_to_flux has found a colour flagged M.'
      print*, 'these should be removed by using read_star.'
      stop
    end if

    mag=col%data
    err=col%err

    if (-mag/2.5 >= log10(huge(col%err)/10.0)) then 
      col%data=huge(col%err)/10.0
    else
      col%data=10.0**(-mag/2.5)
    end if
    ! Evaluating it this way avoids floating point overflows
    ! if err and mag are both large.  Such a condition happens
    ! when the flux is small.
    if ((err-mag)/2.5 >= log10(huge(col%err)/10.0)) then
      col%err=huge(col%err)/10.0
    else
      col%err= 10.0**((err-mag)/2.5) - col%data
    end if

    if (col%neg_flux) col%data=-1.0*col%data

    col%neg_flux=.false.

  end subroutine mag_to_flux

  subroutine flagconv(aflag)

    ! Converts the old numerical flags into the new character ones,
    ! assuming the old flags have been read as characters.

    character(len=2), intent(inout) :: aflag

    character, dimension(0:9) :: convert=(/'O', 'N', 'E', 'B', 'S', &
    'I', 'V', 'A', 'F', 'M'/)
    integer :: i, ichr
    character :: test

    do ichr=1, 2
      do i=0, 9
        write(test,'(i1)') i
        if (aflag(ichr:ichr) == test) aflag(ichr:ichr)=convert(i)
      end do
    end do

  end subroutine flagconv


  subroutine sort_star(star, nstars, sort, icol, ipoint)

    ! Sorts a star structure by various keys.

    ! The "increasing_col(....." keys are actually decreasing
    ! flux, thus allowing for negative fluxes.

    ! ipoint is the original position of the star in the list before
    ! sorting.

    type(a_star), dimension(:), intent(inout) :: star
    integer, intent(in) :: nstars
    character(len=*), intent(in) :: sort
    integer, optional :: icol
    integer, optional, dimension(:) :: ipoint

    integer :: i, j, l, ir, istar
    type(a_star) :: swap
    real, dimension(:), allocatable :: sort_key
    real :: sort_swap
    integer :: ipoint_swap
    
    type(a_colour), dimension(:), allocatable :: work

    if (nstars < 1) then

      print*, ' Error in s/r sort.  Number of stars is ', nstars
      stop

    else if (nstars > 1) then

      ! If nstars==1, then the data are already sorted!

      allocate(sort_key(nstars))
      if (sort == 'increasing_col(1)%data') then
        allocate(work(nstars))
        work=star%col(1)
        do istar=1, nstars
          call mag_to_flux(work(istar))
        end do
        sort_key=-1.0*work%data
      else if (sort == 'increasing_col(2)%data') then
        allocate(work(nstars))
        work=star%col(2)
        do istar=1, nstars
          call mag_to_flux(work(istar))
        end do
        sort_key=-1.0*work%data
      else if (sort == 'increasing_col(?)%data') then
        if (present(icol)) then
          allocate(work(nstars))
          work=star%col(icol)
          do istar=1, nstars
            call mag_to_flux(work(istar))
          end do
          sort_key=-1.0*work%data
        else
          print*, 'increasing_col(?)%data needs fourth argument.'
          stop
        end if
      else if (sort == 'increasing_ra') then
        sort_key = real(star%ra-star(1)%ra)
      else if (sort == 'increasing_dec') then
        sort_key = real(star%dec-star(1)%dec)
      else
        print*, 'Error in s/r sort.  Undefined sort key ', trim(sort)
        stop
      end if

      if (present(ipoint)) then
        do i=1, nstars
          ipoint(i)=i
        end do
      end if

      ir=nstars
      l=ir/2+1
      do
        if (l > 1)then
          l=l-1

          swap=star(l)
          sort_swap = sort_key(l)
          if (present(ipoint)) ipoint_swap=ipoint(l)

        else

          swap=star(ir)
          sort_swap=sort_key(ir)
          if (present(ipoint)) ipoint_swap=ipoint(ir)

          star(ir)=star(1)
          sort_key(ir)=sort_key(1)
          if (present(ipoint)) ipoint(ir)=ipoint(1)

          IR=IR-1
          if (ir == 1) then
            star(1)=swap
            sort_key(1)=sort_swap
            if (present(ipoint)) ipoint(1)=ipoint_swap
            exit
          end if
        end if
        i=l
        j=l+l
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(sort_key(j) < sort_key(j+1)) j=j+1
          ENDIF
          IF(sort_swap < sort_key(j))THEN

            star(i)=star(j)
            sort_key(i)=sort_key(j)
            if (present(ipoint)) ipoint(i)=ipoint(j)
            
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
          GO TO 20
        ENDIF

        star(i)=swap
        sort_key(i)=sort_swap
        if (present(ipoint)) ipoint(i)=ipoint_swap
        
      end do

      deallocate(sort_key)
      if (allocated(work)) deallocate(work)

    end if

  end subroutine sort_star

  subroutine nxtcls_in(ifname)

    character(len=*) :: ifname

    pmt_nxtcls_in=.false.
    nam_nxtcls_in=ifname

  end subroutine nxtcls_in


  integer function start_reading_cluster_file(iunit, nstars, ncol, colstr)

    ! Opens a cluster format file, and reads the header.
    ! The header ends up in cluster_header.

    ! The unit to be used (one day this should be made clever -- see below).
    integer, intent(out) :: iunit
    ! The number of stars, and the number of colours.
    integer, intent(out) :: nstars, ncol
    ! The names of the colours.
    character(len=*), dimension(:), intent(out) :: colstr

    logical :: there, done_header
    integer :: iostat, i, itrim, ipos, istart, iend
    character(len=10), dimension(mcol) :: unformatted_colstr
    character(len=100) :: name, binary_file

    character(len=200) :: work_header

    character(len=80), dimension(36) :: header
    character(len=20) :: ttype
    integer :: ihed, tfields, ifield, idummy

    ! Set iunit.
    iunit=21

    if (pmt_nxtcls_in) then
      get_name: do 
        print*, '> Give input cluster file name.'
        read(*,'(a)') nam_nxtcls_in
        inquire(file=nam_nxtcls_in, exist=there)
        if (there) exit get_name
        print*, 'Cannot find file ', trim(nam_nxtcls_in)
        print*, 'Try again.'
      end do get_name
    else
      pmt_nxtcls_in=.true.
    end if
    open(iunit,file=nam_nxtcls_in,status='old',iostat=iostat,form='formatted')
    start_reading_cluster_file=iostat
    if (iostat /= 0) return

    ! Three lines of header.
    ncol=-1
    read(iunit,*,iostat=iostat) ncol

    ! Now we can establish if this is an unformatted file or not.
    if (iostat/=0 .or. ncol<0 .or. ncol>mcol) then
      rewind(iunit)
      read(iunit, '(a200)') work_header
      if (work_header(1:6) =='SIMPLE') then
        print*, trim(nam_nxtcls_in), ' is a FITS table.'
        close(iunit)
        cluster_header(2)=' '
        open(iunit, file=nam_nxtcls_in,status='old',form='unformatted', access='stream')
        read(iunit) header
        ncol=0
        each_block: do
          read(iunit) header
          do ihed=1, 36
            if (trim(header(ihed)) == 'END') exit each_block
            if (header(ihed)(1:5) == 'TTYPE') then
              read(header(ihed)(11:30),*) ttype
              if (len_trim(ttype) > 2) then
                if (ttype(len_trim(ttype)-2:len_trim(ttype)) == 'MAG') then
                  ncol=ncol+1
                  colstr(ncol)=ttype(1:len_trim(ttype)-4)
                  do ipos=1, len(colstr(ncol))
                    if (colstr(ncol)(ipos:ipos) == '_') colstr(ncol)(ipos:ipos) = '-'
                  end do
                end if
              end if
            else if (header(ihed)(1:7) == 'TFIELDS') then
              read(header(ihed)(11:30),*) tfields
              if (allocated(tform)) deallocate(tform)
              if (allocated(tzero)) deallocate(tzero)
              allocate(tform(tfields), tzero(tfields))
              tzero=0.0
            else if (header(ihed)(1:6) == 'NAXIS2') then
              read(header(ihed)(11:30),*) nstars
            else if (header(ihed)(1:5) == 'TFORM') then
              read(header(ihed)(6:8),*) ifield
              read(header(ihed)(11:30),*) tform(ifield)
            else if (header(ihed)(1:5) == 'TZERO') then
              read(header(ihed)(6:8),*) ifield
              read(header(ihed)(11:30),*) tzero(ifield)
            else if (header(ihed)(1:8) == 'CLUSTER1') then
              read(header(ihed)(11:80),*) cluster_header(1)
            else if (header(ihed)(1:7)=='FILDATE' .or. trim(header(ihed)(1:8))=='DATE') then
              read(header(ihed)(11:80),*) cluster_header(2)
              cluster_header(2)='# File written at '//cluster_header(2)
            else if (header(ihed)(1:8) == 'DATE-CRE') then
              ! A data creation date, overwrite the header item.
              read(header(ihed)(11:80),*) cluster_header(2)
              cluster_header(2)='# Data created at '//cluster_header(2)
            else if (header(ihed)(1:8) == 'CLUSTER3') then
              read(header(ihed)(11:80),*) cluster_header(3)
            end if
          end do
        end do each_block
      else
        read(iunit, '(a200)') work_header
        if (work_header(2:8) == 'VOTABLE') then
          print*, trim(nam_nxtcls_in), ' is a VO table.'
          ! Find the magnitude and colour names.
          ncol=0
          done_header=.false.
          do
            read(iunit, '(a200)', iostat=iostat) work_header
            if (iostat < 0) exit
            if (work_header(1:6) == '<FIELD') then
              ipos=index(work_header, 'name')
              istart=index(work_header(  ipos:len(work_header)), '"')+ipos-1
              iend  =index(work_header(istart+1:len(work_header)), '"')+istart
              read(work_header(istart+1:iend-1),*) name
              if (len_trim(name) > 2) then
                if (name(len_trim(name)-2:len_trim(name)) == 'MAG') then
                  ncol=ncol+1
                  colstr(ncol)=name(1:len_trim(name)-4)
                  do ipos=1, len(colstr(ncol))
                    if (colstr(ncol)(ipos:ipos) == '_') colstr(ncol)(ipos:ipos) = '-'
                  end do
                end if
              end if
            else if (work_header(1:6) == '<TABLE') then 
              ipos=index(work_header, 'nrows')
              istart=index(work_header(  ipos:len(work_header)), '"')+ipos-1
              iend  =index(work_header(istart+1:len(work_header)), '"')+istart
              read(work_header(istart+1:iend-1),*) nstars
            else if (work_header(1:13) == '<STREAM href=') then 
              istart=index(work_header, '"')
              iend  =index(work_header(istart+1:len(work_header)), '"')+istart
              read(work_header(istart+1:iend-1),*) binary_file          
            else if (work_header(1:13)=='<DESCRIPTION>' .and. .not. done_header) then 
              istart=index(work_header, '>')
              iend  =index(work_header(istart+1:len(work_header)), '|')+istart
              cluster_header(1)=work_header(istart+1:iend-1)
              cluster_header(2)=' '
              istart=iend
              iend  =index(work_header(istart+1:len(work_header)), '<')+istart
              cluster_header(3)=work_header(istart+1:iend-1)            
              done_header=.true.
            end if
          end do
          close(iunit)
          ! Now, you may be accessing this file from a different directory, so find any UNIX
          ! directory name, and add it on to the binary file name.
          ipos=index(nam_nxtcls_in, '/', .true.)
          if (ipos > 0) binary_file=nam_nxtcls_in(1:ipos)//binary_file
          print*, 'Will read binary VO table data from ', trim(binary_file)
          open(unit=iunit, file=binary_file,  access='stream', form='unformatted', status='old')
        else
          ! Try unformatted.
          print*, 'Trying unformatted read for file ', trim(nam_nxtcls_in)
          close(iunit)
          open(iunit,file=nam_nxtcls_in,status='old',form='unformatted')
          read(iunit,iostat=iostat) cluster_header(2)
          if (iostat /= 0) then
            start_reading_cluster_file=iostat
            return
          end if
          read(iunit) cluster_header(1)
          read(iunit) cluster_header(3)
          read(iunit) ncol
          read(iunit) (unformatted_colstr(i),i=1,ncol)
          colstr(1:ncol)=unformatted_colstr(1:ncol)
          read(iunit) nstars
        end if
      end if
    else
      rewind(iunit)
      read(iunit,'(a200)') work_header
      read(work_header,*,iostat=iostat) ncol
      if (iostat < 0) ncol=4
      ! Now clear out any comments about the number of colours.
      itrim=index(work_header, '#')
      if (itrim < 1) then
        itrim=index(work_header, 'colours were created.')+21
      end if
      if (itrim < 0) itrim=0
      cluster_header(1)=work_header(itrim+1:len(work_header))
      ! The next line should have the names of the colours on it, but
      ! it may not.
      read(iunit,'(a80)') cluster_header(2)
      read(cluster_header(2),*,iostat=iostat) (colstr(i),i=1,ncol)
      if (iostat < 0) colstr=' '
      read(iunit,'(a80)') cluster_header(3)
      nstars=0
      count_them: do 
        read(iunit,*,iostat=iostat) idummy
        if (iostat < 0) exit count_them
        nstars=nstars+1
      end do count_them
      rewind(iunit)
      read(iunit,'(/,/)')
    end if

    nam_lstcls_in=nam_nxtcls_in

  end function start_reading_cluster_file

  subroutine lstcls_in(ifname)

    character(len=*), intent(out) :: ifname

    ifname=nam_lstcls_in

  end subroutine lstcls_in


  integer function read_cluster_file(nstars, ncol, colstr, star) 

    ! Reads in a cluster format file.

    integer, intent(out) :: nstars, ncol
    character(len=*), dimension(:), intent(out) :: colstr
    type(a_star), dimension(:), allocatable, intent(out) :: star
    type(a_star) :: test_star

    integer :: iostat, istar, iunit
    character(len=15) :: form

    read_cluster_file=0

    iostat=start_reading_cluster_file(iunit, nstars, ncol, colstr)

    if (iostat /= 0) then
      read_cluster_file=iostat
      return
    end if

    if (allocated(star)) deallocate(star)
    allocate(star(nstars))
    do istar=1, nstars
      iostat=read_star(iunit, star(istar), ncol)
      if (iostat /= 0) then
        print*, 'Error reading that file, iostat is', iostat
        print*, 'For star record ', istar
      end if
    end do

    close(iunit)

  end function read_cluster_file

  subroutine nxtcls_out(ifname)

    character(len=*) :: ifname

    pmt_nxtcls_out=.false.
    nam_nxtcls_out=ifname

  end subroutine nxtcls_out

  subroutine write_cluster_file(ncol, colstr, star, nstars)

    integer, intent(in) :: ncol
    character(len=*), dimension(:), intent(in) :: colstr
    type(a_star), dimension(:), intent(in) :: star
    integer, intent(in) :: nstars

    integer :: iunit, istar

    call start_cluster_file(ncol, colstr, iunit, nstars) 

    do istar=1, nstars
      call write_star(iunit, star(istar), ncol)
    end do

    call stop_cluster_file(iunit)

  end subroutine write_cluster_file

  subroutine start_cluster_file(ncol, colstr, iunit, nstars) 

    ! Opens a cluster file for writing and writes out the header.
    ! Prompts for a cluster file name if nxtcls_in has not been called.

    ! The file will be created as unformatted, if the optional argument
    ! nstars is there, and larger than the operating system variable 
    ! CLUSTER_FORMATTED_MAX_STARS, or set to zero.  In either of these 
    ! cases you must call stop_cluster_file to fill in the final padding
    ! and the number of lines in the header.

    ! Need to add the ability to write extra information in lines
    ! 1 and 3 (optional arguments).

    integer, intent(in) :: ncol
    character(len=*), dimension(:), intent(in) :: colstr
    integer, intent(out) :: iunit
    ! If nstars is not given, it will never write the file unformatted.
    integer, optional :: nstars

    ! Locals.
    logical :: open
    character(len=8) :: date
    character(len=10) :: time
    character(len=80) ::date_time
    integer :: i, iostat, icol, ipos, jpos
    character(len=80) :: outstr
    character(len=10), dimension(mcol) :: unformatted_colstr
    integer :: cfms
    character(len=50) :: wrkstr
    character(len=50), dimension(mcol) :: colstr_out
    character(len=100) :: filnam

    ! Get a free unit number.
    find: do i=12, 100
      if (i == 100) then
        print*,'s/r start_cluster_file cannot find a free unit number.'
        stop
      end if
      iunit=i
      inquire(unit=iunit, opened=open)
      if (.not. open) exit find
    end do find

    if (pmt_nxtcls_out) then
      print*, '> Give output file name for catalogue.'
      read(*,'(a)') nam_nxtcls_out
    else
      pmt_nxtcls_out=.true.
    end if
    call date_and_time(date, time)
    date_time=date(1:4)//'-'//date(5:6)//'-'//date(7:8) &
    //'T'//time(1:2)//':'//time(3:4)//':'//time(5:10)

    file_format='A'
    call getenv('CLUSTER_FORMATTED_MAX_STARS', outstr)
    read(outstr, *, iostat=iostat) cfms
    if (iostat /= 0) cfms=huge(cfms)
    if (present(nstars)) then
      if (nstars > cfms) then
        call getenv('CLUSTER_BINARY_FORMAT', outstr)
        if (trim(outstr) == 'V') then
          file_format='V'
        else if (trim(outstr) == 'N') then
          file_format='N'
        else
          file_format='F'
        end if
      end if
    end if

    if (file_format=='V') then

      open(unit=iunit, file=nam_nxtcls_out, form='formatted', status='replace')
      write(iunit,10) date(1:4), date(5:6), date(7:8), &
      time(1:2), time(3:4), time(5:10)
10    format("<?xml version='1.0'?>", /, &
      '<VOTABLE version="1.1"', /, &
      ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"', /, &
      ' xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.1 http://www.ivoa.net/xml/VOTable/v1.1"', /, &
      ' xmlns="http://www.ivoa.net/xml/VOTable/v1.1">', /, &
      '<!--', /, &
      ' !  VOTable written by ARK define_star', /, &
      ' !  at ', a4, '-', a2, '-', a2, 'T', a2, ':', a2, ':', a2, /, &
      ' !-->', /, &
      '<RESOURCE>')
      if (nstars > 0) then
        write(wrkstr,*) nstars
        write(iunit, '(a)') '<TABLE nrows="'//trim(adjustl(wrkstr))//'">'
      else
        write(iunit, '(a)') '<TABLE>'
      end if
      
      write(iunit, '(a)') '<DESCRIPTION>'//trim(cluster_header(1))//'|'// &
      trim(cluster_header(3))//'</DESCRIPTION>'

      ! Remove the minus sign from the colours.
      do icol=1, ncol
        colstr_out(icol)=colstr(icol)
        do ipos=1, len(colstr_out(icol))
          if (colstr_out(icol)(ipos:ipos) == '-') colstr_out(icol)(ipos:ipos) = '_'
        end do
      end do

      if (ncol == 1) then
        write(iunit,20) trim(colstr_out(1)), trim(colstr_out(1))
      else
        write(iunit,20) trim(colstr_out(2)), trim(colstr_out(2))
        write(iunit,20) trim(colstr_out(1)), trim(colstr_out(1))
        if (ncol > 2) then
          do icol=3, ncol
            write(iunit,20) trim(colstr_out(icol)), trim(colstr_out(icol))
          end do
        end if
      end if
20    format('<FIELD datatype="float" name="', a, '_MAG_CLEAN">', /, &
      '<DESCRIPTION>Cleaned ', a, ' magnitude sample (excluded values null)</DESCRIPTION>', /, '</FIELD>') 

      write(iunit,120)
120   format( &
      '<FIELD datatype="int" name="FIELD">', /, &
      '<DESCRIPTION>Field number</DESCRIPTION>', /, &
      '</FIELD>', /, &
      '<FIELD datatype="int" name="CCD">', /, &
      '<DESCRIPTION>CCD number</DESCRIPTION>', /, &
      '</FIELD>', /, &
      '<FIELD datatype="int" name="STAR_ID">', /, &
      '<DESCRIPTION>Star identification number</DESCRIPTION>', /, &
      '</FIELD>', /, &
      '<FIELD datatype="double" name="RA" ucd="POS_EQ_RA_MAIN" unit="degrees">', /, &
      '<DESCRIPTION>Right Ascension</DESCRIPTION>', /, &
      '</FIELD>', /, &
      '<FIELD datatype="double" name="DEC" ucd="POS_EQ_DEC_MAIN" unit="degrees">', /, &
      '<DESCRIPTION>Declination</DESCRIPTION>', /, &
      '</FIELD>', /, &
      '<FIELD datatype="float" name="XPOS">', /, &
      '<DESCRIPTION>X position on CCD</DESCRIPTION>', /, &
      '</FIELD>', /, &
      '<FIELD datatype="float" name="YPOS">', /, &
      '<DESCRIPTION>Y position on CCD</DESCRIPTION>', /, &
      '</FIELD>')
      
      do icol=1, ncol
        write(iunit,50) trim(colstr_out(icol)), trim(colstr_out(icol))
        write(iunit,70) trim(colstr_out(icol)), trim(colstr_out(icol))
        write(iunit,90) trim(colstr_out(icol)), trim(colstr_out(icol))
        write(iunit,14) trim(colstr_out(icol)), trim(colstr_out(icol))
      end do

50    format('<FIELD datatype="float" name="', a, '_MAG">', /, &
      '<DESCRIPTION>', a, ' magnitude</DESCRIPTION>', /, '</FIELD>')
70    format('<FIELD datatype="float" name="', a, '_UNCERT">', /, &
      '<DESCRIPTION>', a, ' magnitude uncertainty</DESCRIPTION>', /, '</FIELD>')
90    format('<FIELD arraysize="2" datatype="char" name="', a, '_FLAG">', /, &
      '<DESCRIPTION>', a, ' magnitude quality flag</DESCRIPTION>', /, '</FIELD>')
14    format('<FIELD datatype="boolean" name="', a, '_NEG_FLUX">', /, &
      '<DESCRIPTION>Is ', a, ' flux negative (and magnitude from its absolute value)</DESCRIPTION>', /, '</FIELD>')
      
      write(iunit,100)
100   format('<DATA>', /, '<BINARY>')
      ipos=index(nam_nxtcls_out, '.', .true.)
      if (ipos == 0) ipos=len(trim(nam_nxtcls_out))+1
      ! Now, the problem with these VO files having reference to the data file is that you
      ! might not be accessing this from the directory in which it lives.  So, remove any
      ! remove any UNIX directory name.
      jpos=index(nam_nxtcls_out, '/', .true.)
      write(iunit,210) trim(nam_nxtcls_out(jpos+1:ipos-1)//'.bin')
210   format('<STREAM href="', a, '"/>')
      write(iunit,220) 
220   format('</BINARY>', /, '</DATA>', /, '</TABLE>', /, '</RESOURCE>', /, '</VOTABLE>')

      close(iunit)
      print*, 'Will write the VO table binary data to ', trim(nam_nxtcls_out(1:ipos-1)//'.bin')
      open(unit=iunit, file=trim(nam_nxtcls_out(1:ipos-1)//'.bin'), access='stream', &
      form='unformatted', status='replace')

    else if (file_format == 'F') then

      call write_cluster_fits_header(iunit, nam_nxtcls_out, ncol, colstr, &
      nstars, date_time, cluster_header(1), cluster_header(3))
      bytes_written(iunit)=0
      lines_written(iunit)=0

    else if (file_format == 'N') then

      ! This is the old unformatted form.
      open(unit=iunit, file=nam_nxtcls_out, status='unknown', &
      form='unformatted')
      outstr=' '
      write(iunit) ' # File written at '//date_time
      write(iunit) cluster_header(1)
      write(iunit) cluster_header(3)
      write(iunit) ncol
      unformatted_colstr(1:ncol)=colstr(1:ncol)
      write(iunit) unformatted_colstr(1:ncol)
      write(iunit) nstars

    else

      open(unit=iunit, file=nam_nxtcls_out, status='replace')!, recl=66+ncol*24)
      write(iunit,*) ncol, 'colours were created. # ', cluster_header(1)
      if (len_trim(creation_date) > 0) then
        write(iunit,*) (trim(colstr(i))//' ',i=1,ncol), ' # Data created at '//trim(creation_date)
      else
        write(iunit,*) (trim(colstr(i))//' ',i=1,ncol), ' # File written at '//date_time
      end if
      write(iunit,*) cluster_header(3)

    end if

  end subroutine start_cluster_file

  subroutine write_cluster_fits_header(iunit, ofname, ncol, colstr, nstars, &
  date, clust_hed1, clust_hed3)

    character(len=*), intent(in) :: ofname, date, clust_hed1, clust_hed3
    integer, intent(in) :: iunit, ncol, nstars
    character(len=*), dimension(:) :: colstr
    character(len=10), dimension(mcol) :: colstr_out

    integer :: icol, ipos

    ! Write the first header block.'
    call clear_header()
    call put_header_c('SIMPLE', 'T', 'Standard FITS format', 1)
    call put_header_i('BITPIX',   8, 'Character data', 1)
    call put_header_i('NAXIS',    0, 'No image, just extensions', 1)
    call put_header_c('EXTEND', 'T', 'There are standard extensions', 1)
    call put_header_s('COMMENT', 'Written the ARK FITS routines', ' ', 1)
    open(unit=iunit, file=ofname, form='unformatted', access='stream', &
    status='replace')
    write(iunit) dross(1:36)

    ! Now the real header.
    call clear_header()
    call put_header_s('XTENSION',  'BINTABLE', 'Binary table extension',    1)
    call put_header_i('BITPIX',             8, '8-bit bytes',               1)
    call put_header_i('NAXIS',              2, '2-dimensional table',       1)
    call put_header_i('NAXIS1',    44+ncol*15, 'Width of table in bytes',   1)
    call put_header_i('NAXIS2',        nstars, 'Number of rows in table',   1)
    call put_header_i('PCOUNT',             0, 'Size of special data area', 1)
    call put_header_i('GCOUNT',             1, 'One data group',            1)
    call put_header_i('TFIELDS',     7+5*ncol, 'number of columns',         1)
    call put_header_s('EXTNAME', trim(ofname), 'table name',                1)
    call put_header_s('DATE',      trim(date), 'Date/time file written',    1)
    if (len_trim(creation_date) > 0) then
      call put_header_s('DATE-CRE', trim(creation_date), &
      'Date/time data were created',    1)
      creation_date=' '
    else
      call put_header_s('DATE-CRE', trim(date), &
      'Date/time data were created',    1)
    end if
    call put_header_s('CLUSTER1',  clust_hed1, ' ',                         1)
    call put_header_s('CLUSTER3',  clust_hed3, ' ',                         1)

    ! Remove the minus sign from the colours.
    do icol=1, ncol
      colstr_out(icol)=colstr(icol)
      do ipos=1, len(colstr_out(icol))
        if (colstr_out(icol)(ipos:ipos) == '-') &
        colstr_out(icol)(ipos:ipos) = '_'
      end do
      if (len_trim(colstr(icol)) == 0) then
        colstr_out(icol)=''
      else
        colstr_out(icol)=trim(colstr_out(icol))//'_'
      end if
    end do

    if (ncol == 1) then
      call put_header_column(trim(colstr_out(1))//'MAG_CLEAN', 'E', &
      'Cleaned '//trim(colstr(1))//' magnitude sample (excluded values null)')
    else
      call put_header_column(trim(colstr_out(2))//'MAG_CLEAN', 'E', &
      'Cleaned '//trim(colstr(2))//' magnitude sample (excluded values null)')
      call put_header_column(trim(colstr_out(1))//'MAG_CLEAN', 'E', &
      'Cleaned '//trim(colstr(1))//' magnitude sample (excluded values null)')
      if (ncol > 2) then
        do icol=3, ncol
          call put_header_column(trim(colstr_out(icol))//'MAG_CLEAN', 'E', &
          'Cleaned '//trim(colstr(icol))//&
          ' magnitude sample (excluded values null)')
        end do
      end if
    end if

    call put_header_column('FIELD',   'K', 'Field number')
    call put_header_column('CCD  ',   'J', 'Detector number')
    call put_header_column('STAR_ID', 'K', 'Star identification number')
    call put_header_column('RA',      'D', 'Right Ascension')
    call put_header_column('DEC',     'D', 'Declination')
    call put_header_column('XPOS',    'E', 'X position on detector')
    call put_header_column('YPOS',    'E', 'Y position on detector')

    do icol=1, ncol
      call put_header_column(trim(colstr_out(icol))//'MAG',   'E', &
      trim(colstr(icol))//' magnitude')
      call put_header_column(trim(colstr_out(icol))//'UNCERT','E', &
      trim(colstr(icol))//' magnitude uncertainty')
      call put_header_column(trim(colstr_out(icol))//'FLAG',  '2A', &
      trim(colstr(icol))//' magnitude quality flag')
      call put_header_column(trim(colstr_out(icol))//'NEG_FLUX',  'L', &
      'Is '//trim(colstr(icol))//&
      ' flux negative (and magnitude from its absolute value)')
    end do
    if (mod(lndrss(), 36) == 0) then
      write(iunit) dross(1:lndrss())
    else
      write(iunit) dross(1:36*(1+(lndrss()/36)))
    end if

  end subroutine write_cluster_fits_header
    

  subroutine stop_cluster_file(iunit)

    ! Closes a cluster file opened for writing.  
    ! Useful for padding out FITS file.

    integer, intent(in) :: iunit
    character(len=100) :: ofname
    integer :: ibyte
    character(len=80), dimension(36) :: work

    if (file_format == 'F') then
      ! Pad out the remaining bytes for a FITS file.
      if (mod(bytes_written(iunit),2880) > 0) then
        do ibyte=mod(bytes_written(iunit),2880)+1, 2880
          write(iunit) ' '
        end do
      end if

      inquire(unit=iunit, name=ofname)
      close(iunit)

      ! Now fix the length of the file in the header.
      open(unit=iunit, file=ofname, access='direct', recl=2880)
      read(iunit, rec=2) work
      write(work(5)(10:30), '(i21)') 0
      write(work(5)(10:30), '(i21)') lines_written(iunit)
      write(iunit, rec=2) work

    end if

    close(iunit)

  end subroutine stop_cluster_file


  integer function append_cluster_file(filnam, iunit)

    ! Allows you to append to a binary (or ASCII) cluster format
    ! file.  Note that for binary files the number of star records
    ! has to have been be set for the final number that will be in there.

    implicit none

    character(len=*) :: filnam
    integer, intent(in) :: iunit

    integer :: iostat, istart, iend
    character(len=200) :: work_header, binary_file
    character :: test_byte

    open(iunit, file=filnam, status='old', iostat=iostat, form='formatted')
    append_cluster_file=iostat
    if (iostat /= 0) return

    ! Three lines of header.
    read(iunit, '(a200)') work_header
    if (work_header(1:6) =='SIMPLE') then
      close(iunit)
      open(iunit, file=filnam, status='old', iostat=iostat, access='stream', &
      form='unformatted')
      bytes_written(iunit)=0
      find_length: do 
        read(iunit,iostat=iostat) test_byte
        if (iostat < 0) exit find_length
        bytes_written(iunit)=bytes_written(iunit)+1
      end do find_length
      close(iunit)
      print*, 'Opening FITS file ', trim(filnam), ' which has ', bytes_written(iunit), &
      ' bytes so far.'
      open(iunit, file=filnam, status='old', access='stream', &
      form='unformatted', position='append')
    else 
      read(iunit, '(a200)') work_header
      if (work_header(2:8) == 'VOTABLE') then
        ! This is a VOTABLE, so need to know the stream name.
        open_binary: do 
          read(iunit, '(a200)') work_header
          if (work_header(1:13) == '<STREAM href=') then 
            istart=index(work_header, '"')
            iend  =index(work_header(istart+1:len(work_header)), '"')+istart
            read(work_header(istart+1:iend-1),*) binary_file          
            close(iunit)
            open(iunit, file=binary_file, status='old', iostat=iostat, access='stream', &
            form='unformatted', position='append')
            append_cluster_file=iostat
            if (iostat /= 0) return
            exit open_binary
          end if
        end do open_binary
      else
        close(iunit)
        open(iunit, file=filnam, status='old', iostat=iostat, form='formatted', position='append')
      end if

    end if
    
  end function append_cluster_file

  subroutine sexag2rad(star, alpha, delta)

    type(sexag_pos), intent(in) :: star
    double precision, intent(out) :: alpha, delta

    alpha=twopi*(dble(star%ra_h)+dble(star%ra_m)/60.0d0+&
    dble(star%ra_s)/3600.0d0)/24.0d0

    ! The abs are just in case some cluster programme still signals
    ! negative dec via the signs of the components.
    delta=twopi*(real(abs(star%dc_d))+real(abs(star%dc_m))/60.0d0+&
    dble(abs(star%dc_s))/3600.0d0)/360.0d0

    if (star%dc_sign == '-') delta=-1.0d0*delta

  end subroutine sexag2rad

  subroutine rad2sexag(ra_rad, dec_rad, star)

    ! Changes from radians to sexagesimal with three decimal places of RA seconds
    ! and two of declination seconds, to match cluster formatted output.

    double precision, intent(in) :: ra_rad, dec_rad
    type(sexag_pos), intent(inout) :: star

    double precision :: ra, dec

    ra=24.0d0*ra_rad/twopi
    ! Round to three significant figures in RA seconds.
    ra=dble(nint(1000.0d0*3600.0d0*ra))/(1000.0d0*3600.0d0)
    star%ra_h=int(ra)
    ra=60.0d0*(ra-dble(star%ra_h))
    ! And check that the above multiplication has not wrecked our rounding.
    ra=dble(nint(1000.0d0*60.0d0*ra))/(1000.0d0*60.0d0)
    star%ra_m=int(ra)
    star%ra_s=real(60.0d0*(ra-dble(star%ra_m)))
    
    dec=360.0d0*dec_rad/twopi
    star%dc_sign='+'
    if (dec < 0.0d0) then
      star%dc_sign='-'
      dec=abs(dec)
    end if
    ! Round to two significant figures in Dec seconds.
    dec=dble(nint(100.0d0*3600.0d0*dec))/(100.0d0*3600.0d0)
    star%dc_d=int(dec)
    dec=60.0d0*(dec-dble(star%dc_d))
    ! And check that the above multiplication has not wrecked our rounding.
    dec=dble(nint(100.0d0*60.0d0*dec))/(100.0d0*60.0d0)
    star%dc_m=int(dec)
    star%dc_s=real(60.0d0*(dec-dble(star%dc_m)))

  end subroutine rad2sexag

  subroutine set_cluster_header(cluster_header_in)

    character(len=*), dimension(3), intent(in) :: cluster_header_in

    cluster_header=cluster_header_in

  end subroutine set_cluster_header

  integer function get_real(iunit, ifield, data)

    ! Reads the next field, and puts it into real, whether or not it was written
    ! as a real.  Uses tform to see how many bytes to read.

    integer, intent(in) :: iunit
    integer, intent(inout) :: ifield
    real, intent(out) :: data

    integer :: bytes, iostat
    character :: i1
    integer :: i4
    integer*8 :: i8

    bytes=4

    if (allocated(tform)) then
      if (trim(tform(ifield)) == 'B') then
        bytes=1
      else if (trim(tform(ifield)) == 'E') then
        bytes=4
      else if (trim(tform(ifield)) == 'D') then
        bytes=8
      else
        print*, 'Cannot cope with form ', trim(tform(ifield)), ' for field ', ifield
        stop
      end if
    end if

    
    if (bytes == 1) then
      read(iunit, iostat=iostat) i1
      data=real(ichar(i1))
    else if (bytes == 4) then 
      read(iunit, iostat=iostat) i4
      data=unswap_r4(i4)
    else if (bytes == 8) then
      read(iunit, iostat=iostat) i8
      data=unswap_r8(i8)
    end if

    if (allocated(tzero)) data=data+tzero(ifield)
    ifield=ifield+1
    get_real=iostat

  end function get_real

  integer function get_integer(iunit, ifield, data)

    ! Reads the next field, and puts it into an integer, whether or not it was written
    ! as a 4-byte integer.  Uses tform to see how many bytes to read.

    integer, intent(in) :: iunit
    integer, intent(inout) :: ifield
    integer, intent(out) :: data

    integer :: bytes, iostat
    character :: i1
    integer*2 :: i2
    integer :: i4
    integer*8 :: i8

    bytes=4

    if (allocated(tform)) then
      if (trim(tform(ifield)) == 'B') then
        bytes=1
      else if (trim(tform(ifield)) == 'I') then
        bytes=2
      else if (trim(tform(ifield)) == 'J') then
        bytes=4
      else if (trim(tform(ifield)) == 'K') then
        bytes=8
      else
        print*, 'Cannot cope with form ', trim(tform(ifield)), ' for field ', ifield
        stop
      end if
    end if

    if (bytes == 1) then
      read(iunit, iostat=iostat) i1
      data=ichar(i1)
    else if (bytes == 2) then
      read(iunit, iostat=iostat) i2
      data=swapped_i2(i2)
    else if (bytes == 4) then
      read(iunit, iostat=iostat) i4
      data=swapped_i4(i4)
    else if (bytes == 8) then
      read(iunit, iostat=iostat) i8
      data=swapped_i8(i8)
    end if

    if (allocated(tzero)) data=data+nint(tzero(ifield))
    ifield=ifield+1
    get_integer=iostat

  end function get_integer

  integer function get_long_integer(iunit, ifield, data)

    ! Reads the next field, and puts it into a long integer, whether or not it was written
    ! as a 8-byte integer.  Uses tform to see how many bytes to read.

    integer, intent(in) :: iunit
    integer, intent(inout) :: ifield
    integer(kind=8), intent(out) :: data

    integer :: bytes, iostat
    character :: i1
    integer*2 :: i2
    integer :: i4
    integer*8 :: i8

    bytes=8

    if (allocated(tform)) then
      if (trim(tform(ifield)) == 'B') then
        bytes=1
      else if (trim(tform(ifield)) == 'I') then
        bytes=2
      else if (trim(tform(ifield)) == 'J') then
        bytes=4
      else if (trim(tform(ifield)) == 'K') then
        bytes=8
      else
        print*, 'Cannot cope with form ', trim(tform(ifield)), ' for field ', ifield
        stop
      end if
    end if

    if (bytes == 1) then
      read(iunit, iostat=iostat) i1
      data=ichar(i1)
    else if (bytes == 2) then
      read(iunit, iostat=iostat) i2
      data=swapped_i2(i2)
    else if (bytes == 4) then
      read(iunit, iostat=iostat) i4
      data=swapped_i4(i4)
    else if (bytes == 8) then
      read(iunit, iostat=iostat) i8
      data=swapped_i8(i8)
    end if

    if (allocated(tzero)) data=data+nint(tzero(ifield))
    ifield=ifield+1
    get_long_integer=iostat

  end function get_long_integer

  integer function swapped_r4(r4)
    real, intent(in) :: r4
    integer :: i
    do i=0, 24, 8
      call mvbits(transfer(r4, 0_4),  i, 8, swapped_r4, 24-i)
    end do
  end function swapped_r4

  integer*8 function swapped_r8(r8)
    double precision, intent(in) :: r8
    integer :: i
    do i=0, 56, 8
      call mvbits(transfer(r8, 0_8),  i, 8, swapped_r8, 56-i)
    end do
  end function swapped_r8
  
  integer function swapped_i2(i2)
    integer*2, intent(in) :: i2
    integer :: i
    do i=0, 16, 8
      call mvbits(transfer(i2, 1),  16-i, 8, swapped_i2, i)
    end do
  end function swapped_i2

  integer function swapped_i4(i4)
    integer, intent(in) :: i4
    integer :: i
    do i=0, 24, 8
      call mvbits(transfer(i4, 1),  24-i, 8, swapped_i4, i)
    end do
  end function swapped_i4

  integer*8 function swapped_i8(i8)
    integer*8, intent(in) :: i8
    integer :: i
    do i=0, 56, 8
      call mvbits(transfer(i8, 0_8),  56-i, 8, swapped_i8, i)
    end do
  end function swapped_i8

  real function unswap_r4(i4)
    integer, intent(in) :: i4
    integer :: i, work
    do i=0, 24, 8
      call mvbits(i4,  i, 8, work, 24-i)
      unswap_r4=transfer(work, 0.0)
    end do
  end function unswap_r4

  double precision function unswap_r8(i8)
    integer*8, intent(in) :: i8
    integer :: i
    integer*8 :: work
    do i=0, 56, 8
      call mvbits(i8,  i, 8, work, 56-i)
    end do
    unswap_r8=transfer(work, 0.0d0)
  end function unswap_r8
  
  character(len=1) function boolean(logical)
    logical, intent(in) :: logical
    if (logical) then
      boolean='T'
    else
      boolean='F'
    end if
  end function boolean

  integer function clean_star_data(col)
    type(a_colour), intent(in) :: col
    real :: work=-1.0
    if (col%err<0.1 .and. col%flg=='OO' .and. .not. col%neg_flux) then
      clean_star_data = swapped_r4(col%data)
    else
      clean_star_data = swapped_r4(sqrt(work))
    end if
  end function clean_star_data

  subroutine bracket(array, value, min, max)

    ! This routine will give you the numbers of the array elements
    ! which bracket the value.  The values of the bracketing elements
    ! may be equal to the value, or above or below.

    double precision, dimension(:), intent(in) :: array
    double precision, intent(in) :: value
    integer, intent(out) :: min, max

    integer :: iclose, test

    min=1
    max=size(array,1)
    iclose=max-min
    do
      test=(max+min)/2
      if (value >= array(test)) then
        min=test
      else if(value <= array(test)) then
        max=test
      end if
      if (max-min == iclose) exit
      iclose=max-min
    end do

  end subroutine bracket

  subroutine init_match(stars_in, nstars_in)

    ! This routine prepares the catalogue to be matches against, by 
    ! storing a sorted version of it.

    type(a_star), dimension(:), intent(in) :: stars_in
    integer, intent(in) :: nstars_in

    integer :: istars

    nstars=nstars_in
    if (allocated(stars)) deallocate(stars)
    if (allocated(ipoint)) deallocate(ipoint)
    allocate(stars(nstars), ipoint(nstars))
    stars=stars_in
    call sort_star(stars, nstars, 'increasing_dec', 0, ipoint)

    min_alp=minval(stars%ra)
    max_alp=maxval(stars%ra)
    min_del=stars(1)%dec
    max_del=stars(nstars)%dec

    !print*, 'Range of RA is ', min_alp, max_alp
    !print*, 'Range of Dec is ', min_del, max_del

  end subroutine init_match

  subroutine faster_match(box, rad, matches, n_match, jpoint)

    ! The star record for the central position.
    type(a_star), intent(in) :: box
    ! The radius in radians
    real, intent(in) :: rad
    type(a_star), dimension(:), intent(out) :: matches
    integer, intent(out) :: n_match
    integer, intent(out), optional, dimension(:) :: jpoint

    if (present(jpoint)) then
      call quicker_match(box%ra, box%dec, rad, matches, n_match, jpoint)
    else
      call quicker_match(box%ra, box%dec, rad, matches, n_match)
    end if

  end subroutine faster_match

  subroutine quicker_match(box_alp, box_del, rad, matches, n_match, jpoint)

    ! This routine matches the error circle box, against the stored
    ! catalogue, using a radius (in radians) of rad.  There are n_match
    ! stars which fall within this radius returned in the array matches.

    ! jpoint is the ordering of the matching stars as they were supplied to
    ! init_match (so may, for example) reflect their magnitude.

    double precision, intent(in) :: box_alp, box_del
    real, intent(in) :: rad
    type(a_star), dimension(:), intent(out) :: matches
    integer, intent(out) :: n_match
    integer, intent(out), optional, dimension(:) :: jpoint

    integer :: istars, i_match
    real :: delta_del, delta_alp
    integer :: min, max, idummy

    n_match=0
    if (box_alp+rad >= min_alp) then
      if (box_alp-rad <= max_alp) then
        if (box_del+rad >= min_del) then
          if (box_del-rad <= max_del) then
            ! Find a bracketting element from low dec.
            call bracket(stars%dec, box_del-rad, min, idummy)
            ! And from above.
            call bracket(stars%dec, box_del+rad, idummy, max)
            each_star: do istars=min, max
              delta_del=real(box_del-stars(istars)%dec)
              if (abs(delta_del) > rad) cycle each_star
              delta_alp=real((box_alp-stars(istars)%ra)*cos(stars(istars)%dec))
              if (abs(delta_alp) > rad) cycle each_star
              if (delta_del*delta_del+delta_alp*delta_alp < rad*rad) then
                n_match=n_match+1
                if (n_match > size(matches,1)) then
                  print*, 'Found more than ', n_match-1, ' matches.  The array'
                  print*, 'matches passed to faster_match needs expanding.'
                  print*, 'Error box at ', box_alp, box_del
                  print*, 'Radius ', rad*57.3*3600.0, ' arcsec.'
                  do i_match=1, n_match-1
                    print*, matches(i_match)%field, matches(i_match)%ccd, &
                    matches(i_match)%id, &
                    matches(i_match)%ra, matches(i_match)%dec
                  end do
                  stop
                end if
                matches(n_match)=stars(istars)
                if (present(jpoint)) jpoint(n_match)=ipoint(istars)
              end if
            end do each_star
          end if
        end if
      end if
    end if

  end subroutine quicker_match

  subroutine stand_cor(alpha, delta, t_alpha, t_delta, q, psi, eta)

    ! Converts from RA and DEC to tanget plane co-ordinates.

    implicit none

    ! The RA and dec to be transformed.
    double precision, intent(in) :: alpha, delta
    ! The tangent point and pincushion/barrel q.
    double precision, intent(in) :: t_alpha, t_delta
    real, intent(in) :: q
    real, intent(out) :: psi, eta

    real :: work

    work=sin(t_delta)*sin(delta)+cos(t_delta)*cos(delta)*cos(alpha-t_alpha)

    ! Long standing bug in this line (t_alpha was written as t_delta)
    ! cured by Timn October 2002.
    psi=cos(delta)*sin(alpha-t_alpha)/work

    eta=cos(t_delta)*sin(delta)-sin(t_delta)*cos(delta)*cos(alpha-t_alpha)

    eta=eta/work

    ! And now allow for the pincushion/barrel.
    work=1.0+q*(psi**2.0 + eta**2.0)
    psi=work*psi
    eta=work*eta

  end subroutine stand_cor

end module define_star
