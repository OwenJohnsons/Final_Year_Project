module red_ext_interp

  use quad
  use cmdfit_system
  use bolcor_subs

  implicit none

  ! The number of extinctions, filters and gravities covered
  integer, parameter :: n_ext=5

  ! A structure to hold elements from the extinction tables
  type ref_ext_entry
     real :: ref_col
     real, dimension(:,:), allocatable :: ref_ext
     real, dimension(n_ext) :: ebmv
     real, dimension(n_ext) :: ebmv_val
  end type ref_ext_entry

  type spec_ext_values
     integer :: n_entries
     real, dimension(:), allocatable :: temp ! Specific effective temperature
     real, dimension(:,:), allocatable :: red ! Specific reddening
     real, dimension(:,:), allocatable :: ext ! Specific extinction
     real, dimension(n_ext) :: ebmv_val
     real :: logg
  end type spec_ext_values

  type entire_ext_values
    integer :: n_entries
    real, dimension(:), allocatable :: teff ! Effective temperatures
    real, dimension(:,:,:), allocatable :: exts ! Entire library of extinctions
    real, dimension(n_ext) :: ebmv_val
    real :: logg
  end type entire_ext_values

contains


  subroutine get_nom_ebmv(real_ebmv_in, nom_ebmv_out, photon, extrap, file)

    ! Size of array and define size of arrays which hold E(B-V) info
    integer, parameter :: n_array = 11
    real, dimension(n_array) :: nom_ebmv, act_ebmv

    ! Input measured E(B-V) and integration method
    real, intent(in) :: real_ebmv_in
    logical, intent(in) :: photon

    ! Output nominal E(B-V), have we extrapolated for this value?
    real, intent(out) :: nom_ebmv_out
    logical, intent(out) :: extrap

    character(len=*) :: file
    character(len=50) :: string

    ! Flag
    character(len=2) :: flag

    ! Assume that we have not extrapolated
    extrap = .false.

    ! Is the file already in hot star E(B-V)?
    open(unit=1, file=file, status='old')
    read_hed: do
      read(1,'(a)') string
      if (string(1:1) /= '#') exit read_hed
      if (index(string, 'Hot star') > 0) then
        nom_ebmv_out=real_ebmv_in
        close(1)
        return
      end if
    end do read_hed
    close(1)
        
    print*, 'Converting hot star E(B-V) to nominal E(B-V).'

    ! Simply state the values, save reading in all the time.
    ! Nominal values
    nom_ebmv(1) = 0.0
    nom_ebmv(2) = 0.2
    nom_ebmv(3) = 0.4
    nom_ebmv(4) = 0.6
    nom_ebmv(5) = 0.8
    nom_ebmv(6) = 1.0
    nom_ebmv(7) = 1.2
    nom_ebmv(8) = 1.4
    nom_ebmv(9) = 1.6
    nom_ebmv(10) = 1.8
    nom_ebmv(11) = 2.0

    ! Real values derived using the Kurucz ODFnew atmoshperic models
    ! in conjunction with the Bessell UxBVRI filters.
    ! E(B-V)_nom vs E(B-V)_real
    if(photon) then
      act_ebmv(1) = 0.0
      act_ebmv(2) = 0.19925366
      act_ebmv(3) = 0.39542684
      act_ebmv(4) = 0.58838926
      act_ebmv(5) = 0.77802991
      act_ebmv(6) = 0.9642553
      act_ebmv(7) = 1.1469898
      act_ebmv(8) = 1.32620025
      act_ebmv(9) = 1.50184575
      act_ebmv(10) = 1.67393535
      act_ebmv(11) = 1.8425
    else
      act_ebmv(1) = 0.0
      act_ebmv(2) = 0.2019442
      act_ebmv(3) = 0.40094522
      act_ebmv(4) = 0.59685932
      act_ebmv(5) = 0.78955541
      act_ebmv(6) = 0.97891975
      act_ebmv(7) = 1.16485705
      act_ebmv(8) = 1.3472928
      act_ebmv(9) = 1.52618905
      act_ebmv(10) = 1.7015029
      act_ebmv(11) = 1.8732434
    end if

    ! Interpolate to get the nominal E(B-V) corresponding to the given E(B-V)
    ! Extrap flag F for interpolation and T for extrapolation
    nom_ebmv_out=linintflg(act_ebmv, nom_ebmv, real_ebmv_in, flag)
    if (trim(flag) == '1') extrap = .true.

  end subroutine get_nom_ebmv

!!!!!!!!!!!!!

  subroutine bell_ext(ext_file, colnam, nom_ebmv, logg_in, teff_in, int_col, int_mag, extrap)

    ! This routine reddens stars according to Cameron's tables.

    ! The primary input is the nominal E(B-V), which is the correponding
    ! value derived for a measured E(B-V) assumed to be for a star of
    ! colour zero. Difference in E(B-V) between this colour zero star (Teff~10000 K)
    ! and the bluest star in model atmoshpere grid (Teff~50000 K) is only 0.03
    ! and so will make little difference in resulting isochrone.

    ! But, to make the extinction absolutely correct you also need to supply
    ! the effective temperature of the atmosphere, and its gravity.

    ! Input nominal E(B-V), stellar gravity, intrinsic colour and magnitude
    character(len=*), intent(in) :: ext_file
    character(len=*), dimension(2), intent(in) :: colnam
    real, intent(in) :: nom_ebmv, logg_in

    ! Check if first time/run
    logical, save :: first=.true.
    character(len=10), dimension(2), save :: prev_colnam =(/'          ', '          '/)

    real :: int_col, int_mag, teff_in

    ! Output reddening and extinction
    ! Have we extrapolated for these values
    real :: red_out, ext_out    
    logical, intent(out) :: extrap

    ! Counters
    integer :: itemp, ilogg

    character(len=2) :: flag

    ! Structures to hold table entries
    type(entire_ext_values), allocatable, dimension(:), save :: all_ext
    type(spec_ext_values), allocatable, dimension(:), save :: logg_ext
    integer, save :: n_logg

    ! For working out colour numbers.
    integer :: n_bc, i_column, i_bc, iminus, iname1, iname2, icol
    integer, save :: n_column
    character(len=10), dimension(:), allocatable :: s_column, names
    real, dimension(:,:), allocatable :: work

    ! If this is the first iteration, read in the extinction data
    if (first) then
      print*, 'Reading extinctions from file ', trim(ext_file)
      call read_head(ext_file, 2, n_column, s_column)
      close(2)
      ! How many of these columns are extinctions?
      n_bc=0
      do i_column=1, n_column
        if (s_column(i_column)(1:2) == 'A_') n_bc=n_bc+1
      end do
      allocate(names(n_bc))
      i_bc=0
      do i_column=1, n_column
        if (s_column(i_column)(1:2) == 'A_') then
          iminus=index(s_column(i_column), '-')
          i_bc=i_bc+1
          if (iminus == 0) then
            ! Simply remove the BC at the beginning.
            names(i_bc)=s_column(i_column)(3:len(s_column(i_column)))
          else
            ! Recall that colours are negative of BCs.
            names(i_bc)=s_column(i_column)(iminus+1:len_trim(s_column(i_column)))//'-'//&
            s_column(i_column)(3:iminus-1)
          end if
        end if
      end do

      call read_extin_info(ext_file, all_ext, logg_ext, n_logg, n_bc)

      do ilogg=1, n_logg
        logg_ext(ilogg)%temp(:)=all_ext(ilogg)%teff(:)
      end do
 
      first=.false.
    end if

    ! Due to bolometric correction nomenclature, the filter order is counter intuitive
    ! i.e. (g-i) is filter 4(i) - filter 2(g)
    do icol=1, 2
      if (trim(prev_colnam(icol)) .ne. trim(colnam(icol))) then
        if (.not. find_bc_combination(colnam(icol), names, n_bc, iname1, iname2)) then
          print*, 'Programming error led to bell_ext being called in non-existant colour ', &
          trim(colnam(icol))
          print*, 'Amoungst the extinctions ', names(1:n_bc)
          stop
        end if
        allocate(work(size(all_ext(1)%exts,1), size(all_ext(1)%exts,2)))
        do ilogg=1, n_logg
          work=all_ext(ilogg)%exts(:,:,iname1)
          if (iname2 /= 0) work=work-all_ext(ilogg)%exts(:,:,iname2)
          if (icol == 1) then
            logg_ext(ilogg)%ext(:,:)=work
          else
            logg_ext(ilogg)%red(:,:)=work
          end if
          deallocate(work)
        end do
        prev_colnam(icol)=colnam(icol)
      end if
    end do

    ! Now interpolate for the required logg and colour, and the required
    ! nominal E(B-V) to find a resulting colour excess and extinction
    ! This requires a treble interpolation-ooh scary!
    call treble_interp(logg_in, int_col, teff_in, nom_ebmv, &
         n_logg, logg_ext, ext_out, red_out, extrap)

    ! Calculate the reddened colour and extincted magnitude.
    int_col=int_col+red_out
    int_mag=int_mag+ext_out

    !print*, 'Resultant reddening and extinction are', red_out, ext_out
    !print*, 'Final colour and magnitude are', int_col, int_mag

!    print*, col_out, (mag_out+8.72)

  end subroutine bell_ext

!!!!!!!!!!!!!

  subroutine read_extin_info(ext_file, all_ext, logg_ext, n_logg, n_filter)

    ! This subroutine reads all of the extinciton information in

    character(len=*), intent(in) :: ext_file

    ! Structures to hold the information
    type(entire_ext_values), allocatable, dimension(:) :: all_ext
    type(spec_ext_values), allocatable, dimension(:) :: logg_ext
    ! Outputs
    integer, intent(out) :: n_logg
    ! Input
    integer, intent(in) :: n_filter

    ! Counters
    integer :: iext, ilogg, ios, ifilt, ientry

    ! Path names
    character(len=1000) :: path, file, prefix, suffix
    character(len=4), allocatable, dimension(:) :: string_logg
    integer :: ntemp

    real, dimension(512) :: grav_list, temp_list
    real :: dummy
    real :: temp_in, logg_in, ebmv_in
    real, dimension(:), allocatable :: bc_in
    character(len=1) :: test
    integer :: ihead, nhead

    open (unit=22, file=ext_file)
    nhead=0
    count_head: do 
      read(22,*) test
      if (test /= '#') exit count_head 
      nhead=nhead+1
    end do count_head
    rewind(22)

    do ihead=1, nhead
      read(22,*)
    end do
    ! Now find the number of gravities.
    read(22,*) dummy, grav_list(1)
    n_logg=1
    count_grav: do
      read(22,*,iostat=ios) dummy, grav_list(n_logg+1)
      if (ios < 0) exit count_grav
      if (minval(abs(grav_list(1:n_logg)-grav_list(n_logg+1))) > 0.01) n_logg=n_logg+1
      if (n_logg >= size(grav_list,1)) then
        print*, 'grav_list needs expanding.'
        stop
      end if
    end do count_grav

    ! Sort them.
    call sort_grav(grav_list(1:n_logg))
    if (allocated(all_ext)) deallocate(all_ext)
    if (allocated(logg_ext)) deallocate(logg_ext)
    allocate(all_ext(n_logg), logg_ext(n_logg))
    all_ext%logg = grav_list(1:n_logg)
    logg_ext%logg = grav_list(1:n_logg)
      
    ! Read in extinction information.
    do ilogg=1, n_logg
       ! Skip the header
       rewind(22)
       do ihead=1, nhead
         read(22,*)
       end do
       ntemp=0
       count_do: do 
         if (ntemp >= size(temp_list,1)) then
           print*, 'temp_list needs expanding.'
           stop
         end if
          read(22,*,iostat=ios) temp_list(ntemp+1), dummy
          if (ios < 0) exit count_do
          if (abs(dummy - grav_list(ilogg)) > 0.01) cycle count_do
          if (ntemp == 0) then
            ntemp=ntemp+1
          else if (minval(abs(temp_list(1:ntemp)-temp_list(ntemp+1))) > 0.01) then
            ntemp=ntemp+1
          end if
       end do count_do

       call sort_grav(temp_list(1:ntemp))
       all_ext(ilogg)%n_entries=ntemp
       logg_ext(ilogg)%n_entries=ntemp
       allocate(all_ext(ilogg)%teff(ntemp), all_ext(ilogg)%exts(ntemp,n_ext,n_filter))
       allocate(logg_ext(ilogg)%temp(ntemp), logg_ext(ilogg)%red(ntemp,n_ext))
       allocate(logg_ext(ilogg)%ext(ntemp,n_ext))
       all_ext(ilogg)%teff=temp_list(1:ntemp)
       logg_ext(ilogg)%temp=temp_list(1:ntemp)
    end do

    ! Now we hard code in the E(B-V) values for each set of repeated columns
    do iext=1, n_ext
       all_ext(:)%ebmv_val(iext) = (real(iext)/2.0)-0.5
       logg_ext(:)%ebmv_val(iext)= (real(iext)/2.0)-0.5
    end do

    rewind(22)
    do ihead=1, nhead
      read(22,*)
    end do
    allocate(bc_in(n_filter))
    ! Now actually read the data.
    read_do: do 
      read(22,*,iostat=ios) temp_in, logg_in, ebmv_in, bc_in
      if (ios < 0) exit read_do
      ! Speed tests here show no difference in overall monte speed between minloc or
      ! or locate_nearest.
      !ilogg=minloc(abs(logg_in-grav_list),1)
      ilogg=locate_nearest(grav_list(1:n_logg), logg_in)
      !ientry=minloc(abs(all_ext(ilogg)%teff-temp_in),1)
      ientry=locate_nearest(all_ext(ilogg)%teff, temp_in)
      !iext=minloc(abs(all_ext(ilogg)%ebmv_val-ebmv_in),1)
      iext=locate_nearest(all_ext(ilogg)%ebmv_val, ebmv_in)
      !print*, ilogg, ientry, iext
      all_ext(ilogg)%exts(ientry,iext,1:n_filter)=bc_in
    end do read_do

    close(22)

    deallocate(bc_in)

  end subroutine read_extin_info

!!!!!!!!!!!!!

  subroutine closest_bound(value, n_array, array, close_int)
    ! This subroutine finds the closest bounds to value in array
    ! It assumes that the values increase with increasing array index.
    real :: value, range_low, range_high
    integer :: n_array
    real, dimension(n_array) :: array
    integer, dimension(2) :: close_int

    ! Find the bounding values
    !close_int(1)=minloc(abs(array-value), 1)
    close_int(1)=locate_nearest(array, value)
    range_low=min(array(1),array(n_array))
    range_high=max(array(1),array(n_array))

    ! Dependent upon the order of the reddening array i.e. this assumes
    ! reddest to bluest as we descend the array
    if (value >= range_low .and. value <= range_high) then
      if ((array(close_int(1))-value)*(array(n_array)-array(1)) > 0.) then
        close_int(2)=close_int(1)
        close_int(1)=close_int(1)-1
      else
        close_int(2)=close_int(1)+1
      end if
      ! Prepare for extrapolation
    else if (value < range_low) then
      close_int(2)=close_int(1)+1
    else if (value > range_high) then
      close_int(2)=close_int(1)-1
    end if

  end subroutine closest_bound

!!!!!!!!!!!!!

  subroutine treble_interp(logg_in, col_in, teff_in, nom_ebmv, n_logg, logg_ext, ext_out, red_out, extrap)
    ! This subroutine finds the required loggs, intrinsic colours and then interpolates for the
    ! required reddening (at the nominal E_B_V)

    ! Have we used extrapolation for the red/ext
    logical :: extrap

    ! Input
    real :: logg_in, col_in, teff_in, nom_ebmv
    integer :: n_logg
    type(spec_ext_values), dimension(n_logg) :: logg_ext

    ! The results
    real :: ext_out, red_out
    integer, dimension(2) :: close_logg
    integer, dimension(2,2) :: close_teff
    integer, dimension(2,2) :: close_ebv
    integer :: iclose

    ! The results of the interpolation
    real, dimension(2,2) :: ebmv_interp, ext_interp
    real, dimension(2) :: intr_interp_red, intr_interp_ext
    real :: logg_interp_red, logg_interp_ext

    ! Counters
    integer :: ilogg, iteff
    character(len=2) :: flag

    ! Assume that we have not extrapolated
    extrap = .false.

    ! We need the bounding logg files
    call closest_bound(logg_in, n_logg, logg_ext(:)%logg, close_logg)

    ! Find the required bounding effective temperatures and finally we find the bounding nominal E(B-V)
    do ilogg=1, 2
       call closest_bound(teff_in, logg_ext(close_logg(ilogg))%n_entries, &
            logg_ext(close_logg(ilogg))%temp, close_teff(:,ilogg))
       call closest_bound(nom_ebmv, n_ext, &
            logg_ext(close_logg(ilogg))%ebmv_val, close_ebv(:,ilogg))
    end do

    ! Now we have the locations in logg_ext of the bounding loggs in close_logg(1 and 2)
    ! We also have the bounding lines of col_in in each close_logg(1,2) logg_ext structure
    ! saved in close_int(1-2 and 1-2).
    ! i.e close_int(1,1) is the location in logg_ext(close_logg(1)) of the lower intrinsic col
    ! The same is true for close_ebv.

    ! Interpolate for the nominal E(B-V)
    ! Add the extrap logical and depending upon which flag is returned, we can say
    ! F (interpolated) or T (extrapolated).
    do ilogg=1, 2
       do iteff=1, 2
          !         print*, ilogg, iteff, '>'
          !         print*, close_logg(ilogg), close_teff(iteff,ilogg)
          ebmv_interp(iteff,ilogg)=linint(logg_ext(close_logg(ilogg))%ebmv_val(close_ebv(:,ilogg)),&
               logg_ext(close_logg(ilogg))%red(close_teff(iteff,ilogg),close_ebv(:,ilogg)), nom_ebmv, flag)
          if (trim(flag) == '1') extrap = .true.
          ext_interp(iteff,ilogg)=linint(logg_ext(close_logg(ilogg))%ebmv_val(close_ebv(:,ilogg)),&
               logg_ext(close_logg(ilogg))%ext(close_teff(iteff,ilogg),close_ebv(:,ilogg)), nom_ebmv, flag)
          if (trim(flag) == '1') extrap = .true.
       end do
    end do

    ! Now interpolate for the efective temperature
    do ilogg=1, 2
       intr_interp_red(ilogg)=linint(logg_ext(close_logg(ilogg))%temp(close_teff(:,ilogg)), &
            ebmv_interp(:,ilogg), teff_in, flag)
       if (trim(flag) == '1') extrap = .true.
       intr_interp_ext(ilogg)=linint(logg_ext(close_logg(ilogg))%temp(close_teff(:,ilogg)), &
            ext_interp(:,ilogg), teff_in, flag)
       if (trim(flag) == '1') extrap = .true.
    end do

    ! Finally, we interpolate for the logg
    logg_interp_red=linint(logg_ext(close_logg(:))%logg, intr_interp_red(:), logg_in, flag)
    if (trim(flag) == '1') extrap = .true.
    logg_interp_ext=linint(logg_ext(close_logg(:))%logg, intr_interp_ext(:), logg_in, flag)
    if (trim(flag) == '1') extrap = .true.

    ext_out=logg_interp_ext
    red_out=logg_interp_red

    !    print*, red_out, ext_out, logg_in, col_in, extrap

  end subroutine treble_interp

end module red_ext_interp
