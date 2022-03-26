      module colteff_subs

      use cmdfit_system
      use quad
      use bolcor_subs

      implicit none

      interface reddening
        module procedure reddening_point
        module procedure reddening_array
      end interface

      contains

      subroutine choose_models(iso_file, bc_file, ext_file, colnam)

        ! This routine reads the setup files with the names (and other details) of
        ! the interior and atmosphere models, and returns the names of the files
        ! to use.
        character(len=*), intent(out) :: iso_file, bc_file, ext_file
        character(len=*), dimension(:), intent(out) :: colnam

        character(len=512), dimension(:), allocatable :: description, files, aux, ext
        character(len=50) :: junk_str
        integer :: n_names, i_name, iostat, i_bc, n_bc, i_column, n_column, iminus
        character(len=10), allocatable, dimension(:) :: s_column, names

        open(unit=23, file=trim(data_dir())//'setup.int')
        n_names=0
        do 
          read(23,*,iostat=iostat) junk_str
          if (iostat < 0) exit
          n_names=n_names+1
        end do

        allocate(description(n_names), files(n_names))

        rewind(23)
        do i_name=1, n_names
          read(23,*) description(i_name), files(i_name)
        end do
        close(23)

        print*, 'The following interior models are available.'
        write(*,30) -1, 'Single age user isochrone file'
        do i_name=1, n_names
          if (len_trim(description(i_name)) > 0) write(*,30) i_name, &
          trim(description(i_name))
        end do
        do
          print*, 'What model number do you want?'
          read(*,*,iostat=iostat) i_name
          if (iostat == 0) exit
          print*, 'Error in input.'
        end do

        if (i_name == -1) then
          print*, '> Give the file name.'
          read(*,*) iso_file
          iso_file='user::'//iso_file
        else
          iso_file=trim(data_dir())//files(i_name)
        end if

        ! Does this model have any colours supplied?
        if (i_name == -1) then
          call read_head(iso_file(7:len(iso_file)), 2, n_column, s_column)
        else
          call read_head(iso_file, 2, n_column, s_column)
        end if
        allocate(names(iso_mcol))
        close(2)
        ! How many of these columns are absolute magnitudes or colours?
        n_bc=0
        do i_column=1, n_column
          iminus=index(s_column(i_column), '-')
          if (s_column(i_column)(1:2)=='M_') then
            n_bc=n_bc+1
            ! Simply remove the M at the beginning.
            names(n_bc)=s_column(i_column)(3:len(s_column(i_column)))
          else if (iminus > 0) then
            n_bc=n_bc+1
            names(n_bc)=s_column(i_column)
          end if
        end do

        if (n_bc > 0) then
          print*, 'This isochrone comes with the following colours ', names(1:n_bc)
        end if

        deallocate(description, files)

        open(unit=23, file=trim(data_dir())//'setup.bc')
        n_names=0
        do 
          read(23,*,iostat=iostat) junk_str
          if (iostat < 0) exit
          n_names=n_names+1
        end do

        allocate(description(n_names), files(n_names), aux(n_names), ext(n_names))

        rewind(23)
        do i_name=1, n_names
          read(23,*) description(i_name), files(i_name), ext(i_name), aux(i_name)
        end do
        close(23)

        print*, 'The following atmospheric models are available.'
        write(*,30) -1, 'user supplied file'
        do i_name=1, n_names
          if (len_trim(aux(i_name))==0 .or. & 
            trim(data_dir())//trim(aux(i_name))==trim(iso_file)) then
            if (len_trim(description(i_name)) > 0) &
            write(*,30) i_name, trim(description(i_name))
          end if
        end do
        if (n_bc > 0) then
          write(*,30) 0, 'Or the colours supplied with the isochrone'
        end if
        do
          print*, 'What model number do you want?'
          read(*,*,iostat=iostat) i_name
          if (iostat == 0) exit
          print*, 'Error in input.'
        end do
        if (i_name == -1) then
          print*, '> Give the name of the bolometric correction file.'
          read(*,*) bc_file
          print*, '> Give the name of the extinction file (<cr>=none).'
          ext_file=' '
          read(*,'(a)') ext_file
        else if (i_name > 0) then
          bc_file=trim(data_dir())//files(i_name)
          ext_file=trim(data_dir())//ext(i_name)
        else
          ! Must have chosen to use the colours supplied with the isochrone.
          bc_file=' '
          ext_file=trim(data_dir())//'with_isochrones.rv'
        end if

        if (i_name /= 0) then
          ! Now let's work out what colours are available and required.
          call read_head(bc_file, 2, n_column, s_column)
          close(2)
          ! How many of these columns are bolometric corrections?
          n_bc=0
          do i_column=1, n_column
            if (s_column(i_column)(1:3) == 'BC_') n_bc=n_bc+1
          end do
          if (allocated(names)) deallocate(names)
          allocate(names(n_bc))
          i_bc=0
          do i_column=1, n_column
            if (s_column(i_column)(1:3) == 'BC_') then
              iminus=index(s_column(i_column), '-')
              i_bc=i_bc+1
              if (iminus == 0) then
                ! Simply remove the BC at the beginning.
                names(i_bc)=s_column(i_column)(4:len(s_column(i_column)))
              else
                ! Recall that colours are negative of BCs.
                names(i_bc)=s_column(i_column)& 
                (iminus+1:len_trim(s_column(i_column)))//'-'//s_column(i_column)(4:iminus-1)
              end if
            end if
          end do

        end if

        deallocate(description, files, aux, ext)

        print*, 'Available filters are ', names(1:n_bc)
        mag: do
          print*, '> Which magnitude do you require?'
          read(*,*) colnam(1)
          if (find_bc_combination(colnam(1), names, n_bc)) exit mag
          print*, 'Colour not available.'
        end do mag
        col: do
          print*, '> Which colour do you require?'
          read(*,*) colnam(2)
          if (find_bc_combination(colnam(2), names, n_bc)) exit col
          print*, 'Colour not available.'
        end do col


        deallocate(names)
        

30      format(1x, i2, 2x, a)
      

      end subroutine choose_models

      subroutine allocate_chrone(iso)

        type(chrone) :: iso

        if (allocated(iso%mass)) then
          deallocate (iso%pnt, iso%mass, iso%teff, iso%logg, &
          iso%lbol, iso%massflg, iso%teffflg, &
          iso%loggflg, iso%lbolflg)
        end if

        allocate (iso%pnt(iso%npts), iso%mass(iso%npts), &
        iso%teff(iso%npts), iso%logg(iso%npts), iso%lbol(iso%npts), &
        iso%massflg(iso%npts), &
        iso%teffflg(iso%npts), iso%loggflg(iso%npts), iso%lbolflg(iso%npts))

        iso%mass=0.0
        iso%teff=0.0
        iso%logg=0.0
        iso%lbol=0.0

        iso%massflg='?'
        iso%teffflg='?'
        iso%loggflg='?'
        iso%lbolflg='?'

      end subroutine allocate_chrone

      subroutine iso_calc(iso)
      
        ! Given an isochrone model number and colour flag this routine
        ! interpolates the isochrone to a given age, and then pastes
        ! a model atmosphere on top.

        type(chrone), intent(inout) :: iso
        type(chrone), dimension(1) :: iso1
        
        integer :: imod, i, j, k, iflag, icol, ilen
        ! luminosity/bolometric corr. at a point               
        integer :: bcflag
        character(len=30) :: ifname
        character(len=30) :: wrkflg
        character(len=10), dimension(2) :: colnam

        real :: m_bol_sun

        ! Read the appropriate interior model.
        ilen=len_trim(iso%iso_file)
        if (iso%iso_file(1:6) == 'user::') then

          ! The general isochrone reader.  As we are pasting atmospheres on we
          ! call it with colour zero, so it does not attempt to read colours.
          colnam(1)=''
          call read_padova_col(iso%iso_file(7:len_trim(iso%iso_file)), iso%age, colnam, iso1)

          iso%ncol=iso1(1)%ncol
          if (iso%ncol > 0) iso%colnam(1:iso%ncol)=iso1(1)%colnam(1:iso%ncol)

          iso%npts=iso1(1)%npts
          call allocate_chrone(iso)
          iso%lbol=iso1(1)%lbol
          iso%teff=iso1(1)%teff
          iso%mass=iso1(1)%mass
          iso%logg=iso1(1)%logg
          iso%lbolflg=iso1(1)%lbolflg
          iso%teffflg=iso1(1)%teffflg
          iso%massflg=iso1(1)%massflg
          iso%loggflg=iso1(1)%loggflg

          iso%pnt=iso1(1)%pnt

          iso%colreq = iso1(1)%colreq

        else if (iso%iso_file(ilen-2:ilen)=='iso') then
          ! Its an isochrone.
          call read_padgen(iso)
        else if (iso%iso_file(ilen-2:ilen)=='trk') then
          ! Its (new format) tracks.
          call read_pms(iso%iso_file, .true., iso)
        else 
          ! It (old format) tracks.
          call read_pms(iso%iso_file, .false., iso)
        end if

        ! Really this should be done in the read routines (as read_pms does).
        iso%massflg='OK'
        iso%teffflg='OK'
        iso%loggflg='OK'
        iso%lbolflg='OK'

        ! lbol and teff derived for that isochrone, now convert to colour and mag
        ! Only read bolometric corrections if the user does not want any of the colours
        ! supplied with the isochrone.

        if (len_trim(iso%bc_file) > 0) then

          ! Get all the bolometric corrections, and apply them.
          call bolcor(iso, m_bol_sun)

          ! Tidy up any flags.
          create: do imod=1, iso%npts
            do icol=1, iso%ncol
              ! Is this a colour or a magnitude?
              if (index(iso%colnam(icol), '-') > 0) then
                !print*, iso%colnam(icol), ' is a colour'
                if (iso%pnt(imod)%col(icol)%flag == '1') &
                iso%pnt(imod)%col(icol)%flag='Colour-Teff'
              else
                !print*, iso%colnam(icol), ' is a magnitude'
                if (iso%pnt(imod)%col(icol)%flag == '1') &
                iso%pnt(imod)%col(icol)%flag='Bolometric correction'
                iso%pnt(imod)%col(icol)%data=m_bol_sun-iso%pnt(imod)%col(icol)%data-2.5*iso%lbol(imod)
              end if
            end do
          end do create

          !print*, 'imod, mag, iso%magflg, col, iso%colflg, Teff,'//&
          !' lbol, mass logg, bc'
          !do imod=1, iso%npts
          !  print*, imod, iso%pnt(imod)%col(1)%data, trim(iso%pnt(imod)%col(1)%flag), &
          !                iso%pnt(imod)%col(2)%data, trim(iso%pnt(imod)%col(1)%flag), 10.0**iso%teff(imod),&
          !  iso%lbol(imod), iso%mass(imod), iso%logg(imod)
          !  write(51,*) iso%col(imod), iso%mag(imod), iso%mass(imod)
          !end do

        end if
            

      end subroutine iso_calc


! *****************************************************************************
      subroutine reddening_array(colnam, ext_file, dmod, red, col, mag)
! *****************************************************************************
      
        ! Changes an absolute magnitude sequence into a reddened one at
        ! at given distance reddening (dmod, col positve) or the reverse.
    
        ! Note that red is the reddening in the colour (e.g. B-V).
    
        character(len=*), dimension(2), intent(in) :: colnam
        character(len=*) :: ext_file
        real, intent(in) :: dmod, red
        real, dimension(:), intent(inout) :: col, mag

        integer :: icol, nhed, ihed, iostat, ired, jred
        character(len=50) :: string
        integer, parameter :: mred=10
        character(len=10), dimension(mred), save :: colnam1, colnam2
        real, dimension(mred), save :: const, red1, red2, col1, col2
        integer, save :: nred, ipos
        character(len=80), save :: last_file=' '

        if (abs(red) > 2.0*tiny(red)) then
          if (trim(ext_file) /= trim(last_file)) then
            ipos=index(ext_file, '.', .true.)
            if (ext_file(ipos+1:len_trim(ext_file)) /= 'rv') then
              print*, 'Can only cope with an extinction vector (not a table) at this point.'
              stop
            end if
            print*, 'Reading extinction vectors from ', trim(ext_file)
            open(37, file=ext_file, action='read', status='old') 
            nhed=0
            do
              read(37,*,iostat=iostat) string
              if (iostat < 0) exit
              string=adjustl(string)
              if (string(1:1) /= '#') exit
              nhed=nhed+1
            end do

            rewind(37)
            do ihed=1, nhed
              read(37,*)
            end do

            nred=0
            do
              read(37,*,iostat=iostat) colnam1(nred+1), colnam2(nred+1), const(nred+1), &
              red1(nred+1), red2(nred+1), col1(nred+1), col2(nred+1)
              if (iostat < 0) exit
              nred=nred+1
              if (nred > mred-1) then
                print*, 'Array sizes need enlarging in reddening.'
                stop
              end if
            end do
            close(37)
            last_file=ext_file
          end if

          do ired=1, nred
            if (trim(colnam(1))==trim(colnam1(ired)) .and. trim(colnam(2))==trim(colnam2(ired))) then
              jred=ired
              exit
            end if
            if (ired == nred) then
              print*, 'Cannot colours find the colours', trim(colnam(1)), &
              ' ', trim(colnam(2)), ' in file.'
              stop
            end if
          end do
          
          ! Change to un-reddened colour (order important for Bessell).
          if (red < 0.0) col=col+red

          mag = mag + red*(const(ired) + (red1(ired)*abs(red)) + (red2(ired)*(red*red)) + (col1(ired)*col) &
          + (col2(ired)*(col*col)))

          ! Change to un-reddened colour (order important for Bessell).
          if (red > 0.0) col=col+red

        end if

        ! And finally deal with the distance modulus.
        mag=mag+dmod
      

      end subroutine reddening_array

! *****************************************************************************
      subroutine reddening_point(colnam, ext_file, dmod, red, col, mag)
! *****************************************************************************

        character(len=*), dimension(2), intent(in) :: colnam
        character(len=*) :: ext_file
        real, intent(in) :: dmod, red
        real, intent(inout) :: col, mag
        
        real, dimension(1) :: a_col, a_mag
        
        a_col(1)=col; a_mag(1)=mag
        
        call reddening_array(colnam, ext_file, dmod, red, a_col, a_mag)

        col=a_col(1); mag=a_mag(1)
      
      end subroutine reddening_point  


      subroutine read_pms(ifname, trk, iso)

      ! Modelled on Rob's Baraffe routine, this is a more generalisable
      ! routine for reading and interpolating PMS isochrones.

      ! At this stage not very general, designed to read the
      ! Palla and Stahler files Cameron created, but should work for
      ! Baraffe.

      ! Slowly drifting towards reading files in the new format.  At the
      ! moment these are defined as being .trk files.

      character(len=*), intent(in) :: ifname
      logical, intent(in) :: trk
      type(chrone) :: iso

      ! Locals
      integer :: n_masses, i_masses, max_n_ages, i_ages
      integer, allocatable, dimension(:) :: n_ages
      
      real, allocatable, dimension(:,:) :: logAge, logL, logT, log_G
      real, allocatable, dimension(:) :: dmass

      real :: old_mass, new_mass, dummy
      integer :: iostat, nmod, ihead, nhead, i
      character(len=50) :: flag, instring


      ! Tracks are all in one big file
      ! idea is to read in a file for each mass and interpolate to get
      ! log L and log T at the appropriate log(age)

      open (1,file=ifname,&
      form='formatted',access='sequential',status='old')

      print*, 'Reading interior models from file ', trim(ifname)

      if (trk) then
        nhead=count_header(1)
      else
        nhead=3
      end if

      ! Start by counting the number of masses.
      rewind(1)
      do ihead=1, nhead
        read(1,*) 
      end do
      n_masses=0
      old_mass=-1
      do
        read(1,*,iostat=iostat) new_mass
        if (iostat < 0) exit
        if (new_mass > nearest(old_mass,1.0)) then
          ! If the new mass is bigger than the old one.
          n_masses=n_masses+1
          old_mass=new_mass
        end if
      end do

      ! Now count the number of ages at each mass.
      allocate(n_ages(n_masses))
      rewind(1)
      do ihead=1, nhead
        read(1,*) 
      end do
      old_mass=-1.0
      i_masses=0
      each_mass: do
        read(1,*,iostat=iostat) new_mass
        if (iostat < 0) exit each_mass
        if (new_mass > nearest(old_mass,1.0)) then
          ! If the new mass is bigger than the old one.
          old_mass=new_mass
          i_masses=i_masses+1
          n_ages(i_masses)=1
        else
          n_ages(i_masses)=n_ages(i_masses)+1
        end if
      end do each_mass
      max_n_ages=maxval(n_ages,1)

      ! Now read the data in.
      rewind(1)
      do ihead=1, nhead
        read(1,*) 
      end do
      allocate(dmass(n_masses), logAge(n_masses, max_n_ages))
      allocate(logT(n_masses, max_n_ages), logL(n_masses, max_n_ages))
      allocate(log_G(n_masses, max_n_ages))
      do i_masses=1, n_masses

        do i_ages=1, n_ages(i_masses)
          read(1,*) dmass(i_masses),logAge(i_masses,i_ages), &
          logT(i_masses,i_ages), log_G(i_masses,i_ages), &
          logL(i_masses,i_ages)
        end do

        ! Convert to logs
        if (.not. trk) then
          do i_ages=1,n_ages(i_masses)
            logAge(i_masses,i_ages)=log10(logAge(i_masses,i_ages))+9.0
            logT(i_masses,i_ages)=log10(logT(i_masses,i_ages))
          end do
        end if

      end do

      close(1)

      ! Now have the models in the 2d grids
      ! for each row in the grid, need to find the index closest to the
      ! required age.
      ! Unfortunately we have to interpolate twice.  Once to count the
      ! number of elements in the model
      iso%npts=0
      count_masses: do i_masses=1, n_masses
        dummy=quadintflg(logAge(i_masses,1:n_ages(i_masses)), &
        logL(i_masses,1:n_ages(i_masses)), iso%age, flag)
        if (trim(flag) /= 'OK') cycle count_masses
        dummy=quadintflg(logAge(i_masses,1:n_ages(i_masses)), &
        logT(i_masses,1:n_ages(i_masses)), iso%age, flag)
        if (trim(flag) /= 'OK') cycle count_masses
        dummy=quadintflg(logAge(i_masses,1:n_ages(i_masses)), &
        log_G(i_masses,1:n_ages(i_masses)), iso%age, flag)
        if (trim(flag) /= 'OK') cycle count_masses
        iso%npts=iso%npts+1
      end do count_masses


      if (iso%npts == 0) then
        print*, 'Could not find appropriate age.'
        print*, 'Required log10(age in years) is ', iso%age
        stop    
      end if

      call allocate_chrone(iso)

      nmod=0
      do i_masses=1, n_masses
        iso%lbol(nmod+1)=quadintflg(logAge(i_masses,1:n_ages(i_masses)), &
        logL(i_masses,1:n_ages(i_masses)), iso%age, iso%massflg(nmod+1))
        if (trim(iso%massflg(nmod+1)) /= 'OK') cycle
        iso%teff(nmod+1)=quadintflg(logAge(i_masses,1:n_ages(i_masses)), &
        logT(i_masses,1:n_ages(i_masses)), iso%age, iso%teffflg(nmod+1))
        if (trim(iso%teffflg(nmod+1)) /= 'OK') cycle
        iso%logg(nmod+1)=quadintflg(logAge(i_masses,1:n_ages(i_masses)), &
        log_G(i_masses,1:n_ages(i_masses)), iso%age, iso%loggflg(nmod+1))
        if (trim(iso%loggflg(nmod+1)) /= 'OK') cycle
        iso%mass(nmod+1)=dmass(i_masses)
        nmod=nmod+1
        ! If the last mass is uninterpolatable we'd better escape or we'll
        ! overrun the arrays.
        if (nmod == iso%npts) exit
        ! print*, nmod, mass(nmod), lbol(nmod), teff(nmod), logg(nmod)
      end do

      end subroutine read_pms


! ******************************************************************************
      subroutine read_padgen(iso)
! ******************************************************************************

      implicit none

      type(chrone) :: iso

      integer :: nmod, icol
      type(chrone), dimension(3) :: isoc
      integer :: imod, imass, jpt, iage
      real, dimension(3) :: teff_in, lbol_in, logg_in, age_in
      real, dimension(:,:), allocatable :: col_in
      real :: uncer, min_mass, max_mass
      character(len=10), dimension(2) :: colnam

      character(len=30) :: flag

      colnam=''
      call read_padova_col(iso%iso_file, iso%age, colnam, isoc)

      ! What mass scale to use?
      ! At some point the line below got changed to min_mass=maxval(isoc(1)%mass). The older versions
      ! had isoc%mass(1), which is conceptually right, but may not compile.  Changed by Timn Jan 2012.
      min_mass=max(isoc(1)%mass(1), isoc(2)%mass(1))
      min_mass=max(min_mass, isoc(3)%mass(1))
      max_mass=min(isoc(1)%mass(isoc(1)%npts), isoc(2)%mass(isoc(2)%npts))
      max_mass=min(max_mass, isoc(3)%mass(isoc(3)%npts))

      iso%npts=0
      count_mass: do imass=1, isoc(2)%npts
        if (isoc(2)%mass(imass) < min_mass) cycle count_mass
        if (isoc(2)%mass(imass) > max_mass) cycle count_mass
        iso%npts=iso%npts+1
      end do count_mass

      call allocate_chrone(iso)

      iso%ncol=isoc(1)%ncol
      iso%colnam=isoc(1)%colnam
      if (iso%ncol > 0) allocate(col_in(3,iso%ncol))

      nmod=0
      ! At each mass point in the middle isochrone.
      each_mass: do imass=1, isoc(2)%npts
        if (isoc(2)%mass(imass) < min_mass) cycle each_mass
        if (isoc(2)%mass(imass) > max_mass) cycle each_mass
        nmod=nmod+1
        iso%mass(nmod)=isoc(2)%mass(imass)
        
        ! Set up the middle of the three values we will interpolate accross.
        teff_in(2)=isoc(2)%teff(imass)
        lbol_in(2)=isoc(2)%lbol(imass)
        logg_in(2)=isoc(2)%logg(imass)
        age_in(2)=isoc(2)%age
        if (iso%ncol > 0) &
        col_in(2,1:iso%ncol)=isoc(2)%pnt(imass)%col(1:iso%ncol)%data

        ! Now interpolate the other two points.
        each_age: do iage=1, 3, 2

          age_in(iage)=isoc(iage)%age

          teff_in(iage)=linint(isoc(iage)%mass(1:isoc(iage)%npts), &
          isoc(iage)%teff(1:isoc(iage)%npts), iso%mass(nmod), flag)
          lbol_in(iage)=linint(isoc(iage)%mass(1:isoc(iage)%npts), &
          isoc(iage)%lbol(1:isoc(iage)%npts), iso%mass(nmod), flag)
          logg_in(iage)=linint(isoc(iage)%mass(1:isoc(iage)%npts), &
          isoc(iage)%logg(1:isoc(iage)%npts), iso%mass(nmod), flag)

          do icol=1, iso%ncol
            col_in(iage,icol)=linint(isoc(iage)%mass(1:isoc(iage)%npts), &
            isoc(iage)%pnt(1:isoc(iage)%npts)%col(icol)%data, iso%mass(nmod), flag)
          end do

        end do each_age

        ! Now interpolate over the three points to get the right age.
        iso%teff(nmod)=linint(age_in, teff_in, iso%age, flag)
        iso%lbol(nmod)=linint(age_in, lbol_in, iso%age, flag)
        iso%logg(nmod)=linint(age_in, logg_in, iso%age, flag)
        do icol=1, iso%ncol
          iso%pnt(nmod)%col(icol)%data=linint(age_in, col_in(:,icol), iso%age, flag)
          iso%pnt(nmod)%col(icol)%flag='OK'
        end do

      end do each_mass

      end subroutine read_padgen

      subroutine read_padova_col(ifname, age, colnam, isoc)

      ! A very general isochrone reader.  If isoc is of dimension one
      ! the reader expects to find an isochrone at exactly the right 
      ! age.  If it is if dimension three it will return the three
      ! nearest in age.

      ! Define the isochrone by the file to be read from (ifname), the age and
      ! strings for the magnitude and colour (colnam).

      character(len=*) :: ifname
      real, intent(in) :: age
      character(len=*), dimension(2) :: colnam
      type(chrone), dimension(:), intent(out) :: isoc
      type(cols) :: col_num

      integer, parameter :: maxcol=200
      integer :: n_ages, iostat, i_ages, i, jage, j
      real :: old_mass, new_mass, old_age, new_age, junk
      integer, allocatable, dimension(:) :: n_masses
      real, allocatable, dimension(:) :: r_masses
      real, allocatable, dimension(:) :: log10_ages
      integer, allocatable, dimension(:) :: icnt
      real :: teff, logg, mag_bol, lbol
      real :: mag_u, mag_b, mag_v, mag_r, mag_i, mag_j, mag_h, mag_k
      integer :: icol_plus, icol_minus, imag_plus, imag_minus

      character(len=10) :: col_plus, col_minus, mag_plus, mag_minus
      character(len=20), dimension(:), allocatable :: col_head
      real, dimension(maxcol) :: data
      integer :: ihed, nhed, ipos
      character(len=20) work_str
      integer, dimension(iso_mcol) :: jcol
      integer :: icol, iminus

      logical :: first=.true.
      logical, parameter :: debug=.false.

      if (first) then
        print*, 'Using interior model data file ', trim(ifname)
        first=.false.
      end if

      call read_head(ifname, 2, nhed, col_head, col_num)
      if (col_num%ilbol == 0) print*, 'Could not find bolometric luminosity.'
      if (col_num%iage  == 0) print*, 'Could not find age.'
      if (col_num%iteff==0 .and. col_num%ilogte==0) print*, 'Could not find effective temperature.'
      if (col_num%imass == 0) print*, 'Could not find mass.'
      if (col_num%ilogg == 0) print*, 'Could not find logG.'

      ! Now get the names of the column heads.
      each_header: do ihed=1, nhed
        ! Get rid of ' and *.
        ipos=index(col_head(ihed), char(39))
        if (ipos > 0) col_head(ihed)(ipos:ipos)=' '
        ipos=index(col_head(ihed), ' ')
        if (ipos > 0) col_head(ihed)(ipos:ipos)=' '
      end do each_header

      ! First get a feel for the ages available.
      if (debug) print*, 'The column with the ages in is number ', col_num%iage
      n_ages=0
      rewind(2)
      each_line1: do
        iostat=read_line(2, nhed, data)
        if (iostat < 0) then
          exit each_line1
        else if (iostat > 0) then
          print*, 'Error ', iostat, ' reading file in each_line2.'
          stop
        end if
        new_age=data(col_num%iage)
        if (n_ages == 0) then
          n_ages=1
          old_age=new_age
        else if (abs(new_age-old_age) > age_diff_min/2.0) then
          n_ages=n_ages+1
          old_age=new_age
        end if
      end do each_line1
      if (debug) print*, 'Number of ages available is ', n_ages

      ! Now count the number of mass points at each age.
      rewind(2)
      allocate(n_masses(n_ages), log10_ages(n_ages))
      n_masses=0
      i_ages=0
      each_line2: do
        iostat=read_line(2, nhed, data)
        if (iostat < 0) then
          exit each_line2
        else if (iostat > 0) then
          print*, 'Error ', iostat, ' reading file in each_line2.'
          stop
        end if
        new_age=data(col_num%iage)
        new_mass=data(col_num%imass)
        if (i_ages == 0) then
          i_ages=1
          old_age =new_age
          old_mass=new_mass
          n_masses(i_ages)=n_masses(i_ages)+1
          log10_ages(i_ages)=new_age
        else if (abs(new_age-old_age) > age_diff_min/2.0) then
          i_ages=i_ages+1
          old_age =new_age
          old_mass=new_mass
          n_masses(i_ages)=1
          log10_ages(i_ages)=new_age
        else
          if (abs(new_mass-old_mass) > 2.0*(epsilon(new_mass)+epsilon(old_mass))) then
            old_age =new_age
            old_mass=new_mass
            n_masses(i_ages)=n_masses(i_ages)+1
          end if
        end if
      end do each_line2

      ! Now sort the age array and associated number of points by age.
      allocate(r_masses(n_ages))
      r_masses=real(n_masses)
      call sort(log10_ages, r_masses, n_ages)
      n_masses=nint(r_masses)
      deallocate(r_masses)

      ! Find the nearest age.
      jage=minloc(abs(log10_ages-age),1)
      !print*, 'Nearest age is ', log10_ages(jage), &
      !'Myr number ', jage
      if (size(isoc,1) == 1) then
        ! Do we have an exact match?
        if (abs(log10_ages(jage)-age) < 0.0001) then
          isoc(1)%age=age
          ! if (size(isoc) > 1) isoc(2)%age=0.0
          ! if (size(isoc) > 2) isoc(3)%age=0.0
          isoc(1)%npts=n_masses(jage)
          call allocate_chrone(isoc(1))
        else 
          print*, 'read_padova_col called for one isochrone when there '
          print*, 'is not an exact match in the file.'
          stop
        end if
      else
        ! Are we outside the available age range?
        if (log10_ages(1) > age) then
          print*, 'Youngest isochrone available is ', log10_ages(1)
          print*, 'So cannot do an age of ', age
          stop
        else if (log10_ages(n_ages) < age) then
          print*, 'Oldest isochrone available is ', log10_ages(n_ages)
          print*, 'So cannot do an age of ', age
          stop
        end if
        ! If the nearest age is the first one, but its above that age
        ! then take the second age as the mid-point of the interpolation.
        if (jage == 1) jage=2
        ! Likewise, if the nearest age is the oldest in the sequence, move the 
        ! interpolation mid-point.
        if (jage == n_ages) jage=n_ages-1
        isoc%age=log10_ages(jage-1:jage+1)
        isoc%npts=n_masses(jage-1:jage+1)
        do i_ages=1, 3
          if (debug) print*, 'Using an age of ', isoc(i_ages)%age
          call allocate_chrone(isoc(i_ages))
        end do
      end if

      ! Any absolute magnitudes or colours?  Note we are filling these things as
      ! array oprations as isoc could be of dimension 1 or 3.
      isoc%ncol=0
      do ihed=1, nhed
        iminus=index(col_head(ihed), '-')
        if (col_head(ihed)(1:2)=='M_') then
          isoc%ncol=isoc%ncol+1
          ! Simply remove the M at the beginning.
          do j=1, size(isoc,1)
            isoc(j)%colnam(isoc%ncol)=col_head(ihed)(3:len(col_head(ihed)))
          end do
          jcol(isoc(1)%ncol)=ihed
        else if (iminus > 0) then
          isoc%ncol=isoc%ncol+1
          do j=1, size(isoc,1)
            isoc(j)%colnam(isoc%ncol)=col_head(ihed)
          end do
          jcol(isoc(1)%ncol)=ihed
        end if
      end do
      
      ! Now get the data for the one or three ages we want, which are in isoc(1:3)%age.
      rewind(2)
      allocate(icnt(size(isoc)))
      icnt=0
      old_mass=0.0
      each_line3: do
        iostat=read_line(2, nhed, data)
        if (iostat < 0) then
          exit each_line3
        else if (iostat > 0) then
          print*, 'Error ', iostat, ' reading file in each_line3.'
          stop
        end if
        new_age=data(col_num%iage) 
        new_mass=data(col_num%imass)
        !j=minval(minloc(abs(new_age-(log10(isoc%age)+6.0))))
        j=minloc(abs(new_age-isoc%age),1)
        if (abs(new_age-isoc(j)%age) < age_diff_min/2.0) then
          if (abs(old_mass-new_mass) > 2.0*(epsilon(new_mass)+epsilon(old_mass))) then 
            icnt(j)=icnt(j)+1
            if (icnt(j) > isoc(j)%npts) then
              print*, 'There appear to be more points for the isochrone of age ', isoc(j)%age
              print*, 'than I first thought.  Is this file propoerly sorted in age then mass?'
              stop
            end if
            isoc(j)%lbol(icnt(j))=data(col_num%ilbol)
            if (col_num%ilogte > 0) then
              isoc(j)%teff(icnt(j))=data(col_num%ilogte)
            else
              isoc(j)%teff(icnt(j))=log10(data(col_num%iteff))
            end if
            isoc(j)%mass(icnt(j))=data(col_num%imass)
            isoc(j)%logg(icnt(j))=data(col_num%ilogg)
            do icol=1, isoc(j)%ncol
              isoc(j)%pnt(icnt(j))%col(icol)%data=data(jcol(icol))
              isoc(j)%pnt(icnt(j))%col(icol)%flag='OK'
            end do
            isoc(j)%lbolflg(icnt(j))='OK'
            isoc(j)%teffflg(icnt(j))='OK'
            isoc(j)%massflg(icnt(j))='OK'
            isoc(j)%loggflg(icnt(j))='OK'
!           print*, new_age, icnt(j), isoc(j)%mass(icnt(j)), &
!           isoc(j)%mag(icnt(j)), isoc(j)%col(icnt(j))
          end if
          old_mass=new_mass
        end if
      end do each_line3
      ! print*, 'icnt ', icnt
      deallocate(icnt)

      deallocate(n_masses, log10_ages)
      close(2)

      do j=1, size(isoc,1)
        if (len_trim(colnam(1)) /= 0) isoc(j)%colreq=colnam
        isoc(j)%iso_file=ifname
      end do

      end subroutine read_padova_col


    end module colteff_subs
