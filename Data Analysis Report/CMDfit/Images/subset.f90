        program subset

!       Outputs a subset of the input file.
!       Timn, January 1992.
!       Became almost trivial in Fortran 90, February 1998.

        use ark_file_io

        implicit none
 
        integer, dimension(2) :: naxis, ibeg, iend
        integer :: ndim
        real, allocatable, dimension(:,:) :: data, axdata
        character :: axch
        integer :: i, j, ifail, iostat, isafe
        real :: beg, end

        call setbug()
 
        print*, '* Enter "end" to exit program.'
        endless: do

          ifail=inpark(naxis, data, axdata)
          if (ifail == -1) exit endless
          if (naxis(2) == 1) then
            ndim=1
          else
            ndim=2
          end if

          iend(2)=1
          ibeg(2)=1
          prompt: do j=1, ndim
            axch='X'
            if (j .eq. 2) axch='Y'
            reader: do
              print*, '* The range of '//axch//' values is ', &
              axdata(1,j), ' to ', axdata(naxis(j),j)
              print*,'> Give the range of '//axch//' values required.'
              read (*,*,iostat=iostat) beg, end
              if (iostat == 0) exit reader
            end do reader
            ibeg(j)=minloc(abs(beg-axdata(1:naxis(j),j)),1)
            iend(j)=minloc(abs(end-axdata(1:naxis(j),j)),1)
            if (ibeg(j) > iend(j)) then
              isafe=iend(j)
              iend(j)=ibeg(j)
              ibeg(j)=isafe
            end if
          end do prompt

          print*, '* Will subset the ', &
          (iend(1)-ibeg(1)+1)*(iend(2)-ibeg(2)+1), ' points, with'
          explain: do j=1, ndim
            axch='X'
            if (j .eq. 2) axch='Y'
            print*, '* '//axch//' values ', axdata(ibeg(j),j), ' to ', &
            axdata(iend(j),j)
            print*, '* Pixels ', ibeg(j), ' to ', iend(j)
          end do explain

!         Correct the axis arrays.
          do i=1, 2
            do j=1, iend(i)-ibeg(i)+1
              axdata(j,i)=axdata(j+ibeg(i)-1,i)
            end do
          end do

!         A little array operation.
          naxis=iend-ibeg+1

          print*, 'Here'

          ifail=makark(naxis, data(ibeg(1):iend(1),ibeg(2):iend(2)), axdata)

          deallocate(data)
          deallocate(axdata)

        end do endless

        end program subset
