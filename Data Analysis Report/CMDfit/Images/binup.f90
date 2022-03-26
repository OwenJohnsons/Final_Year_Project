        program binup

        use ark_file_io

        implicit none

!       Fortran 90 version by Timn.

        real, allocatable :: data(:,:), axdata(:,:)
        integer :: naxis(2), ndim
 
        integer :: i, j, nred(2)
        logical :: first
        data first /.true./
                   
        call setbug()
 
90      i=inpark(naxis, data, axdata) 
        if (i .eq. -1) goto 900

        if (first) then
          if (naxis(1) .ne. 1) then
5           print*, '> Give AXIS1 reduction factor (1=no reduction).'
            read(*,*,err=5) nred(1)
          else
            nred(1)=1
          end if
          if (naxis(2) .ne. 1) then
100         print*, '> Give AXIS2 reduction factor (1=no reduction).'
            read(*,*,err=100) nred(2)
          else
            nred(2)=1
          end if
          first=.false.
        end if

        call doit(data, naxis, nred)
        naxis(1)=naxis(1)/nred(1)
        naxis(2)=naxis(2)/nred(2)

        ndim=1
        if (naxis(2) .gt. 1) ndim=2
        do i=1, ndim
          do j=1, naxis(i)*nred(i), nred(i)
            axdata(((j-1)/nred(i))+1,i)=(axdata(j,i)+axdata(j+nred(i)-1,i))/2.0
          end do
        end do

        i=makark(naxis, data, axdata)
        deallocate(data)
        deallocate(axdata)
        goto 90                                                          
 
900     end


        subroutine doit(data, naxis, nred)
 
        integer :: naxis(2), nred(2)
        real :: data(naxis(1), naxis(2))

        integer :: is, i, j, k
  
        if (nred(1) .ne. 1) then
          do 30 k=1, naxis(2)
            do 20 i=1, naxis(1)/nred(1)
              is=nred(1)*(i-1) + 1
              data(i,k)=data(is,k)
              do 10 j=is+1, is+(nred(1)-1)
                data(i,k)=data(i,k)+data(j,k)
10            continue
20          continue
30        continue
        end if
 
        if (nred(2) .ne. 1) then
          do 130 k=1, naxis(1)/nred(1)
            do 120 i=1, naxis(2)/nred(2)
              is=nred(2)*(i-1) + 1
              data(k,i)=data(k,is)
              do 110 j=is+1, is+(nred(2)-1)
                data(k,i)=data(k,i)+data(k,j)
110           continue
120         continue
130       continue
        end if
 
500     end
