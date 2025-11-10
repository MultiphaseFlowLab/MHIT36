subroutine read_and_stats(nstep)

use commondata

integer :: nstep
character(len=40) :: namedir, namefile
character(len=8) :: numfile
integer :: i,j,k,m
double precision :: tw


namedir='../../channel/output/'
write(numfile,'(i8.8)') nstep


allocate(theta(nx,ny,nz))

write(*,*) 'Reading step ',nstep,' out of ',nend,' , flow'
  !reading theta
  namefile=trim(namedir)//'theta_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) theta
  close(668,status='keep')


!--- compute statistics ---
allocate(mean(nz),rms(nz),skw(nz),flt(nz))

mean=0.0d0
rms=0.0d0
skw=0.0d0
flt=0.0d0

! remove wall temperature (make symmetric)
!do k=1,nz
!  do i=1,nx
!    do j=1,ny
!      if (k .le. nz/2) theta(i,j,k) = 1.d0-theta(i,j,k)
!      if (k .gt. nz/2) theta(i,j,k) = theta(i,j,k) + 1.d0
!    enddo
!  enddo
!enddo

! mean
do k=1,nz
  do i=1,nx
    do j=1,ny
      mean(k) = mean(k) + theta(i,j,k)
    enddo
  enddo
enddo
mean=mean/(dble(nx*ny))
!write(*,*) "mean", mean

! rms
do k=1,nz
  do i=1,nx
    do j=1,ny
      rms(k) = rms(k) + (theta(i,j,k)-mean(k))**2
    enddo
  enddo
enddo
rms=rms/(dble(nx*ny))
rms=rms**0.5d0

! skw
do k=1,nz
  do i=1,nx
    do j=1,ny
      skw(k)=skw(k)+(theta(i,j,k)-mean(k))**3
    enddo
  enddo
enddo
skw=skw/(dble(nx*ny))


! skw
do k=1,nz
  do i=1,nx
    do j=1,ny
      flt(k)=flt(k)+(theta(i,j,k)-mean(k))**4
    enddo
  enddo
enddo
flt=flt/(dble(nx*ny))

! normalization for SKW and FLT
do k=1,nz
  skw(k)=skw(k)/rms(k)**3
  flt(k)=flt(k)/rms(k)**4
enddo


namefile = 'output/tempstat_'//trim(numfile)//'.dat'
!write(*,*) "name", namefile
open(66,status='replace',file=trim(namefile),form='formatted')
write(66,'(5(a12,2x))') '% z',' mean','rms','skw','flt'
!write(66,*)

do k=1,nz
  write(66,'(f12.5,2x,4(es12.5,2x))') z(k),mean(k),rms(k),skw(k),flt(k)
end do
close(66)

deallocate(theta)
deallocate(mean,rms,skw,flt)

return
end
