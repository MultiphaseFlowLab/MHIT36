subroutine read_and_stats(nstep)

use commondata

integer :: nstep
character(len=40) :: namedir, namefile
character(len=8) :: numfile
integer :: i,j,k,m


namedir='../../channel/output/'
write(numfile,'(i8.8)') nstep


allocate(u(nx,ny,nz))
allocate(v(nx,ny,nz))
allocate(w(nx,ny,nz)) ! add nz+1 to account for staggered grid? the top layer is missing (all zero); interpolate at cell center? to aling with u and v?
allocate(phi(nx,ny,nz))

write(*,*) 'Reading step ',nstep,' out of ',nend,' , flow'
  !reading u
  namefile=trim(namedir)//'u_'//numfile//'.dat'
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) u
  close(666,status='keep')
  !reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) v
  close(667,status='keep')
  !reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) w
  close(668,status='keep')


! interpolate w at the cell center (where u and v are also located)
! k=1 is the bottom node and k=nz is the upper most node (nz+1 is not saved to avoid differen file sizes, it's technically an halo point)
if (uflag .eq. 1) then
  do k=1,nz-1
    do i=1,nx
      do j=1,ny
        w(i,j,k)=0.5d0*(w(i,j,k) + w(i,j,k+1))
      enddo
    enddo
  enddo
  do i=1,nx
    do j=1,ny
      w(i,j,nz)=0.5d0*(w(i,j,nz) + 0.d0) ! assume nz+1 is 0 (no-slip)
    enddo
  enddo
endif


!--- compute statistics ---
allocate(mean(nz,3),rms(nz,3),skw(nz,3),flt(nz,3))

mean=0.0d0
rms=0.0d0
skw=0.0d0
flt=0.0d0

! mean
do k=1,nz
  do i=1,nx
    do j=1,ny
      mean(k,1) = mean(k,1) + u(i,j,k)
      mean(k,2) = mean(k,2) + v(i,j,k)
      mean(k,3) = mean(k,3) + w(i,j,k)
    enddo
  enddo
enddo
mean=mean/(dble(nx*ny))

! rms
do k=1,nz
  do i=1,nx
    do j=1,ny
      rms(k,1) = rms(k,1) + (u(i,j,k)-mean(k,1))**2
      rms(k,2) = rms(k,2) + (v(i,j,k)-mean(k,2))**2
      rms(k,3) = rms(k,3) + (w(i,j,k)-mean(k,3))**2
    enddo
  enddo
enddo
rms=rms/(dble(nx*ny))
rms=rms**0.5d0

! skw
do k=1,nz
  do i=1,nx
    do j=1,ny
      skw(k,1)=skw(k,1)+(u(i,j,k)-mean(k,1))**3
      skw(k,2)=skw(k,2)+(v(i,j,k)-mean(k,2))**3
      skw(k,3)=skw(k,3)+(w(i,j,k)-mean(k,3))**3
    enddo
  enddo
enddo
skw=skw/(dble(nx*ny))


! skw
do k=1,nz
  do i=1,nx
    do j=1,ny
      flt(k,1)=flt(k,1)+(u(i,j,k)-mean(k,1))**4
      flt(k,2)=flt(k,2)+(v(i,j,k)-mean(k,2))**4
      flt(k,3)=flt(k,3)+(w(i,j,k)-mean(k,3))**4
    enddo
  enddo
enddo
flt=flt/(dble(nx*ny))

! normalization for SKW and FLT
do k=1,nz
  do m=1,3
    skw(k,m)=skw(k,m)/rms(k,m)**3
    flt(k,m)=flt(k,m)/rms(k,m)**4
  enddo
enddo


namefile = 'output/stat_'//trim(numfile)//'.dat'
!write(*,*) "name", namefile
open(66,status='replace',file=trim(namefile),form='formatted')
write(66,'(13(a12,2x))') '% z','u mean','v mean','w mean','u rms','v rms','w rms','u skw','v skw','w skw','u flt','v flt','w flt'
!write(66,*)

do k=1,nz
  write(66,'(f12.5,2x,12(es12.5,2x))') z(k),mean(k,1),mean(k,2),mean(k,3),rms(k,1),rms(k,2),rms(k,3), &
                                            skw(k,1),skw(k,2),skw(k,3),   flt(k,1),flt(k,2),flt(k,3)
end do
close(66)

deallocate(u,v,w,phi)
deallocate(mean,rms,skw,flt)

return
end
