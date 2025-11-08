program get_interface

use commondata
implicit none

integer :: i
nx=512
lx=2*pi
dx=lx/(nx-1)
dxi=1.d0/dx
nu=0.000625d0

write(*,*) "dxi", dxi
write(*,*) "nu", nu
write(*,*) "lx", lx

! velocity fields
allocate(u(nx,nx,nx))
allocate(v(nx,nx,nx))
allocate(w(nx,nx,nx))

! derivatives
allocate(dudx(nx,nx,nx))
allocate(dudy(nx,nx,nx))
allocate(dudz(nx,nx,nx))
allocate(dvdx(nx,nx,nx))
allocate(dvdy(nx,nx,nx))
allocate(dvdz(nx,nx,nx))
allocate(dwdx(nx,nx,nx))
allocate(dwdy(nx,nx,nx))
allocate(dwdz(nx,nx,nx))

nstart=0
ndump=1000
nend=35000

! #### MAIN LOOOP
do i=nstart,nend,ndump
  !write(*,*) 'Step ',i,' of ',nend
  call compute_diss(i)
enddo


deallocate(u,v,w)
deallocate(dudx,dudy,dudz)
deallocate(dvdx,dvdy,dvdz)
deallocate(dwdx,dwdy,dwdz)

return
end program








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_diss(nstep)

use commondata
integer :: nstep,i,j,k
double precision :: S11,S22,S33,S12,S13,S23

call read_fields(nstep)

!write(*,*) "max u", maxval(u)


do k=1,nx
  do j=1,nx
    do i=1,nx
       ip=i+1
       jp=j+1
       kp=k+1
       im=i-1
       jm=j-1
       km=k-1
       if (ip .gt. nx) ip=1
       if (jp .gt. nx) jp=1
       if (kp .gt. nx) kp=1   
       if (im .lt. 1) im=nx
       if (jm .lt. 1) jm=nx
       if (km .lt. 1) km=nx 
       dudx(i,j,k)= (u(ip,j,k) - u(im,j,k))*0.5d0*dxi
       dudy(i,j,k)= (u(i,jp,k) - u(i,jm,k))*0.5d0*dxi
       dudz(i,j,k)= (u(i,j,kp) - u(i,j,km))*0.5d0*dxi
       dvdx(i,j,k)= (v(ip,j,k) - v(im,j,k))*0.5d0*dxi
       dvdy(i,j,k)= (v(i,jp,k) - v(i,jm,k))*0.5d0*dxi
       dvdz(i,j,k)= (v(i,j,kp) - v(i,j,km))*0.5d0*dxi
       dwdx(i,j,k)= (w(ip,j,k) - w(im,j,k))*0.5d0*dxi
       dwdy(i,j,k)= (w(i,jp,k) - w(i,jm,k))*0.5d0*dxi
       dwdz(i,j,k)= (w(i,j,kp) - w(i,j,km))*0.5d0*dxi
    enddo
  enddo
enddo

epsilon=0.0d0

do k=1,nx
  do j=1,nx
    do i=1,nx
      S11 = dudx(i,j,k)
      S22 = dvdy(i,j,k)
      S33 = dwdz(i,j,k)
      S12 = 0.5d0*(dudy(i,j,k) + dvdx(i,j,k))
      S13 = 0.5d0*(dudz(i,j,k) + dwdx(i,j,k))
      S23 = 0.5d0*(dvdz(i,j,k) + dwdy(i,j,k))
      epsilon = epsilon + S11**2 + S22**2 + S33**2 + 2.d0*(S12**2 + S13**2 + S23**2)
    enddo
  enddo
enddo


epsilon=2*nu*epsilon/nx/nx/nx

write(*,*)  epsilon

return
end 
