program read_to_paraview
use commondata
implicit none

integer :: i,j,k
logical :: check
character(len=40) :: namefile

!! read input
open(10,file='input_par.inp',form='formatted')
read(10,*) nx
read(10,*) nstart
read(10,*) nend
read(10,*) dump
read(10,*) uflag
read(10,*) vflag
read(10,*) wflag
read(10,*) phiflag

!! overide (if there are problems)
!nstart=0
!dump=200
!nend=200
!nx=64

ny=nx
nz=nx

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

dx=2*pi/(nx-1)
dy=2*pi/(ny-1)
dz=2*pi/(nz-1)

x=0.0d0
y=0.0d0
z=0.0d0

do i=1,nx-1
  x(i+1)=x(i)+dx
enddo
do j=1,ny-1
  y(j+1)=y(j)+dy
enddo
do k=1,nz-1
  z(k+1)=z(k)+dz
enddo

! read fluid data
do i=nstart,nend,dump
 call read_fields(i)
enddo

deallocate(x,y,z)

end program read_to_paraview
