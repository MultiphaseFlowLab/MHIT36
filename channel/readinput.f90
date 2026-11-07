
!##########################################################################
!###########################################################################
subroutine readinput
use velocity
use phase
use temperature
use param
use mpivar
implicit none
integer :: i,j,k
double precision :: zk

open(unit=55,file='input.inp',form='formatted',status='old')
!Time step parameters
read(55,*) restart
read(55,*) tstart
read(55,*) tfin
read(55,*) dump
! Domain size
read(55,*) lx
read(55,*) ly
read(55,*) lz
read(55,*) csi
!Flow parameters
read(55,*) inflow
read(55,*) inphi
read(55,*) intheta
read(55,*) dt
read(55,*) mu
read(55,*) rho
! forcing parameters
read(55,*) gradpx
read(55,*) gradpy
! phase-field parameters
read(55,*) radius
read(55,*) sigma
read(55,*) epsr
! level-set parameters
read(55,*) epsadvr
read(55,*) nreinit
read(55,*) tolreinit
! temperature parameters
read(55,*) kappa
read(55,*) alphag
close(55)

! compute pre-defined constants
twopi=8.0_8*atan(1.0_8)
pi=twopi/2.d0
dx = lx/nx
dy = ly/ny
dz = lz/nz
dxi = 1.d0/dx
dyi = 1.d0/dy
ddxi = 1.d0/dx/dx
ddyi = 1.d0/dy/dy
ddzi = 1.d0/dz/dz
rhoi=1.d0/rho
eps=epsr*dx
epsi=1.d0/eps
enum=1.e-16
!write(*,*) "Check on stability", dt*mu*dzi*dzi

if (rank .eq. 0) then
    !enable/disable for debug check parameters
    write(*,*) "----------------------------------------------"
    write(*,*) "███    ███ ██   ██ ██ ████████ ██████   ██████"  
    write(*,*) "████  ████ ██   ██ ██    ██         ██ ██"       
    write(*,*) "██ ████ ██ ███████ ██    ██     █████  ███████"  
    write(*,*) "██  ██  ██ ██   ██ ██    ██         ██ ██    ██" 
    write(*,*) "██      ██ ██   ██ ██    ██    ██████   ██████"         
    write(*,*) "----------------------------------------------"
    write(*,*) "-------------Channel flow setup---------------"
    write(*,*) 'Grid:', nx, 'x', ny, 'x', nz
    write(*,*) "Restart   ", restart
    write(*,*) "Tstart    ", tstart
    write(*,*) "Tfin      ", tfin
    write(*,*) "Dump      ", dump
    write(*,*) "Inflow    ", inflow
    write(*,*) "Deltat    ", dt
    write(*,*) "Mu        ", mu
    write(*,*) "Rho       ", rho
    write(*,*) "Gradpx    ", gradpx
    write(*,*) "Gradpy    ", gradpy
    write(*,*) "Radius    ", radius
    write(*,*) "Sigma     ", sigma
    write(*,*) "Eps       ", eps
    write(*,*) "Epsadvr   ", epsadvr
    write(*,*) "Nreinit   ", nreinit
    write(*,*) "Tolreinit ", tolreinit
    write(*,*) "Kappa     ", kappa
    write(*,*) "Alphag    ", alphag
    write(*,*) "Lx        ", lx
    write(*,*) "Ly        ", ly
    write(*,*) "Lz        ", lz
    write(*,*) "Z-stretch ", csi
!    write(*,*) "dx", dx
!    write(*,*) "dxi", dxi
!    write(*,*) "ddxi", ddxi
!    write(*,*) "dy", dx
!    write(*,*) "dyi", dyi
!    write(*,*) "ddyi", ddyi
!    write(*,*) "dz", dz
!    write(*,*) "dzi", dzi
!    write(*,*) "ddzi", ddzi
endif



! define wavenumbers and grid points axis
! define grid (then move in readinput)
allocate(x(nx),y(ny),z(0:nz+1),dzci(nz),dzi(nz+1),kx(nx),ky(ny))
! location of the pressure nodes (cell centers)
x(1)=dx/2
do i = 2, nx
   x(i) = x(i-1) + dx
enddo
y(1)=dy/2
do j = 2, ny
   y(j) = y(j-1) + dy
enddo
! stretched grid along z; z axis include also the two ghost nodes located at +/- dz/2 above and below the wall
do k = 1, nz
  zk=(dble(k)-0.5d0)/dble(nz)         
  z(k) = 0.5d0*dble(lz)*(1.d0+tanh(csi*(zk-0.5d0))/tanh(0.5d0*csi))
enddo
z(0)=-z(1)
z(nz+1)= lz+(lz-z(nz))
! compute inverse of dz (between cell faces)
dzci(1) = 2.d0/(z(1)+ z(2))
dzci(nz) = 1.d0/(lz-0.5d0*(z(nz)+z(nz-1)))
do k = 2, nz-1
   dzci(k) = 2.d0/(z(k+1)-z(k-1))
enddo
! compute inverse of dz (between nodes)
dzi(1)=0.5d0/z(1)
!dzi(nz+1)=0.5d0/(lz-z(nz))
do k=1, nz+1
   dzi(k) = 1.d0/(z(k)-z(k-1))
enddo
!write(*,*) "dzi", dzi
! quantities for the level-set boundedness constraints (Wang et al. 2025, extended to the stretched z grid)
! dzmax: largest spacing between interior nodes; lapmax: largest diagonal coefficient of the discrete laplacian
dzmax = 0.d0
do k = 2, nz
   dzmax = max(dzmax, 1.d0/dzi(k))
enddo
lapmax = 0.d0
do k = 1, nz
   lapmax = max(lapmax, 2.d0*ddxi + 2.d0*ddyi + dzci(k)*(dzi(k)+dzi(k+1)))
enddo
if (rank .eq. 0) then
   write(*,*) "Dzmax     ", dzmax
   ! lower bound on the reinitialization diffusion: eps >= max(dx,dy,dz)/2
   if (eps .lt. 0.5d0*max(dx,max(dy,dzmax))) write(*,*) "Warning: eps < max(dx,dy,dz)/2, phi may be unbounded (increase epsr)"
   if (epsadvr .lt. 0.5d0) write(*,*) "Warning: epsadvr < 0.5, phi may be unbounded"
endif
! wavenumber
do i = 1, nx/2
   kx(i) = (i-1)*(twopi/lx)
enddo
do i = nx/2+1, nx
   kx(i) = (i-1-nx)*(twopi/lx)
enddo
do j = 1, ny/2
   ky(j) = (j-1)*(twopi/ly)
enddo
do j = ny/2+1, ny
   ky(j) = (j-1-ny)*(twopi/ly)
enddo
! allocate kx_d and ky_d on the device 
allocate(kx_d, source=kx)
allocate(ky_d, source=ky)



end subroutine


