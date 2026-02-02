program DropCountSize
! Serial program to compute the number and the size of the droplets
! Adaptation of the program 'mass_center' from FLOW36
use param
use flowvars
implicit none
integer :: ii

call read_input()

lx = 6.28
dx = lx/nx
ly = 3.14
dy = ly/ny
lz = 2.00
dz = lz/nz
epsr = 1.0

call print_start()

allocate(phi(nx,ny,nz))

! create output file
open(2,file='area.dat',status='new',form='formatted')
 write(2,'(3(a16,2x))') 'iteration','area'
close(2,status='keep')


do ii=begin,finish,step
 write(*,'(1x,a,i8,a,i8)') 'Step ',ii,' out of ',finish
 call read_phi(ii)


 write(*,*) maxval(phi),minval(phi)

 call get_interface(ii)


end do


deallocate(phi)





end program DropCountSize
