program DropCountSize
! Serial program to compute the number and the size of the droplets
! Adaptation of the program 'mass_center' from FLOW36
use param
use flowvars
implicit none
integer :: ii

call read_input()

lx = twopi
dx = lx/(nx-1)

call print_start()

allocate(phi(nx,nx,nx))

! create output file
open(2,file='drop_count.dat',status='new',form='formatted')
 write(2,'(3(a16,2x))') 'iteration','t^+','drop count'
close(2,status='keep')


do ii=begin,finish,step
 write(*,'(1x,a,i8,a,i8)') 'Step ',ii,' out of ',finish
 call read_phi(ii)


 write(*,*) maxval(phi),minval(phi)

 call get_interface(ii)


end do


deallocate(phi)





end program DropCountSize
