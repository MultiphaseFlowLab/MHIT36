subroutine read_input
 use param 
 implicit none
 integer :: ios


open(unit=55,file='input.inp',form='formatted',status='old')
!Time step parameters
read(55,*) nx
read(55,*) ny
read(55,*) nz
read(55,*) begin
read(55,*) finish
read(55,*) step
read(55,*) csi

rootpath='../../channel/output/'

close(55)

 return
end subroutine
!*********************************************
subroutine print_start()
 use param
 implicit none
 
 write(*,'(a)') repeat('-', 70)
 write(*,'(a)') '           Count Droplets and Compute Their Size'
 write(*,'(a)') repeat('-', 70)
 
 write(*,'(a,i6)') 'Grid size      : nx = ', nx
 write(*,'(a,f10.4)') 'Domain size  :  lx = ', lx
 write(*,'(a,i6)') 'Grid size      : ny = ', ny
 write(*,'(a,f10.4)') 'Domain size  :  ly = ', ly
 write(*,'(a,i6)') 'Grid size      : nz = ', nz
 write(*,'(a,f10.4)') 'Domain size  :  lz = ', lz
 write(*,'(a,f10.4)') 'Z-stretch      : csi = ', csi

 write(*,'(a,i8,a,i8,a,i6)') 'Time steps     : from ', begin, ' to ', finish, ' with step ', step
 
 write(*,'(a)') repeat('-', 70)
 write(*,*)
 
 return
end subroutine print_start
!*********************************************
subroutine build_grid()
 ! stretched grid along z, same as channel/readinput.f90
 use param
 implicit none
 integer :: k
 double precision :: zk

 allocate(z(0:nz+1),dzci(nz))
 do k = 1, nz
   zk=(dble(k)-0.5d0)/dble(nz)
   z(k) = 0.5d0*lz*(1.d0+tanh(csi*(zk-0.5d0))/tanh(0.5d0*csi))
 enddo
 z(0)=-z(1)
 z(nz+1)= lz+(lz-z(nz))
 ! inverse of the cell height (between cell faces)
 dzci(1) = 2.d0/(z(1)+ z(2))
 dzci(nz) = 1.d0/(lz-0.5d0*(z(nz)+z(nz-1)))
 do k = 2, nz-1
    dzci(k) = 2.d0/(z(k+1)-z(k-1))
 enddo

 return
end subroutine build_grid
!*********************************************
