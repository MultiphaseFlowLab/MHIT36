subroutine read_input
 use param 
 implicit none
 integer :: ios


open(unit=55,file='input.inp',form='formatted',status='old')
!Time step parameters
read(55,*) nx
read(55,*) begin
read(55,*) finish
read(55,*) step

rootpath='../../hit/output/'
  
close(1)

 return
end subroutine
!*********************************************
subroutine print_start()
 use param
 implicit none
 
 write(*,'(a)') repeat('-', 70)
 write(*,'(a)') '           Compute area'
 write(*,'(a)') repeat('-', 70)
 
 write(*,'(a,i6)') 'Grid size      : nx = ', nx
 write(*,'(a,f10.4)') 'Domain size  :  lx = ', lx

 write(*,'(a,i8,a,i8,a,i6)') 'Time steps     : from ', begin, ' to ', finish, ' with step ', step
 
 write(*,'(a)') repeat('-', 70)
 write(*,*)
 
 return
end subroutine print_start
!*********************************************
