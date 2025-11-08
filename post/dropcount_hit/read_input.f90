subroutine read_input
 use param 
 implicit none
 integer :: ios

 open(unit=1, file='../in_post.txt', status='old', action='read', iostat=ios)
 if (ios /= 0) then
     print *, 'Error opening file.'
     stop
 end if
  
 read(1, *, iostat=ios) nx
 if (ios /= 0) then
     print *, 'Error reading number of points nx.'
     stop
 end if
  
 read(1, *, iostat=ios) begin
 if (ios /= 0) then
     print *, 'Error reading Begin.'
     stop
 end if
  
 read(1, *, iostat=ios) finish
 if (ios /= 0) then
     print *, 'Error reading End.'
     stop
 end if
  
 read(1, *, iostat=ios) step
 if (ios /= 0) then
     print *, 'Error reading Step.'
     stop
 end if

 read(1, *, iostat=ios) rootpath
 if (ios /= 0) then
     print *, 'Error reading Path to Results.'
     stop
 end if

 ! Example for 'in_post.txt'
 !512             ! Line 1: nx
 !0               ! Line 2: Begin (first saving to be processed)
 !140000             ! Line 3: End (last saving to be processed)
 !2000               ! Line 4: Step (step for processing: 1= every saving, 2=every other saving, ...)
 !"/leonardo_work/IscrB_SONORA/RUN1/multi/output"              !Line 4 : rootpath for results

 close(1)

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

 write(*,'(a,i8,a,i8,a,i6)') 'Time steps     : from ', begin, ' to ', finish, ' with step ', step
 
 write(*,'(a)') repeat('-', 70)
 write(*,*)
 
 return
end subroutine print_start
!*********************************************
