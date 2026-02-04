subroutine get_interface(nstep)
  use param
  use flowvars
  implicit none
  integer :: nstep
  integer :: i,j,k,id,jd,kd
  integer :: top(nx,ny,nz),s_drop(nx,ny,nz)
  integer :: drop_count
 

  ! Binarize the Phase-Field
  do k=1,nz
    do j=1,ny
      do i=1,nx
        if(phi(i,j,k).ge.0.5d0)then
          top(i,j,k)=1
        else
          top(i,j,k)=0
        endif
      end do
    end do
  end do

  ! Initialize the count to zero
  drop_count=0

  ! flood fill algorithm
  do kd=1,nz
   do jd=1,ny
    do id=1,nx
     if(top(id,jd,kd).gt.0)then
      drop_count=drop_count+1
      write(*,'(2x,a,i3,a)') 'New drop, ',drop_count,' drops'
      ! single drop part
      s_drop=0
      s_drop(id,jd,kd)=1
      ! flood fill algorithm
      call flood_fill(top,s_drop,id,jd,kd)
      ! remove drops already done from top
      top=top-s_drop
      write(*,*) "loop"
      ! new drop calculation
      ! Compute diameter
      call calculate_deq(s_drop,nstep)
     endif
    enddo
   enddo
  enddo

  write(*,'(2x,a,i4)') 'Number of drops: ',drop_count
  write(*,*)
  
  open(2,file='drop_count.dat',access='append',form='formatted',status='old')
   write(2,'(i16,2x,i16)') nstep,drop_count
  close(2,status='keep')
  
  return
end

