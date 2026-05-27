subroutine get_interface(nstep)
  use param
  use flowvars
  implicit none
  integer :: nstep
  integer :: i, j, k
  integer, allocatable :: top(:,:,:), label(:,:,:), stack(:)
  integer :: drop_count, cluster_size

  allocate(top(nx,nx,nx), label(nx,nx,nx), stack(nx*nx*nx))

  ! Binarize the phase field
  do k = 1, nx
    do j = 1, nx
      do i = 1, nx
        if (phi(i,j,k) >= 0.5d0) then
          top(i,j,k) = 1
        else
          top(i,j,k) = 0
        end if
      end do
    end do
  end do

  label      = 0
  drop_count = 0

  ! Connected-component labeling via iterative 26-connected flood fill.
  ! label==0 means unvisited; no per-cluster array resets needed.
  do k = 1, nx
    do j = 1, nx
      do i = 1, nx
        if (top(i,j,k) == 1 .and. label(i,j,k) == 0) then
          drop_count = drop_count + 1
          write(*,'(2x,a,i8,a)') 'New drop, ', drop_count, ' drops'
          call flood_fill_iter(top, label, stack, drop_count, i, j, k, cluster_size)
          call calculate_deq(cluster_size, nstep)
        end if
      end do
    end do
  end do

  write(*,'(2x,a,i8)') 'Number of drops: ', drop_count
  write(*,*)

  open(2, file='drop_count.dat', access='append', form='formatted', status='old')
    write(2,'(i16,2x,i16)') nstep, drop_count
  close(2, status='keep')

  deallocate(top, label, stack)

end subroutine get_interface
