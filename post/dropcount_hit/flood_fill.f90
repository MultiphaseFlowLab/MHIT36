! Iterative 26-connected flood fill using an explicit DFS stack.
! Flat index encoding: flat = (k-1)*nx*nx + (j-1)*nx + (i-1), range [0, nx^3-1].
! stack workspace (size >= nx*nx*nx) is allocated once by the caller and reused.
subroutine flood_fill_iter(top, label, stack, current_label, i0, j0, k0, cluster_size)
  use param
  implicit none
  integer, intent(in)    :: top(nx,nx,nx)
  integer, intent(inout) :: label(nx,nx,nx)
  integer, intent(inout) :: stack(*)
  integer, intent(in)    :: current_label, i0, j0, k0
  integer, intent(out)   :: cluster_size

  integer :: sp, flat, i, j, k, ni, nj, nk, d

  ! 26-connectivity offsets: 6 faces + 12 edges + 8 corners
  integer, parameter :: ndirs = 26
  integer :: ox(ndirs), oy(ndirs), oz(ndirs)

  data ox / 1,-1, 0, 0, 0, 0,                    &  ! faces
            1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,  &  ! edges
            1, 1, 1, 1,-1,-1,-1,-1 /                 ! corners
  data oy / 0, 0, 1,-1, 0, 0,                    &
            1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,  &
            1, 1,-1,-1, 1, 1,-1,-1 /
  data oz / 0, 0, 0, 0, 1,-1,                    &
            0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,  &
            1,-1, 1,-1, 1,-1, 1,-1 /

  label(i0,j0,k0) = current_label
  sp           = 1
  stack(1)     = (k0-1)*nx*nx + (j0-1)*nx + (i0-1)
  cluster_size = 0

  do while (sp > 0)
    flat = stack(sp)
    sp   = sp - 1
    cluster_size = cluster_size + 1

    ! decode flat (0-based) to 1-based indices
    i = mod(flat,       nx) + 1
    j = mod(flat/nx,    nx) + 1
    k =     flat/(nx*nx)   + 1

    do d = 1, ndirs
      ni = mod(i - 1 + ox(d) + nx, nx) + 1
      nj = mod(j - 1 + oy(d) + nx, nx) + 1
      nk = mod(k - 1 + oz(d) + nx, nx) + 1
      if (top(ni,nj,nk) == 1 .and. label(ni,nj,nk) == 0) then
        label(ni,nj,nk) = current_label
        sp       = sp + 1
        stack(sp) = (nk-1)*nx*nx + (nj-1)*nx + (ni-1)
      end if
    end do
  end do

end subroutine flood_fill_iter
