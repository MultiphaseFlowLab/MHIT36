subroutine flood_fill(top, s_drop, id, jd, kd)
  use param
  implicit none
  integer :: top(nx,ny,nz), s_drop(nx,ny,nz)
  integer :: id, jd, kd
  integer :: stack_i(nx*ny*nz), stack_j(nx*ny*nz), stack_k(nx*ny*nz)
  integer :: topi, topj, topk
  integer :: ni, nj, nk
  integer :: sp
  integer :: ddx(18), ddy(18), ddz(18)
  integer :: n

  ! 18-connectivity: 6 faces + 12 edges
  data ddx / &
      1,-1, 0, 0, 0, 0, &      ! faces
      1, 1,-1,-1, 1, 1,-1,-1, & ! edges (xy)
      1,-1, 0, 0 &             ! edges (xz & yz) (filled below)
  /
  data ddy / &
      0, 0, 1,-1, 0, 0, &      ! faces
      1,-1, 1,-1, 0, 0, 0, 0, & ! edges (xy & xz)
      0, 0, 1,-1 &
  /
  data ddz / &
      0, 0, 0, 0, 1,-1, &      ! faces
      0, 0, 0, 0, 1,-1, 1,-1, & ! edges (xz & yz)
      1,-1, 1,-1 &
  /

  ! If starting voxel is not pore: nothing to do
  if (top(id,jd,kd) /= 1) return

  ! Init stack
  sp = 1
  stack_i(sp) = id
  stack_j(sp) = jd
  stack_k(sp) = kd

  s_drop(id,jd,kd) = 1

  ! DFS using explicit stack
  do while (sp > 0)

    topi = stack_i(sp)
    topj = stack_j(sp)
    topk = stack_k(sp)
    sp = sp - 1

    do n = 1, 18
      ! periodic in X and Y
      ni = 1 + mod(topi - 1 + ddx(n) + nx, nx)
      nj = 1 + mod(topj - 1 + ddy(n) + ny, ny)

      ! Non-periodic in Z
      nk = topk + ddz(n)
      if (nk < 1 .or. nk > nz) cycle

      if (top(ni,nj,nk) == 1 .and. s_drop(ni,nj,nk) == 0) then
        s_drop(ni,nj,nk) = 1
        sp = sp + 1
        stack_i(sp) = ni
        stack_j(sp) = nj
        stack_k(sp) = nk
      end if
    end do

  end do

end subroutine flood_fill

