subroutine flood_fill(top, s_drop, id, jd, kd)
  use param, only: nx
  implicit none
  integer, intent(in)    :: top(nx,nx,nx)
  integer, intent(inout) :: s_drop(nx,nx,nx)
  integer, intent(in)    :: id, jd, kd
  integer :: stack_i(nx*nx*nx), stack_j(nx*nx*nx), stack_k(nx*nx*nx)
  integer :: topi, topj, topk
  integer :: ni, nj, nk
  integer :: sp

  integer :: dx, dy, dz
  integer :: manh

  ! quick rejects
  if (id < 1 .or. id > nx) return
  if (jd < 1 .or. jd > nx) return
  if (kd < 1 .or. kd > nx) return
  if (top(id,jd,kd) /= 1) return
  if (s_drop(id,jd,kd) /= 0) return

  ! init explicit stack (DFS)
  sp = 1
  stack_i(sp) = id
  stack_j(sp) = jd
  stack_k(sp) = kd
  s_drop(id,jd,kd) = 1

  do while (sp > 0)

    topi = stack_i(sp)
    topj = stack_j(sp)
    topk = stack_k(sp)
    sp = sp - 1

    ! 18-connectivity: |dx|+|dy|+|dz| = 1 or 2 (exclude corners where it is 3)
    do dz = -1, 1
      do dy = -1, 1
        do dx = -1, 1

          if (dx == 0 .and. dy == 0 .and. dz == 0) cycle

          manh = abs(dx) + abs(dy) + abs(dz)
          if (manh > 2) cycle

          ! periodic in X, Y, Z (1-based indices)
          ni = 1 + mod(topi - 1 + dx + nx, nx)
          nj = 1 + mod(topj - 1 + dy + nx, nx)
          nk = 1 + mod(topk - 1 + dz + nx, nx)

          if (top(ni,nj,nk) == 1 .and. s_drop(ni,nj,nk) == 0) then
            s_drop(ni,nj,nk) = 1
            sp = sp + 1
            if (sp > nx*nx*nx) stop "flood_fill: explicit stack overflow"
            stack_i(sp) = ni
            stack_j(sp) = nj
            stack_k(sp) = nk
          end if

        end do
      end do
    end do

  end do

end subroutine flood_fill

