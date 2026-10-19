! Tridiagonal solvers (TDMA) along z for the pressure Poisson equation, applied after the FFTs along x and y.
! For each (kx,ky) column, row k (k=1..nz, cell centers) of the system is:
!    a(k)*p(k-1) + (-a(k)-c(k)+lam)*p(k) + c(k)*p(k+1) = d(k),   lam = modified wavenumber (kx,ky)
! with Neumann BCs p(0)=p(1) and p(nz+1)=p(nz); for the (kx,ky)=(0,0) mode p(nz+1)=0 pins the pressure.
! The matrix does not change in time: all the factorizations are done once in the init routines (on the CPU).
! Two versions:
!  - tdma_init / tdma_solve:           z-pencil (kx,ky,z), each rank has full columns (needs the Y<->Z transposes)
!  - tdma_dist_init / tdma_dist_solve: y-pencil (ky,z,kx), z distributed over the ranks of zcomm (no Y<->Z transposes)
module tdma
use cudafor
use iso_c_binding
use param
implicit none

! coefficients of the system (both versions)
double precision, allocatable :: tdma_a(:), tdma_c(:)

! z-pencil TDMA: local sizes (kx,ky), 1/pivot and c' per column and row, top-BC factor per column
integer :: ntx, nty
double precision, allocatable :: tdma_inv(:,:,:), tdma_cp(:,:,:), tdma_ftop(:,:)

! distributed TDMA: y-pencil sizes (memory order y,z,kx) and offsets, z-communicator
integer :: ncy, nzy, nkxy, zoffy, kxoffy, zcomm, zrank, zsize, nchunk
! modified wavenumbers, interior Thomas factor c' and response p to the first row (per column and row)
double precision, allocatable :: lamx(:), lamy(:), dtdma_cp(:,:,:), dtdma_p(:,:,:)
! reduced system (2 rows per z-rank, 2*zsize rows per column): sub-diagonal, 1/pivot, c'
double precision, allocatable :: red_sub(:,:,:), red_inv(:,:,:), red_cp(:,:,:)
! RHS/solution buffers of the reduced system, exchanged with MPI_Alltoall: plain device memory (not managed),
! so that CUDA-aware MPI sends them directly from the GPU; red_sptr/red_rptr are host-side aliases of their
! device addresses, used only as MPI buffer arguments (never dereferenced on the host); c_devloc gives a c_devptr,
! converted with transfer to a plain c_ptr so that the standard c_f_pointer can build a host pointer
complex(8), device, allocatable :: red_sbuf(:,:,:), red_rbuf(:,:,:)
complex(8), pointer :: red_sptr(:), red_rptr(:)

contains

subroutine tdma_coeff
   integer :: k
   allocate(tdma_a(nz), tdma_c(nz))
   do k = 1, nz
      tdma_a(k) = 2.0d0*(dzi(k)**2.d0*dzi(k+1))/(dzi(k)+dzi(k+1))
      tdma_c(k) = 2.0d0*(dzi(k)*dzi(k+1)**2.d0)/(dzi(k)+dzi(k+1))
   enddo
end subroutine tdma_coeff




subroutine tdma_init
   ! z-pencil, memory order (kx,ky,z): forward elimination done once, layout (il,jl,k) so that
   ! consecutive GPU threads (il) access consecutive memory in the solve
   use mpi
   use mpivar, only: rank, ierr
   use cudecompvar, only: piZ_d2z
   integer :: il, jl, ig, jg, k, xoff, yoff
   double precision :: lam, den, cpk, atop, btop
   call tdma_coeff
   if (piZ_d2z%order(1) /= 1 .or. piZ_d2z%order(2) /= 2) then
      if (rank.eq.0) write(*,*) "tdma: z-pencil must be stored as (kx,ky,z)"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   endif
   ntx = piZ_d2z%shape(1)
   nty = piZ_d2z%shape(2)
   xoff = piZ_d2z%lo(1) - 1
   yoff = piZ_d2z%lo(2) - 1
   allocate(tdma_inv(ntx,nty,nz), tdma_cp(ntx,nty,nz), tdma_ftop(ntx,nty))
   do jl = 1, nty
      do il = 1, ntx
         ig = xoff + il
         jg = yoff + jl
         lam = 2.d0*(cos(kx(ig)*dx)-1.d0)*dxi*dxi + 2.d0*(cos(ky(jg)*dy)-1.d0)*dyi*dyi
         ! Neumann BC at bottom ghost (a=0, b=-1, c=1, d=0) -> c'(0) = -1, d'(0) = 0
         cpk = -1.d0
         do k = 1, nz
            den = -tdma_a(k) - tdma_c(k) + lam - tdma_a(k)*cpk
            tdma_inv(il,jl,k) = 1.d0/den
            cpk = tdma_c(k)/den
            tdma_cp(il,jl,k) = cpk
         enddo
         ! Neumann BC at top ghost (a=1, b=-1, c=0, d=0); pressure pinned for the (kx,ky)=(0,0) mode
         atop =  1.d0
         btop = -1.d0
         if (ig == 1 .and. jg == 1) then
            atop = 0.d0
            btop = 1.d0
         endif
         tdma_ftop(il,jl) = atop/(btop - atop*cpk)
      enddo
   enddo
end subroutine tdma_init




subroutine tdma_solve(psi)
   ! Thomas algorithm with precomputed coefficients, one thread per (kx,ky) column, solved in place
   complex(8), device :: psi(ntx,nty,nz)
   integer :: il, jl, k
   complex(8) :: dpk
   !$acc parallel loop collapse(2) gang vector private(dpk)
   do jl = 1, nty
      do il = 1, ntx
         ! Forward sweep: d'(k) = (d(k) - a(k)*d'(k-1)) / pivot(k), with d'(0) = 0 (bottom Neumann)
         dpk = (0.d0,0.d0)
         !$acc loop seq
         do k = 1, nz
            dpk = (psi(il,jl,k) - tdma_a(k)*dpk)*tdma_inv(il,jl,k)
            psi(il,jl,k) = dpk
         enddo
         ! Top ghost node: x(nz+1) = (d(nz+1) - a(nz+1)*d'(nz))/pivot(nz+1), d(nz+1) = 0
         dpk = -tdma_ftop(il,jl)*dpk
         ! Back substitution: x(k) = d'(k) - c'(k)*x(k+1)
         !$acc loop seq
         do k = nz, 1, -1
            dpk = psi(il,jl,k) - tdma_cp(il,jl,k)*dpk
            psi(il,jl,k) = dpk
         enddo
      enddo
   enddo
end subroutine tdma_solve




subroutine tdma_dist_init
   ! y-pencil, memory order (y,z,kx), z split over the ranks with the same kx range (zcomm).
   ! Each rank owns rows 1..nzy: the interior rows 2..nzy-1 are solved locally as x = z - x(1)*p - x(nzy)*q,
   ! with z = T_int^-1 d, T_int p = a(2) e_2, T_int q = c(nzy-1) e_nzy-1 (T_int = interior block).
   ! Rows 1 and nzy of all z-ranks then form a tridiagonal reduced system of 2*zsize rows per (kx,ky) column,
   ! which is solved by one z-rank per column (columns split in zsize chunks, exchanged with MPI_Alltoall).
   use mpi
   use mpivar, only: rank, ierr
   use cudecompvar, only: piY_d2z
   integer :: i, j, k, kg, col, dest, cc, src, row
   double precision :: lam, den, cpk, pk, qk, p2, q2, pm1, qm1, a1, b1, c1, am, bm, cm
   double precision, allocatable :: redc_s(:,:,:,:), redc_r(:,:,:,:)
   double precision, device, allocatable :: redc_sd(:,:,:,:), redc_rd(:,:,:,:)
   double precision, pointer :: redc_sptr(:), redc_rptr(:)
   call tdma_coeff
   if (piY_d2z%order(1) /= 2 .or. piY_d2z%order(2) /= 3) then
      if (rank.eq.0) write(*,*) "tdma_dist: y-pencil must be stored as (y,z,kx)"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   endif
   ncy = piY_d2z%shape(1)
   nzy = piY_d2z%shape(2)
   nkxy = piY_d2z%shape(3)
   zoffy = piY_d2z%lo(2) - 1
   kxoffy = piY_d2z%lo(3) - 1
   if (nzy < 3) then
      write(*,*) "tdma_dist: at least 3 z-planes per rank are needed, rank", rank, "has", nzy
      call mpi_abort(MPI_COMM_WORLD, 1, ierr)
   endif
   call mpi_comm_split(MPI_COMM_WORLD, kxoffy, zoffy, zcomm, ierr)
   call mpi_comm_rank(zcomm, zrank, ierr)
   call mpi_comm_size(zcomm, zsize, ierr)
   ! columns (j,i) -> linear index col = j + (i-1)*ncy, split in zsize chunks: chunk dest is solved by z-rank dest
   nchunk = (ncy*nkxy + zsize - 1)/zsize
   allocate(lamx(nkxy), lamy(ncy), dtdma_cp(ncy,nzy,nkxy), dtdma_p(ncy,nzy,nkxy))
   allocate(redc_s(3,2,nchunk,0:zsize-1), redc_r(3,2,nchunk,0:zsize-1))
   allocate(red_sub(2,nchunk,0:zsize-1), red_inv(2,nchunk,0:zsize-1), red_cp(2,nchunk,0:zsize-1))
   allocate(red_sbuf(2,nchunk,0:zsize-1), red_rbuf(2,nchunk,0:zsize-1))
   do i = 1, nkxy
      lamx(i) = 2.d0*(cos(kx(kxoffy+i)*dx)-1.d0)*dxi*dxi
   enddo
   do j = 1, ncy
      lamy(j) = 2.d0*(cos(ky(j)*dy)-1.d0)*dyi*dyi
   enddo
   ! padding columns of the last chunk: identity rows (sub, diag, sup) = (0, 1, 0)
   redc_s(1,:,:,:) = 0.d0
   redc_s(2,:,:,:) = 1.d0
   redc_s(3,:,:,:) = 0.d0
   red_sbuf = (0.d0,0.d0)
   red_rbuf = (0.d0,0.d0)
   call c_f_pointer(transfer(c_devloc(red_sbuf), c_null_ptr), red_sptr, [2*nchunk*zsize])
   call c_f_pointer(transfer(c_devloc(red_rbuf), c_null_ptr), red_rptr, [2*nchunk*zsize])
   do i = 1, nkxy
      do j = 1, ncy
         lam = lamx(i) + lamy(j)
         ! interior rows 2..nzy-1 (never a wall row): forward elimination, c' and p' (response to a(2) in row 2)
         cpk = 0.d0
         do k = 2, nzy-1
            kg = zoffy + k
            den = -tdma_a(kg) - tdma_c(kg) + lam - tdma_a(kg)*cpk
            if (k == 2) then
               pk = tdma_a(kg)/den
            else
               pk = -tdma_a(kg)*pk/den
            endif
            cpk = tdma_c(kg)/den
            dtdma_cp(j,k,i) = cpk
            dtdma_p(j,k,i) = pk
         enddo
         ! back substitution: p(k) = p'(k) - c'(k)*p(k+1); q(nzy-1) = c'(nzy-1), q(k) = -c'(k)*q(k+1)
         pm1 = dtdma_p(j,nzy-1,i)
         qm1 = dtdma_cp(j,nzy-1,i)
         pk = pm1
         qk = qm1
         do k = nzy-2, 2, -1
            pk = dtdma_p(j,k,i) - dtdma_cp(j,k,i)*pk
            dtdma_p(j,k,i) = pk
            qk = -dtdma_cp(j,k,i)*qk
         enddo
         p2 = pk
         q2 = qk
         ! first row of this rank (global zoffy+1), bottom Neumann on the first rank
         kg = zoffy + 1
         a1 = tdma_a(kg)
         c1 = tdma_c(kg)
         b1 = -a1 - c1 + lam
         if (kg == 1) then
            b1 = b1 + a1
            a1 = 0.d0
         endif
         ! last row of this rank (global zoffy+nzy), top Neumann (pinned for the (0,0) mode) on the last rank
         kg = zoffy + nzy
         am = tdma_a(kg)
         cm = tdma_c(kg)
         bm = -am - cm + lam
         if (kg == nz) then
            if (.not. (kxoffy+i == 1 .and. j == 1)) bm = bm + cm
            cm = 0.d0
         endif
         ! reduced rows (sub, diag, sup) for x(1) and x(nzy); unknowns ordered (x1,xm) of z-rank 0, 1, ...
         col = j + (i-1)*ncy
         dest = (col-1)/nchunk
         cc = col - dest*nchunk
         redc_s(1,1,cc,dest) = a1
         redc_s(2,1,cc,dest) = b1 - c1*p2
         redc_s(3,1,cc,dest) = -c1*q2
         redc_s(1,2,cc,dest) = -am*pm1
         redc_s(2,2,cc,dest) = bm - am*qm1
         redc_s(3,2,cc,dest) = cm
      enddo
   enddo
   ! gather the reduced rows of each column on the z-rank that solves it and factorize (Thomas)
   ! (exchanged through device copies: no MPI on managed buffers)
   allocate(redc_sd(3,2,nchunk,0:zsize-1), redc_rd(3,2,nchunk,0:zsize-1))
   redc_sd = redc_s
   call c_f_pointer(transfer(c_devloc(redc_sd), c_null_ptr), redc_sptr, [6*nchunk*zsize])
   call c_f_pointer(transfer(c_devloc(redc_rd), c_null_ptr), redc_rptr, [6*nchunk*zsize])
   call mpi_alltoall(redc_sptr, 6*nchunk, MPI_DOUBLE_PRECISION, redc_rptr, 6*nchunk, MPI_DOUBLE_PRECISION, zcomm, ierr)
   redc_r = redc_rd
   deallocate(redc_sd, redc_rd)
   do cc = 1, nchunk
      cpk = 0.d0
      do src = 0, zsize-1
         do row = 1, 2
            den = redc_r(2,row,cc,src) - redc_r(1,row,cc,src)*cpk
            cpk = redc_r(3,row,cc,src)/den
            red_sub(row,cc,src) = redc_r(1,row,cc,src)
            red_inv(row,cc,src) = 1.d0/den
            red_cp(row,cc,src) = cpk
         enddo
      enddo
   enddo
   deallocate(redc_s, redc_r)
   if (rank.eq.0) write(*,*) "Distributed TDMA: z-ranks per column", zsize, ", z-planes per rank", nzy
end subroutine tdma_dist_init




subroutine tdma_dist_solve(psi)
   ! Distributed Thomas algorithm with precomputed coefficients, one thread per (ky,kx) column, solved in place
   use mpi
   use mpivar, only: ierr
   complex(8), device :: psi(ncy,nzy,nkxy)
   integer :: i, j, k, kg, col, dest, cc, src, n, row
   double precision :: lam, den, cpk, qk
   complex(8) :: dpk, zf, zl, x1, xm
   ! (1) interior rows 2..nzy-1: z = T_int^-1 d (forward with on-the-fly pivots, backward with precomputed c'),
   !     then RHS of the reduced rows: d(1) - c(1)*z(2) and d(nzy) - a(nzy)*z(nzy-1)
   !$acc parallel loop collapse(2) gang vector private(lam, cpk, den, kg, dpk, zf, zl, col, dest, cc)
   do i = 1, nkxy
      do j = 1, ncy
         lam = lamx(i) + lamy(j)
         cpk = 0.d0
         dpk = (0.d0,0.d0)
         !$acc loop seq
         do k = 2, nzy-1
            kg = zoffy + k
            den = -tdma_a(kg) - tdma_c(kg) + lam - tdma_a(kg)*cpk
            cpk = tdma_c(kg)/den
            dpk = (psi(j,k,i) - tdma_a(kg)*dpk)/den
            psi(j,k,i) = dpk
         enddo
         zl = dpk
         !$acc loop seq
         do k = nzy-2, 2, -1
            dpk = psi(j,k,i) - dtdma_cp(j,k,i)*dpk
            psi(j,k,i) = dpk
         enddo
         zf = dpk
         col = j + (i-1)*ncy
         dest = (col-1)/nchunk
         cc = col - dest*nchunk
         red_sbuf(1,cc,dest) = psi(j,1,i)   - tdma_c(zoffy+1)*zf
         red_sbuf(2,cc,dest) = psi(j,nzy,i) - tdma_a(zoffy+nzy)*zl
      enddo
   enddo
   ! (2) reduced system: rows of each column gathered on one z-rank, Thomas with precomputed factors, sent back
   call mpi_alltoall(red_sptr, 2*nchunk, MPI_DOUBLE_COMPLEX, red_rptr, 2*nchunk, MPI_DOUBLE_COMPLEX, zcomm, ierr)
   !$acc parallel loop gang vector private(dpk, n, row, src)
   do cc = 1, nchunk
      dpk = (0.d0,0.d0)
      !$acc loop seq
      do n = 1, 2*zsize
         row = mod(n-1,2) + 1
         src = (n-1)/2
         dpk = (red_rbuf(row,cc,src) - red_sub(row,cc,src)*dpk)*red_inv(row,cc,src)
         red_rbuf(row,cc,src) = dpk
      enddo
      !$acc loop seq
      do n = 2*zsize-1, 1, -1
         row = mod(n-1,2) + 1
         src = (n-1)/2
         dpk = red_rbuf(row,cc,src) - red_cp(row,cc,src)*dpk
         red_rbuf(row,cc,src) = dpk
      enddo
   enddo
   call mpi_alltoall(red_rptr, 2*nchunk, MPI_DOUBLE_COMPLEX, red_sptr, 2*nchunk, MPI_DOUBLE_COMPLEX, zcomm, ierr)
   ! (3) local reconstruction: x(k) = z(k) - x(1)*p(k) - x(nzy)*q(k), q(nzy-1) = c'(nzy-1), q(k) = -c'(k)*q(k+1)
   !$acc parallel loop collapse(2) gang vector private(x1, xm, qk, col, dest, cc)
   do i = 1, nkxy
      do j = 1, ncy
         col = j + (i-1)*ncy
         dest = (col-1)/nchunk
         cc = col - dest*nchunk
         x1 = red_sbuf(1,cc,dest)
         xm = red_sbuf(2,cc,dest)
         psi(j,1,i) = x1
         psi(j,nzy,i) = xm
         qk = -1.d0
         !$acc loop seq
         do k = nzy-1, 2, -1
            qk = -dtdma_cp(j,k,i)*qk
            psi(j,k,i) = psi(j,k,i) - x1*dtdma_p(j,k,i) - xm*qk
         enddo
      enddo
   enddo
end subroutine tdma_dist_solve

end module tdma
