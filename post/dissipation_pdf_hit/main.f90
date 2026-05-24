! ======================================================================
! main.f90
! ----------------------------------------------------------------------
! Computes the local viscous dissipation
!     eps(x) = 2 * nu(phi) * s_ij s_ij,    s_ij = 1/2 (du_i/dx_j + du_j/dx_i)
! on a staggered (MAC) grid in triply-periodic HIT, then builds PDFs
! conditioned on the phase field phi (inside drops / interface / outside).
!
! Staggered layout assumed:
!     u(i,j,k)   at face (i-1/2, j,   k  )    -- LEFT face
!     v(i,j,k)   at face (i,   j-1/2, k  )    -- BOTTOM face
!     w(i,j,k)   at face (i,   j,   k-1/2)    -- BACK face
!     phi(i,j,k) at cell center (i, j, k)
!
! Output, one set per snapshot, gnuplot-ready columns:
!     eps_pdf_inside_NNNNNNNN.dat
!     eps_pdf_interface_NNNNNNNN.dat
!     eps_pdf_outside_NNNNNNNN.dat
!     eps_pdf_total_NNNNNNNN.dat
! plus a single eps_mean.dat with snapshot-by-snapshot means.
!
! Build (gfortran):
!     gfortran -O3 -march=native -fopenmp main.f90 -o compute_diss
! ======================================================================

module commondata
  integer :: nx
  integer :: nstart, nend, ndump, sdump
  double precision, parameter :: pi = 3.14159265358979d0
  double precision :: re, dt, dx, lx, dxi, epsilon, nu
  double precision, allocatable, dimension(:,:,:) :: u, v, w
  double precision, allocatable, dimension(:,:,:) :: dudx, dudy, dudz
  double precision, allocatable, dimension(:,:,:) :: dvdx, dvdy, dvdz
  double precision, allocatable, dimension(:,:,:) :: dwdx, dwdy, dwdz
  double precision, allocatable, dimension(:) :: x

  ! Path containing input fields (no trailing slash issues; we add '/')
  character(len=*), parameter :: input_dir = '../../hit/output'
end module commondata


module hist_params
  implicit none
  ! Histogram resolution and conditioning thresholds
  integer, parameter :: nbins_lin = 1000
  integer, parameter :: nbins_log = 1000
  double precision, parameter :: phi_in_thr  = 0.95d0  ! phi >= -> inside drop
  double precision, parameter :: phi_out_thr = 0.05d0  ! phi <= -> outside
  double precision, parameter :: eps_floor   = 1.0d-20 ! excluded from log PDF
  ! Phase convention: .true. if phi=1 inside drops, .false. if phi=0 inside.
  logical, parameter :: phi_drop_is_one = .true.
  ! Two-viscosity setup. For matched case set nu_drop = nu_matrix = nu.
  double precision, parameter :: nu_drop   = 0.006d0
  double precision, parameter :: nu_matrix = 0.006d0
end module hist_params


program compute_diss_staggered
  use commondata
  use hist_params
  implicit none
  integer :: nstep, unit_mean

  ! Geometry / physics
  nx    = 512
  lx    = 2.0d0 * pi
  dx    = lx / dble(nx)            ! strictly periodic
  dxi   = 1.0d0 / dx
  nu    = nu_matrix                ! kept for backward compat with old code

  ! Snapshot sweep
  nstart = 200000
  ndump  = 5000
  nend   = 200000

  write(*,'(a,i0,a,i0,a,i0)') '[info] grid: ', nx, ' x ', nx, ' x ', nx
  write(*,'(a,es12.5)') '[info] dx  = ', dx
  write(*,'(a,es12.5,a,es12.5)') '[info] nu_drop = ', nu_drop, &
       '   nu_matrix = ', nu_matrix
  write(*,'(a,a)')      '[info] input_dir = ', trim(input_dir)
  write(*,'(a,i0,a,i0,a,i0)') '[info] snapshots: nstart=', nstart, &
       ' ndump=', ndump, ' nend=', nend

  ! Allocate (matches commondata layout; phi and eps added locally below)
  allocate(u(nx,nx,nx), v(nx,nx,nx), w(nx,nx,nx))
  allocate(dudx(nx,nx,nx), dudy(nx,nx,nx), dudz(nx,nx,nx))
  allocate(dvdx(nx,nx,nx), dvdy(nx,nx,nx), dvdz(nx,nx,nx))
  allocate(dwdx(nx,nx,nx), dwdy(nx,nx,nx), dwdz(nx,nx,nx))

  ! Open the cross-snapshot summary file
  open(newunit=unit_mean, file='eps_mean.dat', status='replace', action='write')
  write(unit_mean,'(a)') '# snapshot  <eps>_total  <eps>_in  <eps>_if  <eps>_out  vf_in  vf_if  vf_out  max|div u|'

  do nstep = nstart, nend, ndump
     call process_snapshot(nstep, unit_mean)
  end do

  close(unit_mean)
  deallocate(u, v, w)
  deallocate(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
  write(*,'(a)') '[done]'

end program compute_diss_staggered


! ----------------------------------------------------------------------
subroutine process_snapshot(nstep, unit_mean)
  use commondata
  use hist_params
  implicit none
  integer, intent(in) :: nstep, unit_mean

  double precision, allocatable :: phi(:,:,:), eps(:,:,:)
  double precision :: eps_min, eps_max, eps_mean, eps_in, eps_if, eps_out
  double precision :: vf_in, vf_if, vf_out, max_div, div_loc
  double precision :: lin_lo, lin_hi, log_lo, log_hi
  integer(kind=8)  :: n_in, n_if, n_out, n_tot, n_pos_in, n_pos_if, n_pos_out, n_pos_tot
  integer(kind=8)  :: hl_in(nbins_lin),  hl_if(nbins_lin),  hl_out(nbins_lin),  hl_tot(nbins_lin)
  integer(kind=8)  :: hg_in(nbins_log),  hg_if(nbins_log),  hg_out(nbins_log),  hg_tot(nbins_log)
  integer :: i, j, k, ip, jp, kp, im, jm, km
  double precision :: nu_loc, sij_sij_loc, s11, s22, s33, s12, s13, s23
  double precision :: phi_d, eps_loc, sum_e, sum_in, sum_if, sum_out

  allocate(phi(nx,nx,nx), eps(nx,nx,nx))

  write(*,'(a,i0)') '[step] reading snapshot ', nstep
  call read_fields_snap(nstep, phi)

  ! --------------------------------------------------------------------
  ! 1) Cell-centered velocity gradient tensor on staggered grid
  ! --------------------------------------------------------------------
  !
  ! u(i,j,k) lives on the LEFT  face  (i-1/2, j,   k  )
  ! v(i,j,k) lives on the BOTTOM face (i,   j-1/2, k  )
  ! w(i,j,k) lives on the BACK   face (i,   j,   k-1/2)
  ! so the cell-centered diagonal derivatives are FORWARD differences:
  !   du/dx|_c(i,j,k) = (u(i+1,j,k) - u(i,j,k)) / dx
  !   dv/dy|_c(i,j,k) = (v(i,j+1,k) - v(i,j,k)) / dy
  !   dw/dz|_c(i,j,k) = (w(i,j,k+1) - w(i,j,k)) / dz
  !
  !$omp parallel do default(none) shared(u,v,w,dudx,dvdy,dwdz,nx,dxi) &
  !$omp private(i,j,k,ip,jp,kp)
  do k = 1, nx
     kp = k + 1; if (kp > nx) kp = 1
     do j = 1, nx
        jp = j + 1; if (jp > nx) jp = 1
        do i = 1, nx
           ip = i + 1; if (ip > nx) ip = 1
           dudx(i,j,k) = (u(ip,j ,k ) - u(i,j,k)) * dxi
           dvdy(i,j,k) = (v(i ,jp,k ) - v(i,j,k)) * dxi
           dwdz(i,j,k) = (w(i ,j ,kp) - w(i,j,k)) * dxi
        end do
     end do
  end do
  !$omp end parallel do

  ! Off-diagonals: store the edge derivatives first, then average.
  ! With u, v, w on the minus-half faces, the natural edge derivative is the
  ! BACKWARD difference, e.g. (u(i,j,k) - u(i,j-1,k)) / dy lives at
  ! the edge (i-1/2, j-1/2, k). The 4-point average to cell center then
  ! looks at index offsets +0 and +1 along face_axis and deriv_axis.

  !$omp parallel do default(none) shared(u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy,nx,dxi) &
  !$omp private(i,j,k,im,jm,km)
  do k = 1, nx
     km = k - 1; if (km < 1) km = nx
     do j = 1, nx
        jm = j - 1; if (jm < 1) jm = nx
        do i = 1, nx
           im = i - 1; if (im < 1) im = nx
           dudy(i,j,k) = (u(i ,j ,k ) - u(i ,jm,k )) * dxi
           dudz(i,j,k) = (u(i ,j ,k ) - u(i ,j ,km)) * dxi
           dvdx(i,j,k) = (v(i ,j ,k ) - v(im,j ,k )) * dxi
           dvdz(i,j,k) = (v(i ,j ,k ) - v(i ,j ,km)) * dxi
           dwdx(i,j,k) = (w(i ,j ,k ) - w(im,j ,k )) * dxi
           dwdy(i,j,k) = (w(i ,j ,k ) - w(i ,jm,k )) * dxi
        end do
     end do
  end do
  !$omp end parallel do

  ! 4-point average of each edge derivative back to the cell center.
  call average_edges_to_center(dudy, 1, 2)   ! du/dy: face_axis=1 (u on x-face), deriv along y (axis 2)
  call average_edges_to_center(dudz, 1, 3)
  call average_edges_to_center(dvdx, 2, 1)
  call average_edges_to_center(dvdz, 2, 3)
  call average_edges_to_center(dwdx, 3, 1)
  call average_edges_to_center(dwdy, 3, 2)

  ! --------------------------------------------------------------------
  ! 2) Build eps and accumulate moments + bin range
  ! --------------------------------------------------------------------
  max_div  = 0.0d0
  eps_min  = 0.0d0
  eps_max  = 100.d0
  sum_e    = 0.0d0; sum_in = 0.0d0; sum_if = 0.0d0; sum_out = 0.0d0
  n_in     = 0; n_if = 0; n_out = 0; n_tot = 0

  !$omp parallel do default(none) &
  !$omp shared(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,phi,eps,nx) &
  !$omp private(i,j,k,s11,s22,s33,s12,s13,s23,sij_sij_loc,nu_loc,phi_d,eps_loc,div_loc) &
  !$omp reduction(max:max_div,eps_max) reduction(min:eps_min) &
  !$omp reduction(+:sum_e,sum_in,sum_if,sum_out,n_in,n_if,n_out,n_tot)
  do k = 1, nx
     do j = 1, nx
        do i = 1, nx
           s11 = dudx(i,j,k)
           s22 = dvdy(i,j,k)
           s33 = dwdz(i,j,k)
           s12 = 0.5d0 * (dudy(i,j,k) + dvdx(i,j,k))
           s13 = 0.5d0 * (dudz(i,j,k) + dwdx(i,j,k))
           s23 = 0.5d0 * (dvdz(i,j,k) + dwdy(i,j,k))
           sij_sij_loc = s11*s11 + s22*s22 + s33*s33 &
                       + 2.0d0 * (s12*s12 + s13*s13 + s23*s23)

           ! Local viscosity: linear mixture rule.
           if (phi_drop_is_one) then
              phi_d = phi(i,j,k)
           else
              phi_d = 1.0d0 - phi(i,j,k)
           end if
           if (phi_d < 0.0d0) phi_d = 0.0d0
           if (phi_d > 1.0d0) phi_d = 1.0d0
           nu_loc  = nu_drop * phi_d + nu_matrix * (1.0d0 - phi_d)

           eps_loc = 2.0d0 * nu_loc * sij_sij_loc
           eps(i,j,k) = eps_loc

           ! Stats
           sum_e = sum_e + eps_loc
           n_tot = n_tot + 1
           if (eps_loc < eps_min) eps_min = eps_loc
           if (eps_loc > eps_max) eps_max = eps_loc

           if (phi_d >= phi_in_thr) then
              n_in = n_in + 1
              sum_in = sum_in + eps_loc
           else if (phi_d <= phi_out_thr) then
              n_out = n_out + 1
              sum_out = sum_out + eps_loc
           else
              n_if = n_if + 1
              sum_if = sum_if + eps_loc
           end if

           ! Divergence check
           div_loc = abs(s11 + s22 + s33)
           if (div_loc > max_div) max_div = div_loc
        end do
     end do
  end do
  !$omp end parallel do

  eps_mean = sum_e / dble(n_tot)
  vf_in  = dble(n_in)  / dble(n_tot)
  vf_if  = dble(n_if)  / dble(n_tot)
  vf_out = dble(n_out) / dble(n_tot)
  if (n_in  > 0) then; eps_in  = sum_in  / dble(n_in);  else; eps_in  = 0.0d0; end if
  if (n_if  > 0) then; eps_if  = sum_if  / dble(n_if);  else; eps_if  = 0.0d0; end if
  if (n_out > 0) then; eps_out = sum_out / dble(n_out); else; eps_out = 0.0d0; end if

  write(*,'(a,i0,a,es11.4,a,es11.4,a,es11.4,a,es11.4)') &
       '[step ', nstep, '] <eps>=', eps_mean, ' in=', eps_in, &
       ' if=', eps_if, ' out=', eps_out
  write(*,'(a,f6.3,a,f6.3,a,f6.3,a,es10.3)') &
       '          vf_in=', vf_in, ' vf_if=', vf_if, ' vf_out=', vf_out, &
       ' max|div|=', max_div

  ! Cross-snapshot one-liner
  write(unit_mean,'(i8,1x,es16.8,1x,es16.8,1x,es16.8,1x,es16.8,1x,f9.6,1x,f9.6,1x,f9.6,1x,es10.3)') &
       nstep, eps_mean, eps_in, eps_if, eps_out, vf_in, vf_if, vf_out, max_div

  ! --------------------------------------------------------------------
  ! 3) Build histograms for this snapshot
  ! --------------------------------------------------------------------
  lin_lo = 0.0d0
  lin_hi = eps_max * 1.001d0
  if (lin_hi <= lin_lo) lin_hi = lin_lo + 1.0d0

  if (eps_min > eps_floor) then
     log_lo = log10(eps_min) - 0.5d0
  else
     log_lo = log10(eps_floor)
  end if
  log_hi = log10(eps_max) + 0.1d0
  if (log_hi <= log_lo) log_hi = log_lo + 1.0d0

  hl_in = 0; hl_if = 0; hl_out = 0; hl_tot = 0
  hg_in = 0; hg_if = 0; hg_out = 0; hg_tot = 0
  n_pos_in = 0; n_pos_if = 0; n_pos_out = 0; n_pos_tot = 0

  call build_hist(eps, phi, &
                  lin_lo, lin_hi, log_lo, log_hi, &
                  hl_in, hl_if, hl_out, hl_tot, &
                  hg_in, hg_if, hg_out, hg_tot, &
                  n_pos_in, n_pos_if, n_pos_out, n_pos_tot)

  ! --------------------------------------------------------------------
  ! 4) Write PDF files
  ! --------------------------------------------------------------------
  call write_pdf_file(nstep, 'inside',    hl_in,  hg_in,  n_in,  n_pos_in,  &
                      lin_lo, lin_hi, log_lo, log_hi, eps_in,  sum_in)
  call write_pdf_file(nstep, 'interface', hl_if,  hg_if,  n_if,  n_pos_if,  &
                      lin_lo, lin_hi, log_lo, log_hi, eps_if,  sum_if)
  call write_pdf_file(nstep, 'outside',   hl_out, hg_out, n_out, n_pos_out, &
                      lin_lo, lin_hi, log_lo, log_hi, eps_out, sum_out)
  call write_pdf_file(nstep, 'total',     hl_tot, hg_tot, n_tot, n_pos_tot, &
                      lin_lo, lin_hi, log_lo, log_hi, eps_mean, sum_e)

  deallocate(phi, eps)
end subroutine process_snapshot


! ----------------------------------------------------------------------
subroutine read_fields_snap(nstep, phi)
  use commondata
  implicit none
  integer, intent(in) :: nstep
  double precision, intent(out) :: phi(nx,nx,nx)
  character(len=512) :: fname
  integer :: io, ierr

  write(fname,'(a,a,a,i8.8,a)') trim(input_dir), '/', 'u_',   nstep, '.dat'
open(newunit=io, file=trim(fname), form='unformatted', &
     access='stream', status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) '[error] cannot open ', trim(fname); stop 1
  end if
  read(io) u; close(io)

  write(fname,'(a,a,a,i8.8,a)') trim(input_dir), '/', 'v_',   nstep, '.dat'
open(newunit=io, file=trim(fname), form='unformatted', &
     access='stream', status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) '[error] cannot open ', trim(fname); stop 1
  end if
  read(io) v; close(io)

  write(fname,'(a,a,a,i8.8,a)') trim(input_dir), '/', 'w_',   nstep, '.dat'
open(newunit=io, file=trim(fname), form='unformatted', &
     access='stream', status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) '[error] cannot open ', trim(fname); stop 1
  end if
  read(io) w; close(io)

  write(fname,'(a,a,a,i8.8,a)') trim(input_dir), '/', 'phi_', nstep, '.dat'
open(newunit=io, file=trim(fname), form='unformatted', &
     access='stream', status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) '[error] cannot open ', trim(fname); stop 1
  end if
  read(io) phi; close(io)
end subroutine read_fields_snap


! ----------------------------------------------------------------------
! Average an array of edge-centered derivatives to cell centers, in place.
!
! face_axis  : axis along which the source velocity is face-staggered
!              (1, 2, or 3 for x, y, z)
! deriv_axis : axis along which the derivative was taken
!
! With the minus-half-face convention, the edge derivative E(i,j,k) lives
! at the edge offset by -1/2 in BOTH face_axis and deriv_axis. The four
! edges surrounding cell center (i,j,k) are at index offsets
!   (0,0), (+1,0), (0,+1), (+1,+1)
! along (face_axis, deriv_axis).
! ----------------------------------------------------------------------
subroutine average_edges_to_center(arr, face_axis, deriv_axis)
  use commondata
  implicit none
  double precision, intent(inout) :: arr(nx,nx,nx)
  integer, intent(in) :: face_axis, deriv_axis

  double precision, allocatable :: tmp(:,:,:)
  integer :: i, j, k
  integer :: i_fp, j_fp, k_fp       ! "+1 along face_axis"
  integer :: i_dp, j_dp, k_dp       ! "+1 along deriv_axis"
  integer :: i_bp, j_bp, k_bp       ! "+1 along both"

  allocate(tmp(nx,nx,nx))

  !$omp parallel do default(none) shared(arr,tmp,nx,face_axis,deriv_axis) &
  !$omp private(i,j,k,i_fp,j_fp,k_fp,i_dp,j_dp,k_dp,i_bp,j_bp,k_bp)
  do k = 1, nx
     do j = 1, nx
        do i = 1, nx
           i_fp = i; j_fp = j; k_fp = k
           i_dp = i; j_dp = j; k_dp = k
           select case (face_axis)
              case (1); i_fp = i + 1; if (i_fp > nx) i_fp = 1
              case (2); j_fp = j + 1; if (j_fp > nx) j_fp = 1
              case (3); k_fp = k + 1; if (k_fp > nx) k_fp = 1
           end select
           select case (deriv_axis)
              case (1); i_dp = i + 1; if (i_dp > nx) i_dp = 1
              case (2); j_dp = j + 1; if (j_dp > nx) j_dp = 1
              case (3); k_dp = k + 1; if (k_dp > nx) k_dp = 1
           end select
           i_bp = i_fp; j_bp = j_fp; k_bp = k_fp
           select case (deriv_axis)
              case (1); i_bp = i_bp + 1; if (i_bp > nx) i_bp = 1
              case (2); j_bp = j_bp + 1; if (j_bp > nx) j_bp = 1
              case (3); k_bp = k_bp + 1; if (k_bp > nx) k_bp = 1
           end select
           tmp(i,j,k) = 0.25d0 * ( arr(i,   j,   k  ) &
                                 + arr(i_fp,j_fp,k_fp) &
                                 + arr(i_dp,j_dp,k_dp) &
                                 + arr(i_bp,j_bp,k_bp) )
        end do
     end do
  end do
  !$omp end parallel do

  arr = tmp
  deallocate(tmp)
end subroutine average_edges_to_center


! ----------------------------------------------------------------------
subroutine build_hist(eps, phi, lin_lo, lin_hi, log_lo, log_hi, &
                      hl_in, hl_if, hl_out, hl_tot, &
                      hg_in, hg_if, hg_out, hg_tot, &
                      n_pos_in, n_pos_if, n_pos_out, n_pos_tot)
  use commondata
  use hist_params
  implicit none
  double precision, intent(in) :: eps(nx,nx,nx), phi(nx,nx,nx)
  double precision, intent(in) :: lin_lo, lin_hi, log_lo, log_hi
  integer(kind=8), intent(inout) :: hl_in(nbins_lin), hl_if(nbins_lin), &
                                    hl_out(nbins_lin), hl_tot(nbins_lin)
  integer(kind=8), intent(inout) :: hg_in(nbins_log), hg_if(nbins_log), &
                                    hg_out(nbins_log), hg_tot(nbins_log)
  integer(kind=8), intent(inout) :: n_pos_in, n_pos_if, n_pos_out, n_pos_tot

  integer :: i, j, k, ib_lin, ib_log
  double precision :: e, phi_d, dlin, dlog, le

  dlin = (lin_hi - lin_lo) / dble(nbins_lin)
  dlog = (log_hi - log_lo) / dble(nbins_log)

  do k = 1, nx
     do j = 1, nx
        do i = 1, nx
           e = eps(i,j,k)
           if (phi_drop_is_one) then
              phi_d = phi(i,j,k)
           else
              phi_d = 1.0d0 - phi(i,j,k)
           end if

           ! Linear bin
           ib_lin = int( (e - lin_lo) / dlin ) + 1
           if (ib_lin < 1) ib_lin = 1
           if (ib_lin > nbins_lin) ib_lin = nbins_lin
           hl_tot(ib_lin) = hl_tot(ib_lin) + 1
           if (phi_d >= phi_in_thr) then
              hl_in(ib_lin) = hl_in(ib_lin) + 1
           else if (phi_d <= phi_out_thr) then
              hl_out(ib_lin) = hl_out(ib_lin) + 1
           else
              hl_if(ib_lin) = hl_if(ib_lin) + 1
           end if

           ! Log bin (only positive eps)
           if (e > eps_floor) then
              le = log10(e)
              ib_log = int( (le - log_lo) / dlog ) + 1
              if (ib_log < 1) ib_log = 1
              if (ib_log > nbins_log) ib_log = nbins_log
              hg_tot(ib_log) = hg_tot(ib_log) + 1
              n_pos_tot = n_pos_tot + 1
              if (phi_d >= phi_in_thr) then
                 hg_in(ib_log) = hg_in(ib_log) + 1
                 n_pos_in = n_pos_in + 1
              else if (phi_d <= phi_out_thr) then
                 hg_out(ib_log) = hg_out(ib_log) + 1
                 n_pos_out = n_pos_out + 1
              else
                 hg_if(ib_log) = hg_if(ib_log) + 1
                 n_pos_if = n_pos_if + 1
              end if
           end if
        end do
     end do
  end do
end subroutine build_hist


! ----------------------------------------------------------------------
subroutine write_pdf_file(nstep, tag, hl, hg, n_samples, n_pos_samples, &
                          lin_lo, lin_hi, log_lo, log_hi, eps_mean, eps_sum)
  use commondata
  use hist_params
  implicit none
  integer, intent(in) :: nstep
  character(len=*), intent(in) :: tag
  integer(kind=8), intent(in) :: hl(nbins_lin), hg(nbins_log)
  integer(kind=8), intent(in) :: n_samples, n_pos_samples
  double precision, intent(in) :: lin_lo, lin_hi, log_lo, log_hi
  double precision, intent(in) :: eps_mean, eps_sum

  character(len=256) :: fname
  integer :: io, b, nrows
  double precision :: dlin, dlog, c_lin, c_log, p_lin, p_log

  write(fname,'(a,a,a,i8.8,a)') 'eps_pdf_', trim(tag), '_', nstep, '.dat'
  open(newunit=io, file=trim(fname), status='replace', action='write')

  dlin = (lin_hi - lin_lo) / dble(nbins_lin)
  dlog = (log_hi - log_lo) / dble(nbins_log)

  write(io,'(a)')           '# Conditional dissipation PDF'
  write(io,'(a,a)')         '# subset : ', trim(tag)
  write(io,'(a,i0)')        '# nstep  = ', nstep
  write(io,'(a,i0)')        '# samples       = ', n_samples
  write(io,'(a,i0)')        '# samples (eps>0) = ', n_pos_samples
  write(io,'(a,es16.8)')    '# <eps>        = ', eps_mean
  write(io,'(a,es16.8,a,es16.8)') '# eps bins lin = [', lin_lo, ', ', lin_hi
  write(io,'(a,f10.5,a,f10.5)')   '# log bins     = [', log_lo, ', ', log_hi
  write(io,'(a)') '# col 1: linear bin center (eps)'
  write(io,'(a)') '# col 2: linear PDF       p(eps)'
  write(io,'(a)') '# col 3: log10 bin center log10(eps)'
  write(io,'(a)') '# col 4: log10 PDF        p(log10 eps)'
  write(io,'(a)') '# col 5: count linear bin'
  write(io,'(a)') '# col 6: count log10 bin'

  nrows = max(nbins_lin, nbins_log)
  do b = 1, nrows
     if (b <= nbins_lin) then
        c_lin = lin_lo + (dble(b) - 0.5d0) * dlin
        if (n_samples > 0) then
           p_lin = dble(hl(b)) / (dble(n_samples) * dlin)
        else
           p_lin = 0.0d0
        end if
     else
        c_lin = 0.0d0; p_lin = 0.0d0
     end if
     if (b <= nbins_log) then
        c_log = log_lo + (dble(b) - 0.5d0) * dlog
        if (n_pos_samples > 0) then
           p_log = dble(hg(b)) / (dble(n_pos_samples) * dlog)
        else
           p_log = 0.0d0
        end if
     else
        c_log = 0.0d0; p_log = 0.0d0
     end if
     write(io,'(es16.8,1x,es16.8,1x,es16.8,1x,es16.8,1x,i12,1x,i12)') &
        c_lin, p_lin, c_log, p_log, &
        merge(hl(b),  0_8, b <= nbins_lin), &
        merge(hg(b),  0_8, b <= nbins_log)
  end do

  close(io)
end subroutine write_pdf_file
