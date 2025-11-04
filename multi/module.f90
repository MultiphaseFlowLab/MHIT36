module param
    integer, parameter :: nx=512
    integer :: ny=nx,nz=nx
    double precision :: pi,lx,dx,dxi,ddxi,rhoi,twopi
    integer :: restart,tstart,tfin,dump
    double precision :: gamma, normag
    double precision :: dt,mu,rho !flow parameters
    integer :: inflow, inphi
    double precision :: f1,f2,f3,k0 ! forcing parameters
    double precision :: radius, sigma, epsr, eps, pos, val, epsi, enum ! phase-field parameters
    double precision :: times,timef
    ! other variables (wavenumber, grid location)
    real(8), allocatable :: x(:), kx(:)
    real(8), allocatable :: x_ext(:)
    double precision :: yinf, ysup ! Inf and Sup of y in PiX (no halo)
    double precision :: zinf, zsup ! Inf and Sup of z in PiX (no halo)
    double precision :: lyloc, lzloc ! Inf and Sup of y in PiX
  
end module param


module mpivar
   ! MPI variables
   integer :: rank, ranks, ierr
   integer :: localRank, localComm
   integer :: nidp1y, nidm1y, nidp1z, nidm1z
end module mpivar


module cudecompvar
   use cudecomp
   integer :: npx, npy, npz
   type(cudecompHandle) :: handle
   type(cudecompGridDesc) :: grid_desc,grid_descD2Z
   type(cudecompGridDescConfig) :: config
   type(cudecompGridDescAutotuneOptions) :: options
   integer :: pdims(2) ! pr x pc pencils
   integer :: gdims(3) ! global grid dimensions
   integer :: halo(3) ! halo extensions
   integer :: halo_ext ! 0 no halo, 1 means 1 halo
   type(cudecompPencilInfo) :: piX, piY, piZ ! size of the pencils in x- y- and z-configuration
   type(cudecompPencilInfo) :: piX_d2z, piY_d2z, piZ_d2z  ! size of the pencils in x- y- and z-configuration for D2Z
   type(cudecompPencilInfo) :: piX_Poiss
   integer(8) :: nElemX, nElemY, nElemZ, nElemWork, nElemWork_halo,nElemWork_halo_d2z
   integer(8) :: nElemX_d2z, nElemY_d2z, nElemZ_d2z, nElemWork_d2z
   logical :: halo_periods(3)
   integer :: pix_yoff,pix_zoff
end module cudecompvar


module velocity
   double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   double precision, allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:)
   double precision, allocatable :: rhsu_o(:,:,:), rhsv_o(:,:,:), rhsw_o(:,:,:)
   complex(8), allocatable :: rhsp_complex(:,:,:)
   double precision, allocatable :: rhsp(:,:,:), p(:,:,:)
   double precision, allocatable :: rhspp(:,:,:), pp(:,:,:)
   double precision, allocatable :: div(:,:,:)
   double precision :: uc, vc, wc, umax, gumax=1.0d0, cou, alpha, beta
   double precision :: h11, h12, h13, h21, h22, h23, h31, h32, h33
   double precision :: umean, vmean, wmean, gumean, gvmean, gwmean
   double precision, allocatable :: mysin(:), mycos(:)
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:), psidi(:,:,:), tanh_psi(:,:,:)
   double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
   double precision, allocatable :: normx_f(:,:,:), normy_f(:,:,:), normz_f(:,:,:)
   double precision :: curv
   double precision, allocatable :: fxst(:,:,:), fyst(:,:,:), fzst(:,:,:)
   double precision, allocatable :: rhsphik1(:,:,:),rhsphik2(:,:,:), rhsphik3(:,:,:), rhsphik4(:,:,:), phi_tmp(:,:,:)
   double precision :: normx_xm,normx_xp,normx_ym,normx_yp,normx_zm,normx_zp
   double precision :: normy_xm,normy_xp,normy_ym,normy_yp,normy_zm,normy_zp
   double precision :: normz_xm,normz_xp,normz_ym,normz_yp,normz_zm,normz_zp
   double precision :: rn_01,rn_11,rn_02,rn_12,rn_03,rn_13
   double precision :: sharpxm,sharpxp,sharpym,sharpyp,sharpzm,sharpzp
end module phase




