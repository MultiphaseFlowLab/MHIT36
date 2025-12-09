#define CHECK_CUDECOMP_EXIT(f) if (f /= CUDECOMP_RESULT_SUCCESS) call exit(1)

program main
use cudafor
use cudecomp
use cufft
use mpi
use velocity 
use phase
use temperature
use param
use mpivar
use cudecompvar
use nvtx


implicit none
! timer for scaling test
real :: t_start, t_end, elapsed
! grid dimensions
integer :: comm_backend
integer :: pr, pc
! cudecomp
! cuFFT
integer :: planXf, planXb, planY
integer :: batchsize
integer :: status
integer :: i,j,k,il,jl,kl,ig,jg,kg,t,stage
integer :: im,ip,jm,jp,km,kp,last,idx,kgm
! TDMA variables
double precision, allocatable :: a(:), b(:), c(:)
double complex, allocatable :: d(:), sol(:)
! working arrays
double complex, allocatable :: psi(:)
double precision, allocatable :: ua(:,:,:)
double precision, allocatable :: uaa(:,:,:)
double complex, device, allocatable :: psi_d(:)
double precision, device, allocatable :: vel_d(:) ! only used for implicit diffusion in z
double complex, pointer, device, contiguous :: work_d(:), work_halo_d(:), work_d_d2z(:), work_halo_d_d2z(:)
character(len=40) :: namefile
character(len=4) :: itcount
! Code variables
double precision ::err, maxErr, meanp, gmeanp
double complex, device, pointer :: psi3d(:,:,:)
double precision :: k2,maxdiv
!integer :: il, jl, ig, jg
integer :: offsets(3), xoff, yoff
integer :: np(3)
! Alan Williamson classic
double precision, parameter :: alpha(3) = (/ 8.d0/15.d0,   5.d0/12.d0,   3.d0/4.d0 /)
double precision, parameter :: beta(3)  = (/ 0.d0,       -17.d0/60.d0,  -5.d0/12.d0 /)
! Stage coefficients for diffusion-optimized SSP RK3
!real(kind=8), parameter :: alpha(3) = (/ 0.444370493651235d0, 0.555629506348765d0, 1.0d0 /)
!real(kind=8), parameter :: beta(3)   = (/ 0.0d0, -0.122243120495896d0, -0.377756879504104d0 /)

! Enable or disable phase field 
#define phiflag 0
! Enable or disable temperature field
#define thetaflag 0

!########################################################################################################################################
! 1. INITIALIZATION OF MPI AND cuDECOMP AUTOTUNING : START
!########################################################################################################################################
! MPI initialization, put in rank the local MPI rank number and ranks total number
! Same procedura defined in the cuDecomp documentation
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
call mpi_comm_size(MPI_COMM_WORLD, ranks, ierr)
call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, localComm, ierr)
call mpi_comm_rank(localComm, localRank, ierr)
ierr = cudaSetDevice(localRank) !assign GPU to MPI rank

! Define grid and decomposition
call readinput

! domain decomposition (pencils in y and z) 
pr = 1
pc = 2
halo_ext=1

! CuDECOMP initialization and settings 
comm_backend = 0 ! Enable full autotuning
CHECK_CUDECOMP_EXIT(cudecompInit(handle, MPI_COMM_WORLD))
! config is a struct and pr and pc are the number of pencils along the two directions
! create an uninitialized configuration struct and initialize it to defaults using cudecompGridDescConfigSetDefaults. 
CHECK_CUDECOMP_EXIT(cudecompGridDescConfigSetDefaults(config)) ! Initializing to default values is required to ensure no entries are left uninitialized.
pdims = [pr, pc] !pr and pc are the number of pencil along the different directions
config%pdims = pdims
halo = [0, halo_ext, halo_ext] ! no halo along x neeed because is periodic and in physical space i have x-pencil
! for transpositions
config%transpose_comm_backend = comm_backend
config%transpose_axis_contiguous = .true.
! for halo exchanges
config%halo_comm_backend = CUDECOMP_HALO_COMM_MPI
! Setting for periodic halos in all directions (non required to be in config)
halo_periods = [.true., .true., .false.]
! create spectral grid descriptor first to select pdims for optimal transposes
gdims = [nx/2+1, ny, nz]
config%gdims = gdims
! Set up autotuning options for spectral grid (transpose related settings)
CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_DOUBLE_COMPLEX
if (comm_backend == 0) then
   options%autotune_transpose_backend = .true.
   options%autotune_halo_backend = .false.
endif
options%transpose_use_inplace_buffers = .true.
options%transpose_input_halo_extents(:, 1) = halo
options%transpose_output_halo_extents(:, 4) = halo
CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_descD2Z, config, options))
! create physical grid descriptor
! take previous config and modify the global grid (nx instead of nx/2+1)
! reset transpose_comm_backend to default value to avoid picking up possible nvshmem
! transpose backend selection (this impacts how workspaces are allocated)
gdims = [nx, ny, nz]
config%gdims = gdims
config%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
! Set up autotuning options for physical grid (halo related settings)
CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_DOUBLE_COMPLEX
if (comm_backend == 0) then
   options%autotune_halo_backend = .true.
endif
options%halo_extents(:) = halo
options%halo_periods(:) = halo_periods
options%halo_axis = 1
CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))
! Get pencil info for the grid descriptor in the physical space pencil struct (piX, piY or piZ)
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piX, 1, halo))
nElemX = piX%size !<- number of total elments in x-configuratiion (including halo)
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piY, 2))
nElemY = piY%size
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piZ, 3))
nElemZ = piZ%size
! Get workspace sizes for transpose (1st row, not used) and halo (2nd row, used)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_desc, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_desc, 1, halo, nElemWork_halo))
! Get pencil info for the grid descriptor in the complex space 
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piX_d2z, 1,halo))
nElemX_d2z = piX_d2z%size !<- number of total elments in x-configuratiion (include halo)
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piY_d2z, 2))
nElemY_d2z = piY_d2z%size
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piZ_d2z, 3))
nElemZ_d2z = piZ_d2z%size
! Get workspace sizes for transpose (1st row,used) and halo (2nd row, not used)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_descD2Z, nElemWork_d2z))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_descD2Z, 1, halo, nElemWork_halo_d2z))
! End cuDecomp initialization



! CUFFT initialization -- Create Plans (along x anf y only, z not required)
! Forward 1D FFT in X: D2Z
batchSize = piX_d2z%shape(2)*piX_d2z%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planXf, nx, CUFFT_D2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan Forward'
! Backward 1D FFT in X: Z2D
batchSize = piX_d2z%shape(2)*piX_d2z%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planXb, nx, CUFFT_Z2D, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan Backward'
! it's always 2 and 3 because y-pencil have coordinates y,z,x
batchSize = piY_d2z%shape(2)*piY_d2z%shape(3)
status = cufftPlan1D(planY, ny, CUFFT_Z2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Y plan Forward & Backward'
!########################################################################################################################################
! 1. INITIALIZATION AND cuDECOMP AUTOTUNING : END
!########################################################################################################################################





!########################################################################################################################################
! START STEP 2: ALLOCATE ARRAYS
!########################################################################################################################################
! Pressure variable
allocate(rhsp(piX%shape(1), piX%shape(2), piX%shape(3))) 
allocate(p(piX%shape(1), piX%shape(2), piX%shape(3))) 
allocate(psi_d(max(nElemX_d2z, nElemY_d2z, nElemZ_d2z))) 
!NS variables
allocate(u(piX%shape(1),piX%shape(2),piX%shape(3)),v(piX%shape(1),piX%shape(2),piX%shape(3)),w(piX%shape(1),piX%shape(2),piX%shape(3))) !velocity vector
allocate(rhsu(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(rhsu_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw_o(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(div(piX%shape(1),piX%shape(2),piX%shape(3))) 
!TDMA solver
allocate(a(0:nz+1),b(0:nz+1),c(0:nz+1),d(0:nz+1),sol(0:nz+1))
!PFM variables
#if phiflag == 1
allocate(phi(piX%shape(1),piX%shape(2),piX%shape(3)),rhsphi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(psidi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(tanh_psi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(normx(piX%shape(1),piX%shape(2),piX%shape(3)),normy(piX%shape(1),piX%shape(2),piX%shape(3)),normz(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(fxst(piX%shape(1),piX%shape(2),piX%shape(3)),fyst(piX%shape(1),piX%shape(2),piX%shape(3)),fzst(piX%shape(1),piX%shape(2),piX%shape(3))) ! surface tension forces
#endif
!Temperature variables
#if thetaflag == 1
allocate(theta(piX%shape(1),piX%shape(2),piX%shape(3)),rhstheta(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(rhstheta_o(piX%shape(1),piX%shape(2),piX%shape(3)))
#endif
! allocate arrays for transpositions and halo exchanges 
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_d, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_halo_d, nElemWork_halo))
! allocate arrays for transpositions
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_descD2Z, work_d_d2z, nElemWork_d2z))
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_descD2Z, work_halo_d_d2z, nElemWork_halo_d2z)) ! not required 
!########################################################################################################################################
! END STEP2: ALLOCATE ARRAYS
!########################################################################################################################################








!########################################################################################################################################
! START STEP 3: FLOW FIELD, PHASE-FIELD AND TEMPERATURE INIT
!########################################################################################################################################
! 3.1 Read/initialize from data without halo grid points (avoid out-of-bound if reading usin MPI I/O)
! 3.2 Call halo exchnages along Y and Z for u, v, w, phi and theta
if (restart .eq. 0) then !fresh start Taylor Green or read from file in init folder
if (rank.eq.0) write(*,*) "Initialize velocity field (fresh start)"
   if (inflow .eq. 0) then
      if (rank.eq.0) write(*,*) "Initialize zero velocity field"
      do k = 1+halo_ext, piX%shape(3)-halo_ext
         do j = 1+halo_ext, piX%shape(2)-halo_ext
            do i = 1, piX%shape(1)
               u(i,j,k) =  0.0d0
               v(i,j,k) =  0.0d0 
               w(i,j,k) =  0.0d0 
            enddo
         enddo
      enddo
   endif
   if (inflow .eq. 1) then
      if (rank.eq.0) write(*,*) "Initialize laminar flow (x) + 3D perturbation"
      do k = 1+halo_ext, piX%shape(3)-halo_ext
         kg = piX%lo(3) + k - 1 - halo_ext                   
         do j = 1+halo_ext, piX%shape(2)-halo_ext
            jg = piX%lo(2) + j - 1 - halo_ext
            do i = 1, piX%shape(1)
               amp=3.d0
               mx= 3.03d0
               my= 2.02d0
               mz= 4.00d0
               !3D divergence free flow with fluctuations that satisfies the boundary conditions
               u(i,j,k) =  20.d0*(1.d0 - ((2*z(kg) - lz)/lz)**2) !
               u(i,j,k) =  u(i,j,k) - amp*cos(twopi*mx*x(i)/lx)*sin(twopi*my*y(jg)/ly)*2.d0*twopi/lz*sin(twopi*z(kg)/lz)*cos(twopi*z(kg)/lz)
               u(i,j,k) =  u(i,j,k) + amp*sin(twopi*mx*x(i)/lx)*(-twopi*my/ly)*sin(2.d0*twopi*my*y(jg)/ly)*sin(twopi*z(kg)/lz)*sin(twopi*z(kg)/lz)
               v(i,j,k) = -amp*cos(twopi*my*y(jg)/ly)*(twopi*mx/lx)*cos(twopi*mx*x(i)/lx)*sin(twopi*z(kg)/lz)*sin(twopi*z(kg)/lz)
               w(i,j,k) =  amp*cos(twopi*mx*x(i)/lx)*(twopi*mx/lx)*sin(twopi*my*y(jg)/ly)*sin(twopi*z(kg)/lz)*sin(twopi*z(kg)/lz)
            enddo
         enddo
      enddo
   endif
   if (inflow .eq. 2) then
      if (rank.eq.0)  write(*,*) "Initialize from data"
      call readfield(1)
      call readfield(2)
      call readfield(3)
   endif
endif
if (restart .eq. 1) then !restart, ignore inflow and read the tstart field 
   if (rank.eq.0)  write(*,*) "Initialize velocity field (from output folder), iteration:", tstart
   call readfield_restart(tstart,1)
   call readfield_restart(tstart,2)
   call readfield_restart(tstart,3)
endif

! update halo cells along y and z directions (enough only if pr and pc are non-unitary)
!$acc host_data use_device(u)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
!$acc host_data use_device(v)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
!$acc host_data use_device(w)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 

! initialize phase-field
#if phiflag == 1
if (restart .eq. 0) then
if (rank.eq.0) write(*,*) 'Initialize phase field (fresh start)'
   if (inphi .eq. 0) then
   if (rank.eq.0) write(*,*) 'Spherical drop'
      do k = 1+halo_ext, piX%shape(3)-halo_ext
      kg = piX%lo(3) + k - 1 - halo_ext
         do j = 1+halo_ext, piX%shape(2)-halo_ext
         jg = piX%lo(2) + j - 1 - halo_ext
            do i = 1, piX%shape(1)
                pos=(x(i)-lx/2)**2d0 +  (y(jg)-ly/2)**2d0 + (z(kg)-lz/2)**2d0
                phi(i,j,k) = 0.5d0*(1.d0-tanh((sqrt(pos)-radius)/2/eps))
            enddo
        enddo
    enddo
   endif
   if (inphi .eq. 1) then
      if (rank.eq.0)  write(*,*) "Initialize phase-field from data"
      call readfield(5)
   endif
endif
if (restart .eq. 1) then
    write(*,*) "Initialize phase-field (restart, from output folder), iteration:", tstart
    call readfield_restart(tstart,5)
endif
! update halo
!$acc host_data use_device(phi)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
#endif

! initialize temperature field
#if thetaflag == 1
if (restart .eq. 0) then
if (rank.eq.0) write(*,*) 'Initialize temperature field (fresh start)'
   if (intheta .eq. 0) then
   if (rank.eq.0) write(*,*) 'Uniform temperature field'
      do k = 1+halo_ext, piX%shape(3)-halo_ext
         kg = piX%lo(3) + k - 1 - halo_ext
         do j = 1+halo_ext, piX%shape(2)-halo_ext
            do i = 1, piX%shape(1)
               theta(i,j,k) = 0.d0!1.d0 - z(kg)  ! uniform temperature
            enddo
         enddo
      enddo
   endif
   if (intheta .eq. 1) then
      if (rank.eq.0) write(*,*) "Initialize temperature from data"
      call readfield(6)
   endif
endif
if (restart .eq. 1) then
    write(*,*) "Initialize temperature (restart, from output folder), iteration:", tstart
    call readfield_restart(tstart,6)
endif
! update halo
!$acc host_data use_device(theta)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, theta, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, theta, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
#endif

!Save initial fields (only if a fresh start)
if (restart .eq. 0) then
   if (rank.eq.0) write(*,*) "Save initial fields"
   call writefield(tstart,1)
   call writefield(tstart,2)
   call writefield(tstart,3)
   call writefield(tstart,4)
   #if phiflag == 1
   call writefield(tstart,5)
   #endif
   #if thetaflag == 1
   call writefield(tstart,6)  ! temperature
   #endif
endif
!########################################################################################################################################
! END STEP 3: FLOW FIELD, PHASE-FIELD AND TEMP INIT FIELD INIT
!########################################################################################################################################






! ########################################################################################################################################
! START TEMPORAL LOOP: STEP 4 to 9 REPEATED AT EVERY TIME STEP
! ########################################################################################################################################
! First step use Euler
gumax=1.d0
tstart=tstart+1
gamma=1.d0*gumax
!$acc data copyin(piX)
#if thetaflag == 1
!$acc data create(rhsu_o, rhsv_o, rhsw_o, rhstheta_o)
#else
!$acc data create(rhsu_o, rhsv_o, rhsw_o)
#endif
!$acc data copyin(mysin, mycos)
call cpu_time(t_start)
! Start temporal loop
do t=tstart,tfin
    ! Create custom label for each marker
    write(itcount,'(i4)') t
    ! Range with custom  color
    call nvtxStartRange("Iteration "//itcount,t)

    if (rank.eq.0) write(*,*) "Time step",t,"of",tfin
    call cpu_time(times)

   call nvtxStartRange("Phase-field")
   !########################################################################################################################################
   ! START STEP 4: PHASE-FIELD SOLVER (EXPLICIT)
   !########################################################################################################################################
   #if phiflag == 1
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            ! compute distance function psi (used to compute normals)
            val = min(phi(i,j,k),1.0d0) ! avoid machine precision overshoots in phi that leads to problem with log
            psidi(i,j,k) = eps*log((val+enum)/(1.d0-val+enum))
            ! compute here the tanh of distance function psi (used in the sharpening term) to avoid multiple computations of tanh
            tanh_psi(i,j,k) = tanh(0.5d0*psidi(i,j,k)*epsi)
         enddo
      enddo
   enddo
   !$acc end kernels

   gamma=1.d0*gumax
   !$acc parallel loop tile(16,4,2)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ! 4.1 RHS computation
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            kg=piX%lo(3) + k - 1 - halo_ext
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            ! convective (first three lines) and diffusive (last three lines)
            rhsphi(i,j,k) =   &
                  - (u(ip,j,k)*0.5d0*(phi(ip,j,k)+phi(i,j,k)) - u(i,j,k)*0.5d0*(phi(i,j,k)+phi(im,j,k)))*dxi   &  
                  - (v(i,jp,k)*0.5d0*(phi(i,jp,k)+phi(i,j,k)) - v(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,jm,k)))*dyi   &  
                  - (w(i,j,kp)*0.5d0*(phi(i,j,kp)+phi(i,j,k)) - w(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,j,km)))*dzci(kg)  &  
                        + gamma*(eps*(phi(ip,j,k)-2.d0*phi(i,j,k)+phi(im,j,k))*ddxi + &                   
                                 eps*(phi(i,jp,k)-2.d0*phi(i,j,k)+phi(i,jm,k))*ddyi + &                   
                                 eps*((phi(i,j,kp)-phi(i,j,k))*dzi(kg+1) - (phi(i,j,k) -phi(i,j,km))*dzi(kg))*dzci(kg))     ! first between centers and then betwenn faces                
            ! Compute normals for sharpening term (gradient)
            normx(i,j,k) = 0.5d0*(psidi(ip,j,k) - psidi(im,j,k))*dxi
            normy(i,j,k) = 0.5d0*(psidi(i,jp,k) - psidi(i,jm,k))*dyi
            normz(i,j,k) =       (psidi(i,j,kp) - psidi(i,j,km))/(z(kg+1)-z(kg-1))
         enddo
      enddo
   enddo

   ! Update normx,normy and normz halos, required to then compute normal derivative
   !$acc host_data use_device(normx,normy,normz)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 

   ! 4.1.3. Compute Sharpening term (gradient)
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            normod = 1.d0/(sqrt(normx(i,j,k)*normx(i,j,k) + normy(i,j,k)*normy(i,j,k) + normz(i,j,k)*normz(i,j,k)) + 1.0E-16)
            ! normod = 1.d0/(sqrt(normx(i,j,k)**2d0 + normy(i,j,k)**2d0 + normz(i,j,k)**2d0) + 1.0E-16)
            normx(i,j,k) = normx(i,j,k)*normod
            normy(i,j,k) = normy(i,j,k)*normod
            normz(i,j,k) = normz(i,j,k)*normod
         enddo
      enddo
   enddo
   !$acc end kernels

   ! Compute sharpening term
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            kg=piX%lo(3) + k - 1 - halo_ext
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            rhsphi(i,j,k)=rhsphi(i,j,k)-gamma*((0.25d0*(1.d0-tanh_psi(ip,j,k)*tanh_psi(ip,j,k))*normx(ip,j,k) - &
                                                0.25d0*(1.d0-tanh_psi(im,j,k)*tanh_psi(im,j,k))*normx(im,j,k))*0.5*dxi + &
                                               (0.25d0*(1.d0-tanh_psi(i,jp,k)*tanh_psi(i,jp,k))*normy(i,jp,k) - &
                                                0.25d0*(1.d0-tanh_psi(i,jm,k)*tanh_psi(i,jm,k))*normy(i,jm,k))*0.5*dyi + &
                                               (0.25d0*(1.d0-tanh_psi(i,j,kp)*tanh_psi(i,j,kp))*normz(i,j,kp) - &
                                                0.25d0*(1.d0-tanh_psi(i,j,km)*tanh_psi(i,j,km))*normz(i,j,km))/(z(kg+1)-z(kg-1))) 
         enddo
      enddo
   enddo
    !$acc end kernels

   ! Get phi at n+1 
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            phi(i,j,k) = phi(i,j,k) + dt*rhsphi(i,j,k)
         enddo
      enddo
   enddo
   !$acc end kernels

   !$acc host_data use_device(phi)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   #endif
   !########################################################################################################################################
   ! END STEP 4: PHASE-FIELD SOLVER (EXPLICIT)
   !########################################################################################################################################
   call nvtxEndRange


   !########################################################################################################################################
   ! START STEP 5: ENERGY EQUATION SOLVER (RK3 EXPLICIT)
   !########################################################################################################################################
   #if thetaflag == 1
   ! Temperature solver (RK3 explicit)
   do stage = 1,3
      !$acc parallel loop tile(16,4,2)
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               kg=piX%lo(3) + k - 1 - halo_ext
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               ! convective terms
               rhstheta(i,j,k) = &
                     - (u(ip,j,k)*0.5d0*(theta(ip,j,k)+theta(i,j,k)) - u(i,j,k)*0.5d0*(theta(i,j,k)+theta(im,j,k)))*dxi &
                     - (v(i,jp,k)*0.5d0*(theta(i,jp,k)+theta(i,j,k)) - v(i,j,k)*0.5d0*(theta(i,j,k)+theta(i,jm,k)))*dyi &
                     - (w(i,j,kp)*0.5d0*(theta(i,j,kp)+theta(i,j,k)) - w(i,j,k)*0.5d0*(theta(i,j,k)+theta(i,j,km)))*dzci(kg)
               ! diffusive terms
               rhstheta(i,j,k) = rhstheta(i,j,k) + kappa*((theta(ip,j,k)-2.d0*theta(i,j,k)+theta(im,j,k))*ddxi + &
                                                          (theta(i,jp,k)-2.d0*theta(i,j,k)+theta(i,jm,k))*ddyi + &
                                                          ((theta(i,j,kp)-theta(i,j,k))*dzi(kg+1) - (theta(i,j,k) -theta(i,j,km))*dzi(kg))*dzci(kg))  

               enddo
         enddo
      enddo
      !$acc parallel loop collapse(3)
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               theta(i,j,k) = theta(i,j,k) + dt*alpha(stage)*rhstheta(i,j,k) + dt*beta(stage)*rhstheta_o(i,j,k)
               rhstheta_o(i,j,k)=rhstheta(i,j,k)
            enddo
         enddo
      enddo
      ! update halos 
      !$acc host_data use_device(theta)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, theta, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, theta, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
      ! impose temperature boundary conditions, interpolated as done for the velocity (see node sketch)
      !$acc parallel loop collapse(3)
      do k=1, piX%shape(3)
         do j=1, piX%shape(2)
            do i=1,nx
               kg = piX%lo(3) + k - 1 - halo_ext                   
               if (kg .eq. 1)    theta(i,j,k-1) =  2.d0*( 1.d0) - theta(i,j,k)     ! mean value between kg and kg-1 (top wall) equal to 1 
               if (kg .eq. nz)   theta(i,j,k+1) =  2.d0*(-1.d0) - theta(i,j,k)     ! mean value between kg and kg+1 (bottom wall) equal to -1 
            enddo
         enddo
      enddo
   enddo
   #endif
   !########################################################################################################################################
   ! END STEP 5: ENERGY EQUATION SOLVER (RK3 EXPLICIT)
   !########################################################################################################################################




   call nvtxStartRange("Projection")
   !########################################################################################################################################
   ! START STEP 6: USTAR COMPUTATION (PROJECTION STEP)
   !########################################################################################################################################
   ! 6.1 compute rhs 
   ! 6.2 obtain ustar and store old rhs in rhs_o
   ! 6.3 Call halo exchnages along Y and Z for u,v,w

   ! Projection step
   do stage = 1,3
      !$acc parallel loop tile(16,4,2) 
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               kg=piX%lo(3) + k - 1 - halo_ext
               if (ip .gt. nx) ip=1  
               if (im .lt. 1) im=nx
               !  compute the products (conservative form)
               h11 = 0.25d0*((u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k)))*dxi
               h12 = 0.25d0*((u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k)))*dyi
               h13 = 0.25d0*((u(i,j,kp)+u(i,j,k))*(w(i,j,kp)+w(im,j,kp))   - (u(i,j,k)+u(i,j,km))*(w(i,j,k)+w(im,j,k)))*dzci(kg) ! divide by cell height
               h21 = 0.25d0*((u(ip,j,k)+u(ip,jm,k))*(v(ip,j,k)+v(i,j,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k)))*dxi
               h22 = 0.25d0*((v(i,jp,k)+v(i,j,k))*(v(i,jp,k)+v(i,j,k))     - (v(i,j,k)+v(i,jm,k))*(v(i,j,k)+v(i,jm,k)))*dyi
               h23 = 0.25d0*((w(i,j,kp)+w(i,jm,kp))*(v(i,j,kp)+v(i,j,k))   - (w(i,j,k)+w(i,jm,k))*(v(i,j,k)+v(i,j,km)))*dzci(kg) ! divide by cell height
               h31 = 0.25d0*((w(ip,j,k)+w(i,j,k))*(u(ip,j,k)+u(ip,j,km))   - (w(i,j,k)+w(im,j,k))*(u(i,j,k)+u(i,j,km)))*dxi
               h32 = 0.25d0*((v(i,jp,k)+v(i,jp,km))*(w(i,jp,k)+w(i,j,k))   - (v(i,j,k)+v(i,j,km))*(w(i,j,k)+w(i,jm,k)))*dyi
               h33 = 0.25d0*((w(i,j,kp)+w(i,j,k))*(w(i,j,kp)+w(i,j,k))     - (w(i,j,k)+w(i,j,km))*(w(i,j,k)+w(i,j,km)))*dzi(kg) ! divie by distance between centers
               rhsu(i,j,k)=-(h11+h12+h13)
               rhsv(i,j,k)=-(h21+h22+h23)
               rhsw(i,j,k)=-(h31+h32+h33)
               ! viscous/diffusive terms
               h11 = mu*(u(ip,j,k)-2.d0*u(i,j,k)+u(im,j,k))*ddxi
               h12 = mu*(u(i,jp,k)-2.d0*u(i,j,k)+u(i,jm,k))*ddyi
               h13 = mu*((u(i,j,kp)-u(i,j,k))*dzi(kg+1)-(u(i,j,k)-u(i,j,km))*dzi(kg))*dzci(kg)
               h21 = mu*(v(ip,j,k)-2.d0*v(i,j,k)+v(im,j,k))*ddxi
               h22 = mu*(v(i,jp,k)-2.d0*v(i,j,k)+v(i,jm,k))*ddyi
               h23 = mu*((v(i,j,kp)-v(i,j,k))*dzi(kg+1)-(v(i,j,k)-v(i,j,km))*dzi(kg))*dzci(kg)
               h31 = mu*(w(ip,j,k)-2.d0*w(i,j,k)+w(im,j,k))*ddxi
               h32 = mu*(w(i,jp,k)-2.d0*w(i,j,k)+w(i,jm,k))*ddyi
               h33 = mu*((w(i,j,kp)-w(i,j,k))*dzci(kg+1)-(w(i,j,k)-w(i,j,km))*dzci(kg))*dzi(kg) ! face to face and then center to center
               rhsu(i,j,k)=rhsu(i,j,k)+(h11+h12+h13)*rhoi
               rhsv(i,j,k)=rhsv(i,j,k)+(h21+h22+h23)*rhoi
               rhsw(i,j,k)=rhsw(i,j,k)+(h31+h32+h33)*rhoi               
               ! Pressure driven
               rhsu(i,j,k)=rhsu(i,j,k) - gradpx
               rhsv(i,j,k)=rhsv(i,j,k) - gradpy
            enddo
         enddo
      enddo

      ! Surface tension forces (Obtain surface tension forces evaluated at the center of the cell (same as where phi is located)
      #if phiflag == 1
      !$acc kernels
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               kg=piX%lo(3) + k - 1 - halo_ext
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               curv=0.5d0*(normx(ip,j,k)-normx(im,j,k))*dxi + 0.5d0*(normy(i,jp,k)-normy(i,jm,k))*dyi + (normz(i,j,kp)-normz(i,j,km))/(z(kg+1)-z(kg-1))
               fxst(i,j,k)= -sigma*curv*0.5d0*(phi(ip,j,k)-phi(im,j,k))*dxi
               fyst(i,j,k)= -sigma*curv*0.5d0*(phi(i,jp,k)-phi(i,jm,k))*dyi
               fzst(i,j,k)= -sigma*curv*(phi(i,j,kp)-phi(i,j,km))/(z(kg+1)-z(kg-1))
            enddo
         enddo
      enddo
      !$acc end kernels

      ! Update halo of fxst, fyst and fzst (required then to interpolate at velocity points)
      !$acc host_data use_device(fxst,fyst,fzst)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 
   
      ! Interpolate force at velocity points
      !$acc parallel loop collapse(3)
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               im=i-1
               jm=j-1
               km=k-1
               if (im .lt. 1) im=nx
               #if thetaflag == 0 
               ! Add surface tension force only (ACDI only)
               rhsu(i,j,k)=rhsu(i,j,k) + 0.5d0*(fxst(im,j,k)+fxst(i,j,k))*rhoi
               rhsv(i,j,k)=rhsv(i,j,k) + 0.5d0*(fyst(i,jm,k)+fyst(i,j,k))*rhoi
               rhsw(i,j,k)=rhsw(i,j,k) + 0.5d0*(fzst(i,j,km)+fzst(i,j,k))*rhoi
               #elif thetaflag == 1
               ! Add here also buoyancy force if temperature is active (case ACDI + temperature)
               rhsw(i,j,k)=rhsw(i,j,k) + alphag*0.5d0*(theta(i,j,km)+theta(i,j,k))
               #endif
               u(i,j,k) = u(i,j,k) + dt*alpha(stage)*rhsu(i,j,k) + dt*beta(stage)*rhsu_o(i,j,k)! -dt*(alpha(stage)+beta(stage))*rho*(p(i,j,k)-p(im,j,k))*dxi
               v(i,j,k) = v(i,j,k) + dt*alpha(stage)*rhsv(i,j,k) + dt*beta(stage)*rhsv_o(i,j,k)! -dt*(alpha(stage)+beta(stage))*rho*(p(i,j,k)-p(i,jm,k))*dyi
               w(i,j,k) = w(i,j,k) + dt*alpha(stage)*rhsw(i,j,k) + dt*beta(stage)*rhsw_o(i,j,k)! -dt*(alpha(stage)+beta(stage))*rho*(p(i,j,k)-p(i,j,km))*dzi
               rhsu_o(i,j,k)=rhsu(i,j,k)
               rhsv_o(i,j,k)=rhsv(i,j,k)
               rhsw_o(i,j,k)=rhsw(i,j,k)
            enddo
         enddo
      enddo
      #else
      ! 5.2 find u, v and w star (RK3), single-phase case
      !$acc parallel loop collapse(3)
      do k=1+halo_ext, piX%shape(3)-halo_ext
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ! Add here also buoyancy force if temperature is active (case NS + temperature)
               #if thetaflag == 1
               km=k-1
               rhsw(i,j,k)=rhsw(i,j,k) + alphag*0.5d0*(theta(i,j,km)+theta(i,j,k))
               #endif
               u(i,j,k) = u(i,j,k) + dt*alpha(stage)*rhsu(i,j,k) + dt*beta(stage)*rhsu_o(i,j,k)
               v(i,j,k) = v(i,j,k) + dt*alpha(stage)*rhsv(i,j,k) + dt*beta(stage)*rhsv_o(i,j,k)
               w(i,j,k) = w(i,j,k) + dt*alpha(stage)*rhsw(i,j,k) + dt*beta(stage)*rhsw_o(i,j,k)
               rhsu_o(i,j,k)=rhsu(i,j,k)
               rhsv_o(i,j,k)=rhsv(i,j,k)
               rhsw_o(i,j,k)=rhsw(i,j,k)
            enddo
         enddo
      enddo
      #endif

      ! 8.3 update halos (y direction), required to then compute the RHS of Poisson equation because of staggered grid
      !$acc host_data use_device(u,v,w)
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
      CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
      !$acc end host_data 

      ! impose velocity boundary conditions, w is at the wall, u and v interpolate so that the mean value is zero, no-slip assumted, i.e. u=0, can be extented to any value
      !$acc parallel loop collapse(3) 
      do k=1, piX%shape(3)
         do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               kg = piX%lo(3) + k - 1 -halo_ext 
               ! bottom wall 
               if (kg .eq. 0)     u(i,j,k)=  -u(i,j,k+1)  !  mean value between kg and kg-1 (wall) equal to zero  
               if (kg .eq. 0)     v(i,j,k)=  -v(i,j,k+1)  !  mean value between kg and kg-1 (wall) equal to zero  
               if (kg .eq. 1)     w(i,j,k)=   0.d0           ! w point is at the wall
               ! top wall
               if (kg .eq. nz+1)   u(i,j,k)=  -u(i,j,k-1)  !  mean value between kg and kg+1 (wall) equal to zero 
               if (kg .eq. nz+1)   v(i,j,k)=  -v(i,j,k-1)  !  mean value between kg and kg+1 (wall) equal to zero 
               if (kg .eq. nz+1)   w(i,j,k)=   0.d0           ! w point (nz+1) is at the wall   
            enddo
         enddo
      enddo  

   enddo
   call nvtxEndRange


   !########################################################################################################################################
   ! END STEP 6: USTAR COMPUTATION 
   !########################################################################################################################################
   call nvtxEndRange


   call nvtxStartRange("Poisson")
   !########################################################################################################################################
   ! START STEP 7: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################
   ! initialize rhs and analytical solution
   ! 7.1 Compute rhs of Poisson equation div*ustar: divergence at the cell center 
   ! I've done the halo updates so to compute the divergence at the pencil border i have the *star from the halo
   call nvtxStartRange("compute RHS")
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            kg = piX%lo(3)  + k - 1 - halo_ext
            if (ip > nx) ip=1
            rhsp(i,j,k) =                    (rho*dxi/dt)*(u(ip,j,k)-u(i,j,k))
            rhsp(i,j,k) = rhsp(i,j,k) +      (rho*dyi/dt)*(v(i,jp,k)-v(i,j,k))
            rhsp(i,j,k) = rhsp(i,j,k) + (rho*dzci(kg)/dt)*(w(i,j,kp)-w(i,j,k))
         enddo
      enddo
   enddo
   !$acc end kernels
   call nvtxEndRange


   call nvtxStartRange("FFT forward w/ transpositions")
   !$acc host_data use_device(rhsp)
   status = cufftExecD2Z(planXf, rhsp, psi_d)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X forward error: ', status
   !$acc end host_data
   ! psi(kx,y,z) -> psi(y,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeXToY(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX,piX_d2z%halo_extents, [0,0,0]))
   ! psi(y,z,kx) -> psi(ky,z,kx)
   status = cufftExecZ2Z(planY, psi_d, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y forward error: ', status
   ! psi(ky,z,kx) -> psi(z,kx,ky)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToZ(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX)) 

   call nvtxEndRange
   np(piZ_d2z%order(1)) = piZ_d2z%shape(1)
   np(piZ_d2z%order(2)) = piZ_d2z%shape(2)
   np(piZ_d2z%order(3)) = piZ_d2z%shape(3)
   call c_f_pointer(c_devloc(psi_d), psi3d, piZ_d2z%shape)
   offsets(piZ_d2z%order(1)) = piZ_d2z%lo(1) - 1
   offsets(piZ_d2z%order(2)) = piZ_d2z%lo(2) - 1
   offsets(piZ_d2z%order(3)) = piZ_d2z%lo(3) - 1

   xoff = offsets(1)
   yoff = offsets(2)
   npx = np(1)
   npy = np(2)
   call nvtxStartRange("Solution")
   !$acc parallel loop collapse(2) gang private(a,b,c,d,factor) 
   do jl = 1, npy
      do il = 1, npx
         ! compute global index ig and jg
         jg = yoff + jl
         ig = xoff + il
         ! Set up tridiagonal system for each i and j
         ! Fill diagonals and rhs for each
         ! 0 and nz+1 are the ghost nodes
         do k = 1, nz
            a(k) =  2.0d0*(dzi(k)**2.d0*dzi(k+1))/(dzi(k)+dzi(k+1))
            c(k) =  2.0d0*(dzi(k)*dzi(k+1)**2.d0)/(dzi(k)+dzi(k+1))
            b(k) =  -a(k) - c(k) + 2.d0*(cos(kx_d(ig)*dx)-1.d0)*dxi*dxi + 2.d0*(cos(ky_d(jg)*dy)-1.d0)*dyi*dyi !opt: precompute cosines? 
            d(k) =  psi3d(k,il,jl)
         enddo
         ! Neumann BC at bottom
         a(0) =  0.d0
         b(0) = -1.d0
         c(0) =  1.d0  
         d(0) =  0.d0
         ! Neumann BC at top
         a(nz+1) =  1.d0
         b(nz+1) = -1.d0
         c(nz+1) =  0.d0
         d(nz+1) =  0.d0
         ! Enforce pressure at one point? one interior point, avodig messing up with BC
         if (ig == 1 .and. jg == 1) then
            a(nz+1) = 0.d0
            b(nz+1) = 1.d0
            c(nz+1) = 0.d0
         end if
         ! Forward elimination (Thomas)
         !$acc loop seq
         do k = 1, nz+1
            factor = a(k)/b(k-1)
            b(k) = b(k) - factor*c(k-1)
            d(k) = d(k) - factor*d(k-1)
         end do
         ! Back substitution
         psi3d(nz,il,jl) = (d(nz) - c(nz)*d(nz+1)/b(nz+1))/b(nz)
         !$acc loop seq
         do k = nz-1, 1, -1
            psi3d(k,il,jl) = (d(k) - c(k)*psi3d(k+1,il,jl))/b(k)
         end do
      end do
   end do

   call nvtxStartRange("FFT backwards along x and y w/ transpositions")
   ! psi(z,kx,ky) -> psi(ky,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeZToY(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX))
   ! psi(ky,z,kx) -> psi(y,z,kx)
   status = cufftExecZ2Z(planY, psi_d, psi_d, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y inverse error: ', status
   ! psi(y,z,kx) -> psi(kx,y,z)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToX(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX,[0,0,0], piX_d2z%halo_extents))
   !$acc host_data use_device(p)
   ! psi(kx,y,z) -> p(x,y,z)
   status = cufftExecZ2D(planXb, psi_d, p)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X inverse error: ', status
   !$acc end host_data

   ! normalize pressure (must be done here, not in the TDMA)
   !$acc parallel loop collapse(3)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            p(i,j,k) = p(i,j,k)/dble(nx*ny)
         end do
      end do
   end do
   ! update halo nodes with pressure 
   !$acc host_data use_device(p)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, p, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, p, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !########################################################################################################################################
   ! END STEP 7: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################


   call nvtxStartRange("Correction")
   !########################################################################################################################################
   ! START STEP 8: VELOCITY CORRECTION
   ! ########################################################################################################################################
   ! 8.1 Correct velocity 
   ! 8.2 Call halo update
   ! Correct velocity, pressure has also the halo
   umax=0.d0
   vmax=0.d0
   wmax=0.d0
   !$acc parallel loop collapse(3) reduction(max:umax,vmax,wmax)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i = 1, piX%shape(1) ! equal to nx (no halo on x)
              im=i-1
              jm=j-1
              km=k-1
              kg=piX%lo(3)  + k - 1 - halo_ext
              if (im < 1) im=nx
              if (kg .eq. 1) then
              u(i,j,k)=u(i,j,k) - dt/rho*(p(i,j,k)-p(im,j,k))*dxi
              v(i,j,k)=v(i,j,k) - dt/rho*(p(i,j,k)-p(i,jm,k))*dyi
              else
              u(i,j,k)=u(i,j,k) - dt/rho*(p(i,j,k)-p(im,j,k))*dxi
              v(i,j,k)=v(i,j,k) - dt/rho*(p(i,j,k)-p(i,jm,k))*dyi
              w(i,j,k)=w(i,j,k) - dt/rho*(p(i,j,k)-p(i,j,km))*dzi(kg)
              endif
              !umax=max(umax,u(i,j,k))
              !vmax=max(vmax,v(i,j,k))
              !wmax=max(wmax,w(i,j,k))
          enddo
      enddo
   enddo

   ! 8.3 update halos (y direction), required to then compute the RHS of Poisson equation because of staggered grid
   !$acc host_data use_device(u,v,w)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 

   maxdiv=0.d0
   !$acc parallel loop collapse(3) reduction(max:maxdiv)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i = 1, piX%shape(1) ! equal to nx (no halo on x)
              ip=i+1
              jp=j+1
              kp=k+1
              kg=piX%lo(3)  + k - 1 - halo_ext
              if (ip > nx) ip=1
              div(i,j,k)=dxi*(u(ip,j,k)-u(i,j,k)) + dyi*(v(i,jp,k)-v(i,j,k)) + dzci(kg)*(w(i,j,kp)-w(i,j,k))
              maxdiv=max(maxdiv,abs(div(i,j,k)))
          enddo
      enddo
   enddo
   write(*,*) "Max divergence after correction ", maxdiv

   call MPI_Allreduce(umax,gumax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(vmax,gvmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(wmax,gwmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)

   cflx=gumax*dt*dxi
   cfly=gvmax*dt*dyi
   cflz=gwmax*dt*lz/nz
   cou=max(cflx,cfly)
   cou=max(cou,gcflz)
   if (rank.eq.0) then
      write(*,*) "CFL (max among tasks)", cou
      if (cou .gt. 7) stop
   endif
   
   call cpu_time(timef)
   if (rank.eq.0) print '(" Time elapsed = ",f6.1," ms")',1000*(timef-times)
   !########################################################################################################################################
   ! END STEP 8: VELOCITY CORRECTION  
   !########################################################################################################################################
   call nvtxEndRange


   !########################################################################################################################################
   ! START STEP 9: OUTPUT FIELDS 
   ! ########################################################################################################################################
   if (mod(t,dump) .eq. 0) then
      if (rank .eq. 0) write(*,*) "Saving output files"
         ! write velocity and pressure fiels (1-4)
         call writefield(t,1)
         call writefield(t,2)
         call writefield(t,3)
         call writefield(t,4)
         #if phiflag == 1
         ! write phase-field (5)
         call writefield(t,5)
         #endif
         #if thetaflag == 1
         ! write temperature field (6)
         call writefield(t,6)
         #endif
   endif
   !########################################################################################################################################
   ! END STEP 9: OUTPUT FIELDS N  
   !########################################################################################################################################

call nvtxEndRange
!call nvtxEndRange
enddo
call cpu_time(t_end)
elapsed = t_end-t_start
if (rank .eq. 0) write(*,*)  'Elapsed time (seconds):', elapsed
!$acc end data
!$acc end data
!$acc end data

! Remove allocated variables (add new)
deallocate(u,v,w)
deallocate(tanh_psi, mysin, mycos)
deallocate(rhsu,rhsv,rhsw)
deallocate(rhsu_o,rhsv_o,rhsw_o)
#if phiflag == 1
deallocate(phi,rhsphi,normx,normy,normz)
#endif
#if thetaflag == 1
deallocate(theta,rhstheta,rhstheta_o)
#endif

call mpi_finalize(ierr)

end program main
