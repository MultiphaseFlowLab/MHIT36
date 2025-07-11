#define CHECK_CUDECOMP_EXIT(f) if (f /= CUDECOMP_RESULT_SUCCESS) call exit(1)

program main
use cudafor
use cudecomp
use cufft
use mpi
use velocity 
use phase
use param
use mpivar
use cudecompvar


implicit none
! grid dimensions
integer :: comm_backend
integer :: pr, pc
! cudecomp
! cuFFT
integer :: planX, planY, planZ
integer :: batchsize
integer :: status
! other variables (wavenumber, grid location)
real(8), allocatable :: x(:), kx(:)
integer :: i,j,k,il,jl,kl,ig,jg,kg,t
integer :: im,ip,jm,jp,km,kp,last
!integer, parameter :: Mx = 1, My = 0, Mz = 0
real(8), device, allocatable :: kx_d(:)
! working arrays
complex(8), allocatable :: psi(:), ua(:,:,:)
complex(8), device, allocatable :: psi_d(:)
complex(8), pointer, device, contiguous :: work_d(:), work_halo_d(:)
character(len=40) :: namefile
! Code variables

! Enable or disable phase field (acceleration eneabled by default)
#define phiflag 0

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

! hard coded, then from input
pr = 0
pc = 0
halo_ext=1
comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P

CHECK_CUDECOMP_EXIT(cudecompInit(handle, MPI_COMM_WORLD))

! config is a struct and pr and pc are the number of pencils along the two directions
! gdims is the global grid
! create an uninitialized configuration struct and initialize it to defaults using cudecompGridDescConfigSetDefaults. 
! Initializing to default values is required to ensure no entries are left uninitialized.
CHECK_CUDECOMP_EXIT(cudecompGridDescConfigSetDefaults(config))
pdims = [pr, pc] !pr and pc are the number of pencil along the different directions
config%pdims = pdims
gdims = [nx, ny, nz]
config%gdims = gdims
halo = [0, halo_ext, halo_ext] ! no halo along x neeed because is periodic and in physical space i have x-pencil
! for transpositions
config%transpose_comm_backend = comm_backend
config%transpose_axis_contiguous = .true.
! for halo exchanges
config%halo_comm_backend = CUDECOMP_HALO_COMM_MPI
! Setting for periodic halos in all directions (non required to be in config)
halo_periods = [.true., .true., .true.]

CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_DOUBLE_COMPLEX
if (comm_backend == 0) then
   options%autotune_transpose_backend = .true.
endif

! initialize cuDecomp with the config file 
CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))

! Print information on configuration
if (rank == 0) then
   write(*,"(' Running on ', i0, ' x ', i0, ' process grid ...')") config%pdims(1), config%pdims(2)
   write(*,"(' Using ', a, ' transpose backend ...')") &
            cudecompTransposeCommBackendToString(config%transpose_comm_backend)
   write(*,"(' Using ', a, ' halo backend ...')") &
            cudecompHaloCommBackendToString(config%halo_comm_backend)
endif

! get pencil info
! This function returns a pencil struct (piX, piY or piZ) that contains the shape, global lower and upper index bounds (lo and hi), 
! size of the pencil, and an order array to indicate the memory layout that will be used (to handle permuted, axis-contiguous layouts).
! Additionally, there is a halo_extents data member that indicates the depth of halos for the pencil, by axis.
! Side note:  ! cudecompGetPencilInfo(handle, grid_desc, pinfo_x, 1, [1, 1, 1]) <- in this way the x-pencil also have halo elements
! If no halo regions are necessary, a NULL pointer can be provided in place of this array (or omitted)
! Pencil info in x-configuration present in PiX (shape,lo,hi,halo_extents,size)
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piX, 1, halo))
nElemX = piX%size !<- number of total elments in x-configuratiion (including halo)
! Pencil info in Y-configuration present in PiY
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piY, 2))
nElemY = piY%size
! Pencil info in Z-configuration present in PiZ
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piZ, 3))
nElemZ = piZ%size

! Get workspace sizes for transpose (1st row) and halo (2nd row)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_desc, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_desc, 1, halo, nElemWork_halo))

! show the order 1=x, 2=y, 3=z
! x-pencils are x,y,z
! y-pencils are y,z,x
! z-pencils are z,y,x
! This is to make the setup of cuFFT easier, cuFFT cannot do FFT along inner directions and also for performance (no stride)
! if (rank .eq. 0) then
!   write(*,*) "Order in X-pencil is", piX%order(1),piX%order(2),piX%order(3)
!   write(*,*) "Order in Y-pencil is", piY%order(1),piY%order(2),piY%order(3)
!   write(*,*) "Order in Z-pencil is", piZ%order(1),piZ%order(2),piZ%order(3)
! endif

! CUFFT initialization
! Create plans (forward and backward are the same because is Z2Z!)
batchSize = piX%shape(2)*piX%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planX, nx, CUFFT_Z2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan'

! it's always 2 and 3 because y-pencil have coordinates y,z,x
batchSize = piY%shape(2)*piY%shape(3)
status = cufftPlan1D(planY, ny, CUFFT_Z2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Y plan'

! it's always 2 and 3 because y-pencil have coordinates z,y,x
batchSize = piZ%shape(2)*piZ%shape(3)
status = cufftPlan1D(planZ, nz, CUFFT_Z2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Z plan'

! define grid
allocate(x(nx),kx(nx))
x(1)= 0
do i = 2, nx
   x(i) = x(i-1) + dx
enddo
do i = 1, nx/2
   kx(i) = (i-1)*(twoPi/lx)
enddo
do i = nx/2+1, nx
   kx(i) = (i-1-nx)*(twoPi/LX)
enddo
! allocate k_d on the device (later on remove and use OpenACC + managed memory?)
allocate(kx_d, source=kx)
!########################################################################################################################################
! 1. INITIALIZATION AND cuDECOMP AUTOTUNING : END
!########################################################################################################################################






!########################################################################################################################################
! START STEP 2: ALLOCATE ARRAYS
!########################################################################################################################################
! allocate arrays
allocate(psi(max(nElemX, nElemY, nElemZ))) !largest among the pencil
allocate(psi_d, mold=psi) ! phi on device
allocate(ua(nx, piX%shape(2), piX%shape(3)))
! Pressure variable
allocate(rhsp_complex(piX%shape(1), piX%shape(2), piX%shape(3)))
allocate(rhsp(piX%shape(1), piX%shape(2), piX%shape(3))) 
allocate(p(piX%shape(1), piX%shape(2), piX%shape(3))) 
!allocate variables
!NS variables
allocate(u(piX%shape(1),piX%shape(2),piX%shape(3)),v(piX%shape(1),piX%shape(2),piX%shape(3)),w(piX%shape(1),piX%shape(2),piX%shape(3))) !velocity vector
allocate(rhsu(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(rhsu_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw_o(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(div(piX%shape(1),piX%shape(2),piX%shape(3)))
!PFM variables
#if phiflag == 1
allocate(phi(piX%shape(1),piX%shape(2),piX%shape(3)),rhsphi(piX%shape(1),piX%shape(2),piX%shape(3)),rhsphi_o(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(psidi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(normx(piX%shape(1),piX%shape(2),piX%shape(3)),normy(piX%shape(1),piX%shape(2),piX%shape(3)),normz(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(chempot(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(fxst(piX%shape(1),piX%shape(2),piX%shape(3)),fyst(piX%shape(1),piX%shape(2),piX%shape(3)),fzst(piX%shape(1),piX%shape(2),piX%shape(3))) ! surface tension forces
#endif

! allocate arrays for transpositions and halo exchanges 
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_d, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_halo_d, nElemWork_halo))
!########################################################################################################################################
! END STEP2: ALLOCATE ARRAYS
!########################################################################################################################################









!########################################################################################################################################
! START STEP 3: FLOW AND PHASE FIELD INIT
!########################################################################################################################################
! 3.1 Read/initialize from data without halo grid points (avoid out-of-bound if reading usin MPI I/O)
! 3.2 Call halo exchnages along Y and Z for u,v,w and phi
if (restart .eq. 0) then !fresh start Taylor Green or read from file in init folder
if (rank.eq.0) write(*,*) "Initialize velocity field (fresh start)"
   if (inflow .eq. 0) then
      if (rank.eq.0) write(*,*) "Initialize Taylor-green"
      do k = 1+halo_ext, piX%shape(3)-halo_ext
         kg = piX%lo(3) + k - 1 
         do j = 1+halo_ext, piX%shape(2)-halo_ext
            jg = piX%lo(2) + j - 1 
            do i = 1, piX%shape(1)
               u(i,j,k) =   sin(x(i))*cos(x(jg))*cos(x(kg))
               v(i,j,k) =  -cos(x(i))*sin(x(jg))*cos(x(kg))
               w(i,j,k) =  0.d0
            enddo
         enddo
      enddo
   endif
   if (inflow .eq. 1) then
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
      kg = piX%lo(3) + k - 1 
         do j = 1+halo_ext, piX%shape(2)-halo_ext
         jg = piX%lo(2) + j - 1 
            do i = 1, piX%shape(1)
                pos=(x(i)-lx/2)**2d0 +  (x(jg)-lx/2)**2d0 + (x(kg)-lx/2)**2d0
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

!Save initial fields (only if a fresh start)
if (restart .eq. 0) then
   if (rank.eq.0) write(*,*) "Save initial fields"
   call writefield(tstart,1)
   call writefield(tstart,2)
   call writefield(tstart,3)
   !call writefield(tstart,4)
   #if phiflag == 1
   call writefield(tstart,5)
   #endif
endif
!########################################################################################################################################
! END STEP 3: FLOW AND PHASE FIELD INIT
!########################################################################################################################################










! ########################################################################################################################################
! START TEMPORAL LOOP: STEP 4 to 9 REPEATED AT EVERY TIME STEP
! ########################################################################################################################################
! First step use Euler
alpha=1.0d0
beta=0.0d0
gumax=1.d0
tstart=tstart+1
gamma=1.d0*gumax
! Start temporal loop
do t=tstart,tfin

    if (rank.eq.0) write(*,*) "Time step",t,"of",tfin
    call cpu_time(times)

   !########################################################################################################################################
   ! START STEP 4: PHASE-FIELD SOLVER (EXPLICIT)
   !########################################################################################################################################
   #if phiflag == 1
   ! 4.1 RHS computation, no need of halo update in thoery
   ! 4.1.1 Convective term
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
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               rhsphi(i,j,k) = - (u(ip,j,k)*0.5d0*(phi(ip,j,k)+phi(i,j,k)) - u(i,j,k)*0.5d0*(phi(i,j,k)+phi(im,j,k)))*dxi  &
                               - (v(i,jp,k)*0.5d0*(phi(i,jp,k)+phi(i,j,k)) - v(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,jm,k)))*dxi  &
                               - (w(i,j,kp)*0.5d0*(phi(i,j,kp)+phi(i,j,k)) - w(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,j,km)))*dxi
            enddo
        enddo
   enddo
   !$acc end kernels

   ! 4.1.2 Compute diffusive term 
   gamma=1.0d0*gumax
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
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               rhsphi(i,j,k)=rhsphi(i,j,k)+gamma*(eps*(phi(ip,j,k)-2.d0*phi(i,j,k)+phi(im,j,k))*ddxi + &
                                                  eps*(phi(i,jp,k)-2.d0*phi(i,j,k)+phi(i,jm,k))*ddxi + &         
                                                  eps*(phi(i,j,kp)-2.d0*phi(i,j,k)+phi(i,j,km))*ddxi)
            enddo
        enddo
   enddo
   !$acc end kernels

   ! compute distance function psi (used to compute normals)
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            val = min(phi(i,j,k),1.0d0) ! avoid machine precision overshoots in phi that leads to problem with log
            psidi(i,j,k) = eps*log((val+enum)/(1.d0-val+enum))
         enddo
      enddo
   enddo
   !$acc end kernels

   ! 4.1.3. Compute Sharpening term (gradien)
   ! Substep 1 computer normals
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
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            normx(i,j,k) = (psidi(ip,j,k) - psidi(im,j,k))
            normy(i,j,k) = (psidi(i,jp,k) - psidi(i,jm,k))
            normz(i,j,k) = (psidi(i,j,kp) - psidi(i,j,km)) 
         enddo
      enddo
   enddo 
   !$acc end kernels

   ! Update normx,normy and normz halos, required to then compute normal derivative
   !$acc host_data use_device(normx)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(normy)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(normz)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 

   ! Step 2: Compute normals (1.e-16 is a numerical tollerance to avodi 0/0)
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            normod = 1.d0/(sqrt(normx(i,j,k)**2d0 + normy(i,j,k)**2d0 + normz(i,j,k)**2d0) + 1.0E-16)
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
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               ! ACDI
               rhsphi(i,j,k)=rhsphi(i,j,k)-gamma*((0.25d0*(1.d0-(tanh(0.5d0*psidi(ip,j,k)*epsi))**2)*normx(ip,j,k)- 0.25d0*(1.d0-(tanh(0.5d0*psidi(im,j,k)*epsi))**2)*normx(im,j,k))*0.5*dxi +&
                                                  (0.25d0*(1.d0-(tanh(0.5d0*psidi(i,jp,k)*epsi))**2)*normy(i,jp,k)- 0.25d0*(1.d0-(tanh(0.5d0*psidi(i,jm,k)*epsi))**2)*normy(i,jm,k))*0.5*dxi +&
                                                  (0.25d0*(1.d0-(tanh(0.5d0*psidi(i,j,kp)*epsi))**2)*normz(i,j,kp)- 0.25d0*(1.d0-(tanh(0.5d0*psidi(i,j,km)*epsi))**2)*normz(i,j,km))*0.5*dxi)
            enddo
        enddo
    enddo
    !$acc end kernels

   ! 4.2 Get phi at n+1 using AB2
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
                phi(i,j,k) = phi(i,j,k) + dt*(alpha*rhsphi(i,j,k)-beta*rhsphi_o(i,j,k))
                rhsphi_o(i,j,k)=rhsphi(i,j,k)
            enddo
        enddo
    enddo
   !$acc end kernels

   ! 4.3 Call halo exchnages along Y and Z for phi 
   !$acc host_data use_device(phi)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   #endif
   !########################################################################################################################################
   ! END STEP 4: PHASE-FIELD SOLVER 
   !########################################################################################################################################









   !########################################################################################################################################
   ! START STEP 5: USTAR COMPUTATION (PROJECTION STEP)
   !########################################################################################################################################
   ! 5.1 compute rhs 
   ! 5.2 obtain ustar and store old rhs in rhs_o
   ! 5.3 Call halo exchnages along Y and Z for u,v,w

   ! Projection step, convective terms
   ! 5.1a Convective terms NS
   ! Loop on inner nodes
   !$acc parallel loop collapse(3) private(im,jm,km)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            ! Manual periodicity ony along x (x-pencil), along y and z directions use halos
            if (ip .gt. nx) ip=1  
            if (im .lt. 1) im=nx
            ! compute the products (conservative form)
            h11 = (u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k))
            h12 = (u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
            h13 = (u(i,j,kp)+u(i,j,k))*(w(i,j,kp)+w(im,j,kp))   - (u(i,j,k)+u(i,j,km))*(w(i,j,k)+w(im,j,k))
            h21 = (u(ip,j,k)+u(ip,jm,k))*(v(ip,j,k)+v(i,j,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
            h22 = (v(i,jp,k)+v(i,j,k))*(v(i,jp,k)+v(i,j,k))     - (v(i,j,k)+v(i,jm,k))*(v(i,j,k)+v(i,jm,k))
            h23 = (w(i,j,kp)+w(i,jm,kp))*(v(i,j,kp)+v(i,j,k))   - (w(i,j,k)+w(i,jm,k))*(v(i,j,k)+v(i,j,km))
            h31 = (w(ip,j,k)+w(i,j,k))*(u(ip,j,k)+u(ip,j,km))   - (w(i,j,k)+w(im,j,k))*(u(i,j,k)+u(i,j,km))
            h32 = (v(i,jp,k)+v(i,jp,km))*(w(i,jp,k)+w(i,j,k))   - (v(i,j,k)+v(i,j,km))*(w(i,j,k)+w(i,jm,k))
            h33 = (w(i,j,kp)+w(i,j,k))*(w(i,j,kp)+w(i,j,k))     - (w(i,j,k)+w(i,j,km))*(w(i,j,k)+w(i,j,km))
            ! compute the derivative
            h11=h11*0.25d0*dxi
            h12=h12*0.25d0*dxi
            h13=h13*0.25d0*dxi
            h21=h21*0.25d0*dxi
            h22=h22*0.25d0*dxi
            h23=h23*0.25d0*dxi
            h31=h31*0.25d0*dxi
            h32=h32*0.25d0*dxi
            h33=h33*0.25d0*dxi
            ! add to the rhs
            rhsu(i,j,k)=-(h11+h12+h13)
            rhsv(i,j,k)=-(h21+h22+h23)
            rhsw(i,j,k)=-(h31+h32+h33)
         enddo
      enddo
   enddo

   ! 5.1b Compute viscous terms
   !$acc parallel loop collapse(3) private(im,jm,km)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            h11 = mu*(u(ip,j,k)-2.d0*u(i,j,k)+u(im,j,k))*ddxi
            h12 = mu*(u(i,jp,k)-2.d0*u(i,j,k)+u(i,jm,k))*ddxi
            h13 = mu*(u(i,j,kp)-2.d0*u(i,j,k)+u(i,j,km))*ddxi
            h21 = mu*(v(ip,j,k)-2.d0*v(i,j,k)+v(im,j,k))*ddxi
            h22 = mu*(v(i,jp,k)-2.d0*v(i,j,k)+v(i,jm,k))*ddxi
            h23 = mu*(v(i,j,kp)-2.d0*v(i,j,k)+v(i,j,km))*ddxi
            h31 = mu*(w(ip,j,k)-2.d0*w(i,j,k)+w(im,j,k))*ddxi
            h32 = mu*(w(i,jp,k)-2.d0*w(i,j,k)+w(i,jm,k))*ddxi
            h33 = mu*(w(i,j,kp)-2.d0*w(i,j,k)+w(i,j,km))*ddxi
            rhsu(i,j,k)=rhsu(i,j,k)+(h11+h12+h13)*rhoi
            rhsv(i,j,k)=rhsv(i,j,k)+(h21+h22+h23)*rhoi
            rhsw(i,j,k)=rhsw(i,j,k)+(h31+h32+h33)*rhoi
         enddo
      enddo
   enddo

   ! 5.1c NS forcing
   !$acc kernels
   do k = 1+halo_ext, piX%shape(3)-halo_ext
      kg = piX%lo(3) + k - 1 
      do j = 1+halo_ext, piX%shape(2)-halo_ext
         jg = piX%lo(2) + j - 1 
         do i = 1, piX%shape(1)
            ! ABC forcing
            rhsu(i,j,k)= rhsu(i,j,k) + f3*sin(k0*(x(kg)+dx/2))+f2*cos(k0*(x(jg)+dx/2))
            rhsv(i,j,k)= rhsv(i,j,k) + f1*sin(k0*(x(i) +dx/2))+f3*cos(k0*(x(kg)+dx/2))
            rhsw(i,j,k)= rhsw(i,j,k) + f2*sin(k0*(x(jg)+dx/2))+f1*cos(k0*(x(i) +dx/2))
            ! TG Forcing
            !rhsu(i,j,k)= rhsu(i,j,k) + f1*sin(k0*x(i))*cos(k0*x(j))*cos(k0*x(k))
            !rhsv(i,j,k)= rhsv(i,j,k) - f1*cos(k0*x(i))*sin(k0*x(j))*sin(k0*x(k))
         enddo
      enddo
   enddo
   !$acc end kernels

   ! Surface tension forces
   #if phiflag == 1
   !$acc kernels
   !Obtain surface tension forces evaluated at the center of the cell (same as where phi is located)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            chempot(i,j,k)=phi(i,j,k)*(1.d0-phi(i,j,k))*(1.d0-2.d0*phi(i,j,k))*epsi-eps*(phi(ip,j,k)+phi(im,j,k)+phi(i,jp,k)+phi(i,jm,k)+phi(i,j,kp)+phi(i,j,km)- 6.d0*phi(i,j,k))*ddxi
            fxst(i,j,k)=6.d0*sigma*chempot(i,j,k)*0.5d0*(phi(ip,j,k)-phi(im,j,k))*dxi
            fyst(i,j,k)=6.d0*sigma*chempot(i,j,k)*0.5d0*(phi(i,jp,k)-phi(i,jm,k))*dxi
            fzst(i,j,k)=6.d0*sigma*chempot(i,j,k)*0.5d0*(phi(i,jp,k)-phi(i,jm,k))*dxi
         enddo
      enddo
   enddo
   !$acc end kernels

   ! Update halo of fxst, fyst and fzst (required then to interpolate at velocity points)
   !$acc host_data use_device(fxst)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(fyst)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(fzst)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   
   ! Interpolate force at velocity points
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            im=i-1
            jm=j-1
            km=k-1
            if (im .lt. 1) im=nx
            rhsu(i,j,k)=rhsu(i,j,k) + 0.5d0*(fxst(im,j,k)+fxst(i,j,k))*rhoi
            rhsv(i,j,k)=rhsv(i,j,k) + 0.5d0*(fyst(i,jm,k)+fyst(i,j,k))*rhoi
            rhsw(i,j,k)=rhsw(i,j,k) + 0.5d0*(fzst(i,j,km)+fzst(i,j,k))*rhoi
         enddo
      enddo
   enddo
   !$acc end kernels   
   #endif

   ! 5.2 find u, v and w star (explicit Eulero), only in the inner nodes 
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
          do i=1,nx
              u(i,j,k) = u(i,j,k) + dt*(alpha*rhsu(i,j,k)-beta*rhsu_o(i,j,k))
              v(i,j,k) = v(i,j,k) + dt*(alpha*rhsv(i,j,k)-beta*rhsv_o(i,j,k))
              w(i,j,k) = w(i,j,k) + dt*(alpha*rhsw(i,j,k)-beta*rhsw_o(i,j,k))
              rhsu_o(i,j,k)=rhsu(i,j,k)
              rhsv_o(i,j,k)=rhsv(i,j,k)
              rhsw_o(i,j,k)=rhsw(i,j,k)
          enddo
      enddo
   enddo
   !$acc end kernels

   ! store rhs* in rhs*_o 
   ! After first step move to AB2 
   alpha=1.5d0
   beta= 0.5d0


   ! 5.3 update halos (y and z directions), required to then compute the RHS of Poisson equation because of staggered grid
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
   !########################################################################################################################################
   ! END STEP 5: USTAR COMPUTATION 
   !########################################################################################################################################









   !########################################################################################################################################
   ! START STEP 6: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################
   ! initialize rhs and analytical solution
   ! 6.1 Compute rhs of Poisson equation div*ustar: divergence at the cell center 
   ! I've done the halo updates so to compute the divergence at the pencil border i have the *star from the halo
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
          do i=1,nx
              ip=i+1
              jp=j+1
              kp=k+1
              if (ip > nx) ip=1
              rhsp(i,j,k) =               (rho*dxi/dt)*(u(ip,j,k)-u(i,j,k))
              rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(v(i,jp,k)-v(i,j,k))
              rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(w(i,j,kp)-w(i,j,k))
          enddo
      enddo
   enddo
   !$acc end kernels

   ! Feed rhsp into the Poisson solver, this part can be optimized in the future (less copies of arrays)
   !$acc kernels
   do kl = 1+halo_ext, piX%shape(3)-halo_ext
      !kg = piX%lo(3) + kl - 1 
      do jl = 1+halo_ext, piX%shape(2)-halo_ext
         !jg = piX%lo(2) + jl - 1 
         do i = 1, piX%shape(1)
            rhsp_complex(i,jl,kl) = cmplx(rhsp(i,jl,kl),0.0) !<- use this once solver is ok
            !For Poisson testing (knwon solutions)
            !rhsp_complex(i,jl,kl) = cmplx(sin(Mx*x(i)),0.0) 
            !rhsp_complex(i,jl,kl) = cmplx(sin(Mx*x(i))*sin(My*x(jg))*sin(Mz*x(kg)),0.0)  ! RHS of Poisson equation
         enddo
      enddo
   enddo
   !$acc end kernels

   ! do the FFT3D forward 
   ! psi(x,y,z) -> psi(kx,y,z)
   !$acc host_data use_device(rhsp_complex)
   status = cufftExecZ2Z(planX, rhsp_complex, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X forward error: ', status
   !$acc end host_data
   ! psi(kx,y,z) -> psi(y,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeXToY(handle, grid_desc, psi_d, psi_d, work_d, CUDECOMP_DOUBLE_COMPLEX,Pix%halo_extents, [0,0,0]))
   ! psi(y,z,kx) -> psi(ky,z,kx)
   status = cufftExecZ2Z(planY, psi_d, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y forward error: ', status
   ! psi(ky,z,kx) -> psi(z,kx,ky)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToZ(handle, grid_desc, psi_d, psi_d, work_d, CUDECOMP_DOUBLE_COMPLEX))
   ! psi(z,kx,ky) -> psi(kz,kx,ky)
   status = cufftExecZ2Z(planZ, psi_d, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Z forward error: ', status
   ! END of FFT3D forward

   block
   complex(8), device, pointer :: psi3d(:,:,:)
   real(8) :: k2
   !integer :: il, jl, ig, jg
   integer :: offsets(3), xoff, yoff
   integer :: np(3)
   np(piZ%order(1)) = piZ%shape(1)
   np(piZ%order(2)) = piZ%shape(2)
   np(piZ%order(3)) = piZ%shape(3)
   call c_f_pointer(c_devloc(psi_d), psi3d, piZ%shape)

   ! divide by -K**2, and normalize
   offsets(piZ%order(1)) = piZ%lo(1) - 1
   offsets(piZ%order(2)) = piZ%lo(2) - 1
   offsets(piZ%order(3)) = piZ%lo(3) - 1

   xoff = offsets(1)
   yoff = offsets(2)
   npx = np(1)
   npy = np(2)
   !$acc kernels
   do jl = 1, npy
      do il = 1, npx
         jg = yoff + jl
         ig = xoff + il
         do k = 1, nz
            k2 = kx_d(ig)**2 + kx_d(jg)**2 + kx_d(k)**2
            psi3d(k,il,jl) = -psi3d(k,il,jl)/k2/(int(nx,8)*int(ny,8)*int(nz,8))
         enddo
      enddo
   enddo
   !$acc end kernels
 
   ! specify mean (corrects division by zero wavenumber above)
   if (xoff == 0 .and. yoff == 0) psi3d(1,1,1) = 0.0
   end block

   ! psi(kz,kx,ky) -> psi(z,kx,ky)
   status = cufftExecZ2Z(planZ, psi_d, psi_d, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Z inverse error: ', status
   ! psi(z,kx,ky) -> psi(ky,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeZToY(handle, grid_desc, psi_d, psi_d, work_d, CUDECOMP_DOUBLE_COMPLEX))
   ! psi(ky,z,kx) -> psi(y,z,kx)
   status = cufftExecZ2Z(planY, psi_d, psi_d, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y inverse error: ', status
   ! psi(y,z,kx) -> psi(kx,y,z)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToX(handle, grid_desc, psi_d, psi_d, work_d, CUDECOMP_DOUBLE_COMPLEX, [0,0,0], Pix%halo_extents))
   ! psi(kx,y,z) -> psi(x,y,z)
   !$acc host_data use_device(rhsp_complex)
   status = cufftExecZ2Z(planX, psi_d, rhsp_complex, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X inverse error: ', status
   !$acc end host_data

   ! update halo nodes with pressure (needed for the pressure correction step), using device variable no need to use host-data
   ! Update X-pencil halos in Y and Z direction
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, rhsp_complex, work_halo_d, CUDECOMP_DOUBLE_COMPLEX, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, rhsp_complex, work_halo_d, CUDECOMP_DOUBLE_COMPLEX, piX%halo_extents, halo_periods, 3))

   !D2H transfer
   !psi = psi_d


   ! check against analytical solution
   block
   complex(8), pointer :: psi3(:,:,:)
   real(8) :: err, maxErr = -1.0, err2
   call c_f_pointer(c_loc(psi), psi3, [piX%shape(1), piX%shape(2), piX%shape(3)])

   !!Check errro on complex (ua and psi3 are complex)
   !do kl = 1, piX%shape(3)
   !   do jl = 1, piX%shape(2)
   !      do i = 1, piX%shape(1)
   !         err = abs(ua(i,jl,kl)-psi3(i,jl,kl))
   !         if (err > maxErr) maxErr = err
   !      enddo
   !   enddo
   !enddo

   ! take back p from psi3
   ! i can span also the halo because they have been already updated
   !$acc kernels
   do kl = 1, piX%shape(3)
      do jl = 1, piX%shape(2)
         do i = 1, piX%shape(1)
            p(i,jl,kl) = real(rhsp_complex(i,jl,kl))
         enddo
      enddo
   enddo
   !$acc end kernels

   end block

   !########################################################################################################################################
   ! END STEP 7: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################








   !########################################################################################################################################
   ! START STEP 8: VELOCITY CORRECTION
   ! ########################################################################################################################################
   ! 8.1 Correct velocity 
   ! 8.2 Remove mean velocity if using ABC forcing
   ! 8.3 Call halo exchnages along Y and Z for u,v,w
   ! Correct velocity, pressure has also the halo
   !$acc kernels 
   do k = 1+halo_ext, piX%shape(3)-halo_ext
      do j = 1+halo_ext, piX%shape(2)-halo_ext
         do i = 1, piX%shape(1) ! equal to nx (no halo on x)
              im=i-1
              jm=j-1
              km=k-1
              if (im < 1) im=nx
              u(i,j,k)=u(i,j,k) - dt/rho*(p(i,j,k)-p(im,j,k))*dxi
              v(i,j,k)=v(i,j,k) - dt/rho*(p(i,j,k)-p(i,jm,k))*dxi
              w(i,j,k)=w(i,j,k) - dt/rho*(p(i,j,k)-p(i,j,km))*dxi
          enddo
      enddo
   enddo
   !$acc end kernels 

   ! Remove mean velocity (get local mean of the rank)
   umean=0.d0
   vmean=0.d0
   wmean=0.d0
   !$acc kernels 
   do k = 1+halo_ext, piX%shape(3)-halo_ext
      do j = 1+halo_ext, piX%shape(2)-halo_ext
         do i = 1, piX%shape(1) ! equal to nx (no halo on x)
              umean=umean + u(i,j,k)
              vmean=vmean + v(i,j,k)
              wmean=wmean + w(i,j,k)
          enddo
      enddo
   enddo
   !$acc end kernels 

   ! Divide by total number of points in the pencil
   umean=umean/nx/(piX%shape(2)-2*halo_ext)/(piX%shape(3)-2*halo_ext)
   vmean=vmean/nx/(piX%shape(2)-2*halo_ext)/(piX%shape(3)-2*halo_ext)
   wmean=wmean/nx/(piX%shape(2)-2*halo_ext)/(piX%shape(3)-2*halo_ext)

   ! Find global mean (MPI_SUM and then divide by number of ranks)
   call MPI_Allreduce(umean,gumean,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(vmean,gvmean,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(wmean,gwmean,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

   ! remove mean value
   !$acc kernels 
   u=u-(gumean/ranks)
   v=v-(gvmean/ranks)
   w=w-(gwmean/ranks)
   !$acc end kernels 

   ! 8.3 update halos (y and z directions), required to then compute the RHS of Poisson equation because of staggered grid
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

   ! find local maximum velocity
   uc=maxval(u)
   vc=maxval(v)
   wc=maxval(w)
   umax=max(wc,max(uc,vc))
   call MPI_Allreduce(umax,gumax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)
   cou=gumax*dt*dxi
   if (rank.eq.0) then
      write(*,*) "CFL (max among tasks)", cou
      if (cou .gt. 7) stop
   endif

   call cpu_time(timef)
   if (rank.eq.0) print '(" Time elapsed = ",f8.1," ms")',1000*(timef-times)
   !########################################################################################################################################
   ! END STEP 8: VELOCITY CORRECTION  
   !########################################################################################################################################




   !########################################################################################################################################
   ! START STEP 9: OUTPUT FIELDS 
   ! ########################################################################################################################################
   if (mod(t,dump) .eq. 0) then
      if (rank .eq. 0) write(*,*) "Saving output files"
          ! write velocity and pressure fiels (1-4)
         call writefield(t,1)
         call writefield(t,2)
         call writefield(t,3)
         !call writefield(t,4)
         #if phiflag == 1
         ! write phase-field (5)
         call writefield(t,5)
         #endif
   endif
   !########################################################################################################################################
   ! END STEP 9: OUTPUT FIELDS N  
   !########################################################################################################################################


enddo

deallocate(u,v,w)
deallocate(rhsu,rhsv,rhsw)
deallocate(rhsu_o,rhsv_o,rhsw_o)
deallocate(phi,rhsphi,rhsphi_o,normx,normy,normz)

call mpi_finalize(ierr)

end program main
