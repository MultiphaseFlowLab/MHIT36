!==============================================================================
!  union_find_gpu.f90 – OpenACC Union‑Find connected‑component
!  labelling .
!  ---------------------------------------------------------------------------
!  • Compile with 26‑point connectivity (faces + edges + corners + periodic):
!      nvfortran -acc -Minfo=accel -O3 -D USE_26 -o union26 union_find_gpu.f90
!  • Compile with 6‑point connectivity (faces only + periodic):
!      nvfortran -acc -Minfo=accel -O3 -D USE_6  -o union6  union_find_gpu.f90
!==============================================================================
program union_find_gpu
   use openacc
   implicit none

   !---------------- user parameters --------------------------------------------
   logical, parameter :: readf = .true.          ! read external binary field?
   integer, parameter :: nx = 512                ! cubic grid size
   real   , parameter :: thresh = 0.5            ! threshold for mask
   integer :: first_step = 26000                     ! start processing from
   integer :: last_step = 26000
   integer :: step = 1000                        !
   !-----------------------------------------------------------------------------

   integer :: n, iter, max_iter, min_label, iter2, nn
   integer :: local_count(1)
   logical :: converged
   integer :: cx, cy, cz, r2
   integer :: lblmax

   integer, parameter :: tile_size=16
   integer :: sub_id(0:tile_size+1,0:tile_size+1,0:tile_size+1)
   logical :: sub_mask(0:tile_size+1,0:tile_size+1,0:tile_size+1)
   logical :: not_converged,valid
   integer :: num_not_converged
   integer :: i, j, k, ii, jj, kk, iii, jjj, kkk, new_id

   ! loop over more timesteps
   integer :: it
   character(len=200) :: phi_file, u_file, v_file, w_file

   ! Timing variables
   real :: t_start, t_end, t_elapsed

   ! double precision id array for saving
   real(kind=8), allocatable :: id_double(:,:,:)

   !==================== neighbour stencil =====================================
   #ifdef USE_26
      integer, parameter :: nneigh = 26
      integer, parameter :: di(nneigh) = [  1,-1, 0, 0, 0, 0,   &
            1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0, &
            1, 1, 1, 1,-1,-1,-1,-1 ]
      integer, parameter :: dj(nneigh) = [  0, 0, 1,-1, 0, 0,   &
            1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, &
            1, 1,-1,-1, 1, 1,-1,-1 ]
      integer, parameter :: dk(nneigh) = [  0, 0, 0, 0, 1,-1,   &
            0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, &
            1,-1, 1,-1, 1,-1, 1,-1 ]
   #elif USE_6
      integer, parameter :: nneigh = 6              ! faces only
      integer, parameter :: di(nneigh) = [  1,-1, 0, 0, 0, 0 ]
      integer, parameter :: dj(nneigh) = [  0, 0, 1,-1, 0, 0 ]
      integer, parameter :: dk(nneigh) = [  0, 0, 0, 0, 1,-1 ]
   #endif
  !=============================================================================

   integer :: ios
   character(len=100) :: namefile

   ! field arrays
   logical,  allocatable :: mask(:,:,:)
   integer,  allocatable :: id(:,:,:)

   ! comment if needed !!!!!!!!!!!!!!!!!!!!!!!!
      integer,  allocatable :: new_id_vec(:,:,:)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   real(kind=8), allocatable :: phi(:,:,:)
   real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   logical, allocatable :: valid_subblock(:,:,:)
   ! stats arrays
   integer, allocatable :: droplet_volume(:)
   real    , allocatable :: xcm(:), ycm(:), zcm(:)
   real    , allocatable :: vxcm(:), vycm(:), vzcm(:)
   double precision, allocatable :: sum_sinx(:), sum_cosx(:), sum_siny(:), sum_cosy(:), sum_sinz(:), sum_cosz(:), x(:)! added for periodicity
   double precision :: lx, pi, dx
   integer :: n_droplets
   real    :: xcom, ycom, zcom
   real    :: vx, vy, vz
   integer :: d, ni, nj, nk, num_empty

   allocate(phi(nx,nx,nx), mask(0:nx+1,0:nx+1,0:nx+1), id(0:nx+1,0:nx+1,0:nx+1))
   allocate(u(nx,nx,nx), v(nx,nx,nx), w(nx,nx,nx))
   allocate(valid_subblock(0:nx/tile_size-1,0:nx/tile_size-1,0:nx/tile_size-1))
   phi = 0.0
   u = 0.0
   v = 0.0
   w = 0.0
   valid_subblock=.false.

   ! comment if needed !!!!!!!!!!!!!!!!!!!!!!!!
   allocate(new_id_vec(0:nx+1,0:nx+1,0:nx+1))
   new_id_vec = 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do it = first_step, last_step, step
   write(phi_file, '("../hit/output/phi_", I8.8,".dat")') it
   write(u_file, '("../hit/output/u_", I8.8,".dat")') it
   write(v_file, '("../hit/output/v_", I8.8,".dat")') it
   write(w_file, '("../hit/output/w_", I8.8,".dat")') it
   if (readf) then
!      open(55, file='../../hit/output/phi_00000000.dat', form='unformatted', access='stream', status='old')
      open(55, file=phi_file, form='unformatted', access='stream', status='old')
      read(55) phi
      close(55)
      open(56, file=u_file, form='unformatted', access='stream', status='old')
      read(56) u
      close(56)
      open(57, file=v_file, form='unformatted', access='stream', status='old')
      read(57) v
      close(57)
      open(58, file=w_file, form='unformatted', access='stream', status='old')
      read(58) w
      close(58)
   else
      ! simple synthetic sphere for test
      cx = nx/4; cy = nx/3; cz = nx/2; r2 = (nx/10)**2
      do k=1,nx; do j=1,nx; do i=1,nx
         if ((i-cx)**2+(j-cy)**2+(k-cz)**2 <= r2) phi(i,j,k) = 1.0
      end do; end do; end do
   end if
   
   print *, "Processing timestep:", it

   !$acc data copyin(phi) copy(id) create(mask) copyin(valid_subblock)

   !$acc parallel loop collapse(3) present(phi,mask)
   do k=0,nx+1; do j=0,nx+1; do i=0,nx+1
      ! mask(i,j,k) = .false.
      if (i>0 .and. i<nx+1 .and. j>0 .and. j<nx+1 .and. k>0 .and. k<nx+1) then
         mask(i,j,k) = (phi(i,j,k) >= thresh)
      else  
         mask(i,j,k) = .false.
      end if
   end do; end do; end do

   !$acc kernels
   mask(0,:,:) = mask(nx,:,:)
   mask(:,0,:) = mask(:,nx,:)
   mask(:,:,0) = mask(:,:,nx)
   mask(nx+1,:,:) = mask(1,:,:)
   mask(:,nx+1,:) = mask(:,1,:)
   mask(:,:,nx+1) = mask(:,:,1)
   !$acc end kernels

   !$acc parallel loop collapse(3) present(mask,id)
   do k=0,nx+1; do j=0,nx+1; do i=0,nx+1
      if (mask(i,j,k)) then
         id(i,j,k) = (k)*(nx+1)*(nx+1) + (j)*(nx+1) + i+1
      !   write(*,*) id(i,j,k)
      else
         id(i,j,k) = 0
      end if
   end do; end do; end do

!======================== GPU data region ====================================
   call cpu_time(t_start)
   num_empty=0
   !$acc parallel loop collapse(3) reduction(+:num_empty)
   do k=0,nx-1,tile_size
      do j=0,nx-1,tile_size
         do i=0,nx-1,tile_size
            valid=.false.
            !$acc loop vector collapse(3) reduction (.or.: valid)
            do kk=1,tile_size
               do jj=1,tile_size
                  do ii=1,tile_size
                     valid = valid .or. mask(i+ii,j+jj,k+kk)
                  enddo
               enddo
            enddo
            if (.not. valid) num_empty=num_empty+1
            valid_subblock(i/tile_size,j/tile_size,k/tile_size) = valid
         enddo
      enddo
   enddo
   not_converged = .true.

   iter = 0; max_iter = 10000
   do while (not_converged .and. iter < max_iter)
      iter = iter + 1;
      write(*,*) 'Iteration:', iter
      not_converged = .false.
      local_count(1)=0
      !$acc parallel loop collapse(3) private(sub_mask, sub_id) reduction(.or.:not_converged) &
      !$acc vector_length(min(256,tile_size*tile_size*tile_size))
      do k=0,nx-1,tile_size
         do j=0,nx-1,tile_size
            do i=0,nx-1,tile_size
               if (valid_subblock(i/tile_size,j/tile_size,k/tile_size) ) then
                  !$acc cache(sub_mask,sub_id)
                  !$acc loop vector collapse(3)
                  do kk=0,tile_size+1
                     do jj=0,tile_size+1
                        do ii=0,tile_size+1 
                           sub_mask(ii,jj,kk)=mask(i+ii,j+jj,k+kk)
                           sub_id(ii,jj,kk)=id(i+ii,j+jj,k+kk)
                        enddo
                     enddo
                  enddo
                  num_not_converged = 1 
                  iter2 = 0
                  do while (num_not_converged.gt.0)
                     iter2 = iter2 + 1
                     ! write(*,*) '  Local iteration:', iter2
                     num_not_converged = 0
                     !$acc loop vector collapse(3) reduction(+:num_not_converged)
                     do kk=1,tile_size
                        do jj=1,tile_size
                           do ii=1,tile_size
                              if (sub_mask(ii,jj,kk)) then
                                 min_label = sub_id(ii,jj,kk)
                                 ! write(*,*) '  Local iteration:', i
                                 do n = 1, nneigh
                                    ni = ii+di(n)
                                    nj = jj+dj(n)
                                    nk = kk+dk(n)
                                    if (sub_mask(ni,nj,nk)) then
                                       min_label = min(min_label, sub_id(ni,nj,nk))
                                    endif
                                 end do
                                 if (min_label < sub_id(ii,jj,kk)) then
                                    sub_id(ii,jj,kk) = min_label;
                                    num_not_converged = num_not_converged + 1
                                 end if
                              end if
                           enddo
                        enddo
                     enddo
                     if (num_not_converged.gt.0) not_converged=.true.
                  enddo !while
                  !$acc loop vector collapse(3)
                  do kk=1,tile_size
                     do jj=1,tile_size
                        do ii=1,tile_size
                           mask(i+ii,j+jj,k+kk)=sub_mask(ii,jj,kk)
                           id(i+ii,j+jj,k+kk)=sub_id(ii,jj,kk)
                        enddo
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
      !   call MPI_BUFFER_EX()
      !$acc kernels
      id(0,:,:) = id(nx,:,:)
      id(:,0,:) = id(:,nx,:)
      id(:,:,0) = id(:,:,nx)
      id(nx+1,:,:) = id(1,:,:)
      id(:,nx+1,:) = id(:,1,:)
      id(:,:,nx+1) = id(:,:,1)
      !$acc end kernels
   end do


   call cpu_time(t_end)
   t_elapsed = t_end - t_start
   !$acc end data

   print *, 'Execution time (s):', t_elapsed
   !====================== end GPU region =======================================

   print *, 'Connectivity:', nneigh, ' neighbours; iterations:', iter

   !---------------- host post-processing --------------------------------------
   lblmax = maxval(id)
   allocate(droplet_volume(lblmax)); droplet_volume = 0
   allocate(xcm(lblmax)); xcm = 0.0
   allocate(ycm(lblmax)); ycm = 0.0
   allocate(zcm(lblmax)); zcm = 0.0
   allocate(vxcm(lblmax)); vxcm = 0.0
   allocate(vycm(lblmax)); vycm = 0.0
   allocate(vzcm(lblmax)); vzcm = 0.0
   allocate(sum_sinx(lblmax)); sum_sinx = 0.0
   allocate(sum_cosx(lblmax)); sum_cosx = 0.0
   allocate(sum_siny(lblmax)); sum_siny = 0.0
   allocate(sum_cosy(lblmax)); sum_cosy = 0.0
   allocate(sum_sinz(lblmax)); sum_sinz = 0.0
   allocate(sum_cosz(lblmax)); sum_cosz = 0.0
   allocate(x(nx)); x = 0.0

   pi = 4.0d0*atan(1.0d0)
   lx = 2.d0*pi
   dx = lx/real(nx)
   x(1) = dx/2.d0
   do i=2,nx
      x(i) = x(i-1) + dx
   end do

   write(*,*) "data", x(1), x(2),x(511)

   do k=1,nx; do j=1,nx; do i=1,nx
      n = id(i,j,k)
      if (n > 0) then
         droplet_volume(n) = droplet_volume(n) + 1
         sum_sinx(n) = sum_sinx(n) + sin(2.0d0*pi*x(i)/lx)
         sum_cosx(n) = sum_cosx(n) + cos(2.0d0*pi*x(i)/lx)
         sum_siny(n) = sum_siny(n) + sin(2.0d0*pi*x(j)/lx)
         sum_cosy(n) = sum_cosy(n) + cos(2.0d0*pi*x(j)/lx)
         sum_sinz(n) = sum_sinz(n) + sin(2.0d0*pi*x(k)/lx)
         sum_cosz(n) = sum_cosz(n) + cos(2.0d0*pi*x(k)/lx)
         vxcm(n) = vxcm(n) + u(i, j, k)
         vycm(n) = vycm(n) + v(i, j, k)
         vzcm(n) = vzcm(n) + w(i, j, k)
      end if
   end do; end do; end do

   n_droplets = 0
   do n = 1, lblmax
      if (droplet_volume(n) > 0) then
         n_droplets = n_droplets + 1
         xcm(n) = lx*atan2(sum_sinx(n), sum_cosx(n)) / (2.0*pi)
         ycm(n) = lx*atan2(sum_siny(n), sum_cosy(n)) / (2.0*pi)
         zcm(n) = lx*atan2(sum_sinz(n), sum_cosz(n)) / (2.0*pi)
         ! ensure positive coordinates
         if (xcm(n) < 0.0) xcm(n) = xcm(n) + lx
         if (ycm(n) < 0.0) ycm(n) = ycm(n) + lx
         if (zcm(n) < 0.0) zcm(n) = zcm(n) + lx
         ! normalize the center of mass for the vloume
         !xcom = xcm(n) / droplet_volume(n)
         !ycom = ycm(n) / droplet_volume(n)
         !zcom = zcm(n) / droplet_volume(n)
         vx = vxcm(n) / droplet_volume(n)
         vy = vycm(n) / droplet_volume(n)
         vz = vzcm(n) / droplet_volume(n)
         print *, 'Drop', n_droplets, ': vox', droplet_volume(n), 'centre', xcm(n), ycm(n), zcm(n), 'velocities', vx, vy, vz
      end if
   end do
   print *, 'Total droplets:', n_droplets
   print*, 'empty blocks',num_empty

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! rescale IDs (no needed, remove wioth related variables for performance) !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   new_id = 1
   do n = 1, lblmax
      if (droplet_volume(n) > 0) then
            do k=1,nx; do j=1,nx; do i=1,nx
               if (id(i,j,k) == n) new_id_vec(i,j,k) = new_id
            end do; end do; end do
            new_id = new_id + 1
      end if
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! filepath: /leonardo_scratch/fast/IscrB_SONORA/DIEGO/UNIONFIND/main.f90
! Save id array (no halos) as binary file for post-processing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Convert id array to double precision array and save!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(id_double(1:nx,1:nx,1:nx))

do jjj = 1, nx
   do iii = 1, nx
      do kkk = 1, nx
         id_double(iii,jjj,kkk) = real(new_id_vec(iii,jjj,kkk), kind=8)
      end do
   end do
end do
! id_double = real(id(1:nx,1:nx,1:nx), kind=8)

!namefile = './id_00212000.dat'
!
!open(unit=99, file=namefile, form='unformatted', access='stream', status='replace', iostat=ios)
!if (ios /= 0) then
!   print *, 'Error opening file: ', trim(namefile)
!else
!   write(99) id_double
!   close(99)
!   write(*,*) 'ID array (double) saved to file:', trim(namefile)
!end if
!
deallocate(id_double)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call write_centers(it, droplet_volume, xcm, ycm, zcm, vxcm, vycm, vzcm, lblmax, dx)
! write(*,*) 'ID array saved to file:', trim(namefile)
end do
end program union_find_gpu

subroutine write_centers(timestep, droplet_volume, xcm, ycm, zcm, vxcm, vycm, vzcm, lblmax,dx)
   implicit none
   integer, intent(in) :: timestep, lblmax
   integer, intent(in) :: droplet_volume(lblmax)
   real,    intent(in) :: xcm(lblmax), ycm(lblmax), zcm(lblmax)
   real,    intent(in) :: vxcm(lblmax), vycm(lblmax), vzcm(lblmax)
   double precision,    intent(in) :: dx

   character(len=200) :: fname
   integer :: n, ios, ndrops

   write(fname,'("center_timestep_",I8.8,".dat")') timestep
   open(unit=77, file=fname, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      print *, "Error opening ", trim(fname)
      return
   end if

   write(77,*) "# ID x y z u v w volume"

   ndrops = 0
   do n = 1, lblmax
      if (droplet_volume(n) > 0) then
         ndrops = ndrops + 1
         write(77,'(I8,1X,7E20.10)') ndrops,  &
              xcm(n),  &
              ycm(n),  &
              zcm(n),  &
              vxcm(n)/droplet_volume(n), &
              vycm(n)/droplet_volume(n), &
              vzcm(n)/droplet_volume(n), real(droplet_volume(n))*dx*dx*dx
      end if
   end do
   close(77)

   ! Append number of droplets to drop_count.dat
   open(unit=88, file="drop_count.dat", status="unknown", position="append")
   write(88,'(I12,2X,I8)') timestep, ndrops
   close(88)

   ! Append total volume to time_check.dat
   open(unit=89, file="time_check.dat", status="unknown", position="append")
   write(89,'(I12,2X,I8,2X,E20.10)') timestep, ndrops, real(sum(droplet_volume))
   close(89)

   print *, "Wrote ", ndrops, " droplets to ", trim(fname)

end subroutine write_centers
