subroutine writefield(t,fieldn)
! Output field, file is written in the /output folder 

use velocity
use phase
use mpi
use mpivar
use param
use cudecompvar

implicit none

integer :: g_size(3),p_size(3),fstart(3)
integer :: t,fieldn
character(len=40) :: namefile
integer(mpi_offset_kind) :: offset=0
integer :: f_handle ! file handle
integer :: ftype
double precision, allocatable :: out(:,:,:)

! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means w
! fieldn=4 means p
! fieldn=5 means phi

! define basic quantities to be used later (gloabl and pencil size)
g_size=[nx, ny, nz] ! global size
p_size=[piX%shape(1), piX%shape(2)-2*halo_ext, piX%shape(3)-2*halo_ext] !<- pencil has no halo along x
fstart=[piX%lo(1)-1,piX%lo(2)-1,piX%lo(3)-1]
! for debug
!write(*,*) "g_size", g_size
!write(*,*) "p_size", p_size
!write(*,*) "fstart", fstart
allocate(out(p_size(1),p_size(2),p_size(3))) !<- halo removed 
 
!write(*,*) "in readwrite"

if (fieldn .eq. 1) then
  out=u(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext) !<- out only the inner parts (no halo)
  write(namefile,'(a,i8.8,a)') './output/u_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  call mpi_file_write_all(f_handle,out,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
endif

if (fieldn .eq. 2) then
  out=v(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext) !<- out only the inner parts (no halo)
  write(namefile,'(a,i8.8,a)') './output/v_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  call mpi_file_write_all(f_handle,out,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
endif

if (fieldn .eq. 3) then
  out=w(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext) !<- out only the inner parts (no halo)
  write(namefile,'(a,i8.8,a)') './output/w_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  call mpi_file_write_all(f_handle,out,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
endif

if (fieldn .eq. 4) then
  out=p(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext) !<- out only the inner parts (no halo)
  write(namefile,'(a,i8.8,a)') './output/p_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  call mpi_file_write_all(f_handle,out,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
endif

if (fieldn .eq. 5) then
  out=phi(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext) !<- out only the inner parts (no halo)
  write(namefile,'(a,i8.8,a)') './output/phi_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  call mpi_file_write_all(f_handle,out,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
endif

deallocate(out)

end subroutine










subroutine readfield(fieldn)
! Input field, file is written in the /input folder 

use velocity
use phase
use mpi
use mpivar
use param
use cudecompvar

implicit none

integer :: g_size(3),p_size(3),fstart(3)
integer :: fieldn
character(len=40) :: namefile
integer(mpi_offset_kind) :: offset=0
integer :: f_handle ! file handle
integer :: ftype
double precision, allocatable :: in(:,:,:)

! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means w
! fieldn=4 means p
! fieldn=5 means phi

! define basic quantities to be used later (gloabl and pencil size)
g_size=[nx, ny, nz] ! global size
p_size=[piX%shape(1), piX%shape(2)-2*halo_ext, piX%shape(3)-2*halo_ext] !<- pencil has no halo along x
fstart=[piX%lo(1)-1,piX%lo(2)-1,piX%lo(3)-1] !<- MPI is in C and index start from 0 (not 1)
! for debug
!write(*,*) "g_size", g_size
!write(*,*) "p_size", p_size
!write(*,*) "fstart", fstart
allocate(in(p_size(1),p_size(2),p_size(3))) !<- no halos read
 
!write(*,*) "in readwrite"

if (fieldn .eq. 1) then
  namefile='./input/u.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  u(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 2) then
  namefile='./input/v.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  v(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 3) then
  namefile='./input/w.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  w(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 4) then
  namefile='./input/p.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  p(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 5) then
  namefile='./input/phi.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  phi(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

deallocate(in)

end subroutine










subroutine readfield_restart(t,fieldn)
! Used in case of restart, file is read from the multi/output folder (iteration tstart must be present!)

use velocity
use phase
use mpi
use mpivar
use param
use cudecompvar

implicit none

integer :: g_size(3),p_size(3),fstart(3)
integer :: t,fieldn
character(len=40) :: namefile
integer(mpi_offset_kind) :: offset=0
integer :: f_handle ! file handle
integer :: ftype
double precision, allocatable :: in(:,:,:)

! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means w
! fieldn=4 means p
! fieldn=5 means phi

! define basic quantities to be used later (gloabl and pencil size)
g_size=[nx, ny, nz] ! global size
p_size=[piX%shape(1), piX%shape(2)-2*halo_ext, piX%shape(3)-2*halo_ext] !<- pencil has no halo along x
fstart=[piX%lo(1)-1,piX%lo(2)-1,piX%lo(3)-1] !<- MPI is in C and index start from 0 (not 1)
! for debug
!write(*,*) "g_size", g_size
!write(*,*) "p_size", p_size
!write(*,*) "fstart", fstart
allocate(in(p_size(1),p_size(2),p_size(3))) !<- no halos read
 
!write(*,*) "in readwrite"

if (fieldn .eq. 1) then
  write(namefile,'(a,i8.8,a)') './output/u_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  u(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 2) then
  write(namefile,'(a,i8.8,a)') './output/v_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  v(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 3) then
  write(namefile,'(a,i8.8,a)') './output/w_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  w(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 4) then
  write(namefile,'(a,i8.8,a)') './output/p_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  p(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

if (fieldn .eq. 5) then
  write(namefile,'(a,i8.8,a)') './output/phi_',t,'.dat'
  call mpi_file_open(MPI_COMM_WORLD,namefile,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
  call mpi_type_create_subarray(3,g_size,p_size,fstart,mpi_order_fortran,mpi_double_precision,ftype,ierr)
  call mpi_type_commit(ftype,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
  !call mpi_file_read_all(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_read(f_handle,in,p_size(1)*p_size(2)*p_size(3),mpi_double_precision,mpi_status_ignore,ierr)
  call mpi_file_close(f_handle,ierr)
  phi(1:nx,1+halo_ext:piX%shape(2)-halo_ext,1+halo_ext:piX%shape(3)-halo_ext)=in !<- read only the inner parts (no halo) u has halos; in no halos
endif

deallocate(in)

end subroutine