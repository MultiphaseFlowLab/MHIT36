module commondata
 integer :: nx
 integer :: nstart,nend,ndump,sdump

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt,dx,lx,dxi,epsilon,nu
 double precision, allocatable, dimension(:,:,:) :: u,v,w
 double precision, allocatable, dimension(:,:,:) :: dudx,dudy,dudz
 double precision, allocatable, dimension(:,:,:) :: dvdx,dvdy,dvdz
 double precision, allocatable, dimension(:,:,:) :: dwdx,dwdy,dwdz
 double precision, allocatable,dimension(:) :: x
end module commondata


