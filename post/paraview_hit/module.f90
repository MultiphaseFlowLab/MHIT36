module commondata
 integer :: nx, ny, nz
 integer :: nstart,nend,dump
 integer :: uflag, vflag, wflag, phiflag, nfields
 double precision, parameter :: pi=3.14159265358979
 double precision :: dx,dy,dz
 double precision, allocatable, dimension(:) :: x,y,z
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi
end module commondata
