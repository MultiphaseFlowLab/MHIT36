module commondata
 integer :: nx, ny, nz
 integer :: nstart,nend,dump
 integer :: uflag, phiflag, nfields
 double precision, parameter :: pi=3.14159265358979
 double precision :: dx,dy,dz,zk
 double precision :: lx,ly,lz,csi
 double precision, allocatable, dimension(:) :: x,y,z
 double precision, allocatable, dimension(:,:) :: rms,mean,skw,flt
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi
end module commondata
