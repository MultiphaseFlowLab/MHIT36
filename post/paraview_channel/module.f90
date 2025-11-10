module commondata
 integer :: nx, ny, nz
 integer :: nstart,nend,dump
 integer :: uflag, vflag, wflag, phiflag, thetaflag, nfields
 double precision, parameter :: pi=3.14159265358979
 double precision :: dx,dy,dz,lx,ly,lz,csi,zk
 double precision, allocatable, dimension(:) :: x,y,z
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi,theta
end module commondata
