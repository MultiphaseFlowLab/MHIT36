module param
 implicit none
 double precision, parameter :: pi=4.0d0*atan(1.0d0)
 double precision, parameter :: twopi=8.0d0*atan(1.0d0)
 double precision, parameter :: n6opi=6.0d0/pi
 integer, parameter:: clen=300 !char length
 integer :: nx
 integer :: begin, finish, step
 character(len=clen) :: rootpath
 double precision ::dx,lx
end module param
!*********************************************
module flowvars
 implicit none
 double precision, allocatable, dimension(:,:,:) :: phi
end module flowvars
!*********************************************