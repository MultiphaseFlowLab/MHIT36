subroutine get_interface(nstep)
  use param
  use flowvars
  implicit none
  integer :: nstep
  integer :: i,j,k,id,jd,kd
  integer :: top(nx,ny,nz),s_drop(nx,ny,nz)
  double precision :: area 

  area=0.d0

  do k=1,nz
    do j=1,ny
      do i=1,nx
        if( (phi(i,j,k).ge.0.05d0) .and. (phi(i,j,k).le.0.95d0)) then
          area=area + dx*dy*dz
        endif
      end do
    end do
  end do

  area=area/(5.8888*0.51*dx)
  
  open(2,file='area.dat',access='append',form='formatted',status='old')
write(2,'(i16,2x,f16.6)') nstep, area
  close(2,status='keep')
  
  return
end

