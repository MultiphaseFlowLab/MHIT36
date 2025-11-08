subroutine calculate_deq(s_drop,nstep)
 use param
 use flowvars
 implicit none 
 integer :: nstep
 integer :: s_drop(nx,nx,nx)
 character(len=clen) :: time,namefile
 !double precision :: sdrop_volume(nx,nx,nx)
 double precision :: vol,deq

 !DZ: This is the previous approach
 !    I think this is superfluous. 
 !    The mask s_drop already contains the info of sdrop_volume
 !  sdrop_volume = 0.0d0
 !  do kk=1,nx
 !     do jj=1,nx
 !       do ii=1,nx
 !         if (s_drop(ii,jj,kk).eq.1)then
 !           ! store volume cells
 !           if(phi(ii,jj,kk).ge.0.5d0) sdrop_volume(ii,jj,kk)=1.0d0
 !         endif
 !       enddo
 !     enddo
 !   enddo
 
 !  ! Calculate volume of the single drop
 !  vol=SUM(sdrop_volume)

 vol=sum(dble(s_drop))

 ! Equivalent Diameter; twice the radius of the equivalent sphere ( Vsphere=(4/3)*pi*R^3 )
 ! The volume has to be multiplied by dx^3. Therefore I can just multiply by dx outside.
 deq=dx*(n6opi*vol)**(1.0d0/3.0d0)
    
 write(*,'(1x,a,E16.6)') 'Diameter',deq
 write(*,*)
    
 ! output file
 write(time,'(i8.8)') nstep
 namefile='deq_'//trim(time)//'.dat'
 open(3,file=trim(namefile),access='append',form='formatted',status='unknown')
 write(3,'(E12.5)') deq
 close(3,status='keep')
    

    
    
 return
end subroutine
!***************************************
! OLD stuff
! use commondata
! use fields

! !variables declaration
! integer :: i,j,k
! double precision :: sdrop_volume(nx,nz,ny),sdrop_layer(nx,nz,ny)
! double precision :: sumxy(nz)
! double precision :: dx,dy,vol,req,deq!,aeq,int,area,sph
! character(len=40) :: namefile
! character(len=8) :: time


! sdrop_volume=0.0d0
! !sdrop_layer=0.0d0
! ! Store bubble volume cells
! do j=1,ny
!   do k=1,nz
!     do i=1,nx
!       if (s_drop(i,k,j).eq.1)then
!         ! store volume cells
!         if(phi(i,k,j).ge.0.0d0) sdrop_volume(i,k,j)=1.0d0
!       endif
!     enddo
!   enddo
! enddo

! dx=xl/dble(nx-1)
! dy=yl/dble(ny-1)



! !!! Calculate volume
! vol=SUM(sdrop_volume)
! ! sum in x and y directions
! sumxy=sum(sum(sdrop_volume,3),1)
! ! check if the drop really exist:
! !(if its volume is zero it means that the grid is not able to capture it,
! !so we can remove it, as it is not essential for the calculation sphericity statistics)
! if (sum(sumxy).eq.0.0) then

!   drop_count=drop_count-1

! else

!   ! integrate columns
!   do k=1,nz-1
!     vol=vol+sumxy(k)*(z(k)-z(k+1))*dx*dy
!   enddo
!   ! radius of the equivalent sphere ( Vsphere=(4/3)*pi*R^3 )
!   req=((3.0d0*vol)/(4.0d0*pi))**(1.0d0/3.0d0)

!   ! diameter of the equivalent sphere
!   deq=2.0d0*req
!   ! -> wall units
!   deq=deq*re

!   !!! Calculate equivalent interface area ( Asphere=4*pi*R^2 )
!   !aeq=4.0d0*pi*req**2.0d0

!   !!! Calculate real interface area
!   !int=0.0d0
!   ! sum in x and y directions
!   !sumxy=sum(sum(sdrop_layer,3),1)
!   ! integrate columns
!   !do k=1,nz-1
!   !  int=int+sumxy(k)*(z(k)-z(k+1))*dx*dy
!   !enddo
!   ! interface total area: divide the interface layer volume per its thickness
!   ! (approximate way to calculate the interface area)
!   !area=int/(2.0d0*ch*dsqrt(2.0d0)*atanh(dabs(threshold)))

!   !! Calculate sphericity ( sphericity = A_eq/A_real )
!   !sph=aeq/area
!   !! with final formula ( sphericity = [pi(6V)^2]^(1/3) / A_real  )
!   !sph=(pi*(6.0d0*vol)**2.0d0)**(1.0d0/3.0d0) / area    
