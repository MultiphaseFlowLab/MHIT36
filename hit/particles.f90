module particles
 use sort
 implicit none 
 ! Total number of particles
 integer:: npart
 ! Local (per process) number of particles
 integer :: nploc
 ! Estimated maximum number of particles per process
 integer :: nplocmax 
 ! Array containing the info of the particles (every column is a different particle)
 double precision, allocatable :: part(:,:)
 ! Buffer for particle array
 double precision, allocatable :: partbuff(:,:) 
 ! Number of Info per particle
 integer, parameter :: ninfop = 12 
 !....Indices for particle info
 integer, parameter :: Itag=1                     ! Row 1: Tag of the particle
 integer, parameter :: Ixp=2,Iyp=3,Izp=4          ! Row 2 to 4: Current x,y,z of the particle 
 integer, parameter :: Iup=5,Ivp=6,Iwp=7          ! Row 5 to 7: Current u,v,w of the particle
 integer, parameter :: Iphip=8                    ! Row 8: phi at the current location of the particle
 integer, parameter :: Iup1=9,Ivp1=10,Iwp1=11    ! Row 9 to 11: Previous (step) u,v,w of the particle
 integer, parameter :: Iphip1=12                 ! Row 12: Previous (step) u,v,w of the particle
 
 ! Array for particle ordering
 double precision, allocatable :: vec_p(:)
 integer, allocatable :: order_p(:) 
 integer :: inpart

 ! Number of Particles to be Transferred
 integer::Nsendp1,Nsendm1,Nrecvp1,Nrecvm1

 ! Buffer for particle transfer
 double precision, allocatable,dimension(:,:)::buffvar1,buffvar2
 
 contains
!***********************************************************************
subroutine particlegenerator(inflag)
 use param
 use mpivar
 use cudecompvar
 implicit none
 integer :: inflag
 integer :: i
   
 !....Particle Tag
 do i = 1, nploc
 part(i,Itag) = dble(i)
 part(i,Itag) = dble(i) + rank * nploc
 enddo
    

 if(inflag.eq.1)then
     !....Generate xp position
     call RANDOM_SEED()
     call RANDOM_NUMBER(vec_p)
     do i = 1, nploc
     part(i,Ixp) = vec_p(i) * lx
     enddo
     
     !....Generate yp position
     call RANDOM_SEED()
     call RANDOM_NUMBER(vec_p)
     do i = 1, nploc
     part(i,Iyp) = vec_p(i) * lyloc + yinf 
     end do
     
     !....Generate zp position
     call RANDOM_SEED()
     call RANDOM_NUMBER(vec_p)
     do i = 1, nploc
     part(i,Izp) = vec_p(i) * lzloc + zinf 
     end do
 endif

 return
end subroutine
!***********************************************************************
subroutine linearinterpolation
    use velocity
    use param
    use cudecompvar
    implicit none
    integer :: i 
    double precision :: xpt, ypt, zpt
    integer :: ix, iy, iz, ixplus1
    double precision :: tx, ty, tz

 !$acc parallel loop present(part,u,v,w) private(xpt,ypt,zpt,ix,iy,iz,tx,ty,tz)
  do i = 1, nploc
    ! Particle position
    xpt = part(i,Ixp)
    ypt = part(i,Iyp)
    zpt = part(i,Izp)

    ! Global left-cell indices (uniform, fully periodic)
    ix = floor( xpt / dx )
    iy = floor( ypt / dx )
    iz = floor( zpt / dx )

    ! Local fractions in [0,1)
    tx = xpt/dx - dble(ix)
    ty = ypt/dx - dble(iy)
    tz = zpt/dx - dble(iz)

    ! Map to 1-based global indices and wrap, then to local by offsets
    ix = 1 + ix
    ixplus1 = ix + 1
    ! Periodicity
    if(ixplus1>nx) ixplus1 = 1
    ! Periodicity is already accounted by the halos
    iy = 1 + iy - pix_yoff
    iz = 1 + iz - pix_zoff

    ! 1-D linear per component (MAC: faces along native dir)
    part(i,Iup) = (1.0d0 - tx)*u(ix    , iy    , iz    ) + tx*u(ixplus1  , iy    , iz    )
    part(i,Ivp) = (1.0d0 - ty)*v(ix    , iy    , iz    ) + ty*v(ix       , iy+1  , iz    )
    part(i,Iwp) = (1.0d0 - tz)*w(ix    , iy    , iz    ) + tz*w(ix       , iy    , iz+1  )
  end do
 !$acc end parallel loop
    
 return
end subroutine
!***********************************************************************   
subroutine SortPartY()
 implicit none
 integer :: i, j

 !$acc parallel loop collapse(2)
 do j = 1, ninfop
    do i = 1, nploc
        partbuff(i,j) = part(i,j)
        vec_p(i) = part(i,Iyp)
        order_p(i) = i
    enddo
 enddo
    
 !$acc host_data use_device(vec_p,order_p)
 call fsort(vec_p,order_p,nploc,.true.)
 !$acc end host_data

 !$acc parallel loop collapse(2)
 do j = 1, ninfop
    do i = 1, nploc
        part(i,j) = partbuff(order_p(i),j)
    end do
 end do

 return
end subroutine
!***********************************************************************
subroutine SortPartZ()
 implicit none
 integer :: i, j

 !$acc parallel loop collapse(2)
 do j = 1, ninfop
    do i = 1, nploc
        partbuff(i,j) = part(i,j)
        vec_p(i) = part(i,Izp)
        order_p(i) = i
    enddo
 enddo
    
 !$acc host_data use_device(vec_p,order_p)
 call fsort(vec_p,order_p,nploc,.true.)
 !$acc end host_data

 !$acc parallel loop collapse(2)
 do j = 1, ninfop
    do i = 1, nploc
        part(i,j) = partbuff(order_p(i),j)
    end do
 end do

 return
end subroutine
!***********************************************************************
subroutine CountPartTransfY()
 use param
 implicit none 
 integer :: flag,scanner
 integer :: i
 integer :: addp,addm
 
 Nsendp1=0 !Number of particle to send to nid+1 (in Y)
 Nrecvm1=0 !Number of particle to receive from nid-1 (in Y)
 Nsendm1=0 !Number of particle to send to nid-1 (in Y)
 Nrecvp1=0 !Number of particle to receive from nid+1 (in Y)


 !$acc parallel loop default(present) private(addp,addm) reduction(+:Nsendp1,Nsendm1)
 do i = 1, nploc
  addp=0
  addm=0
  if(part(i,Izp).GE.ysup) addp=1
  if(part(i,Izp).LT.yinf) addm=1
  Nsendp1=Nsendp1+addp
  Nsendm1=Nsendm1+addm
 end do

 return
end subroutine 
!***********************************************************************
subroutine CountPartTransfZ()
 use param
 implicit none 
 integer :: flag,scanner
 integer :: i
 integer :: addp,addm
 
 Nsendp1=0 !Number of particle to send to nid+1 (in Z)
 Nrecvm1=0 !Number of particle to receive from nid-1 (in Z)
 Nsendm1=0 !Number of particle to send to nid-1 (in Z)
 Nrecvp1=0 !Number of particle to receive from nid+1 (in Z)


 !$acc parallel loop default(present) private(addp,addm) reduction(+:Nsendp1,Nsendm1)
 do i = 1, nploc
  addp=0
  addm=0
  if(part(i,Izp).GE.zsup) addp=1
  if(part(i,Izp).LT.zinf) addm=1
  Nsendp1=Nsendp1+addp
  Nsendm1=Nsendm1+addm
 end do

 return
end subroutine 
!***********************************************************************
subroutine SendPartUP(direction)
  use mpi
  use mpivar
  implicit none
  integer :: direction
  integer::NsendMAX
  integer,dimension(MPI_STATUS_SIZE):: status
  integer :: ii, jj
  integer :: nidp1, nidm1
  
  if(direction.eq.2)then
      nidp1 = nidp1y
      nidm1 = nidm1y
  elseif(direction.eq.3)then
      nidp1 = nidp1z
      nidm1 = nidm1z
  endif
  
  !...Proc nid send to proc nid+1 the number of particles 
  !...that nid+1 has to receive from nid
  call MPI_Sendrecv(Nsendp1,1,MPI_INTEGER, nidp1, nidp1,&
  Nrecvm1,1,MPI_INTEGER,nidm1,rank,MPI_COMM_WORLD,status,ierr)
  
  call MPI_Allreduce(Nsendp1,NsendMAX,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  if(NsendMAX.GT.0)then
      if(Nsendp1.GT.0)then
         !$acc parallel loop collapse(2)
         do jj=1,Nsendp1
          do ii=1,ninfop
           buffvar1(ii,jj)=part(nploc-Nsendp1+jj,ii)
          end do
         end do
      end if
    
      !$acc host_data use_device(buffvar1,buffvar2)
      !....Proc nid send to nid+1 and receive from nid-1
      call MPI_Sendrecv(buffvar1(1,1),int(ninfop*Nsendp1),MPI_DOUBLE_PRECISION,&
      nidp1, nidp1,buffvar2(1,1),int(ninfop*Nrecvm1),MPI_DOUBLE_PRECISION,&
      nidm1,rank,MPI_COMM_WORLD,STATUS,ierr)
      !$acc end host_data
    
      if(Nrecvm1.GT.0)then
         !....Each proc place received data in the first block
         !$acc parallel loop collapse(2)
         do jj=1,ninfop
          do ii=1,Nrecvm1
           partbuff(ii,jj)=buffvar2(jj,ii)
          end do
         end do
      end if
  end if
  
  
  !....Adjust central block (not transferred data)
  !$acc parallel loop collapse(2)
  do jj=1,ninfop
   do ii=1,(nploc-Nsendm1-Nsendp1)
    partbuff(Nrecvm1+ii,jj)=part(Nsendm1+ii,jj)
   end do
  end do
  
  return
end subroutine 
!***********************************************************************
subroutine SendPartDOWN(direction)
  use mpi
  use mpivar
  implicit none
  integer :: direction
  integer::NsendMAX
  integer,dimension(MPI_STATUS_SIZE):: status
  integer :: ii, jj
  integer :: nidp1, nidm1

  
  if(direction.eq.2)then
      nidp1 = nidp1y
      nidm1 = nidm1y
  elseif(direction.eq.3)then
      nidp1 = nidp1z
      nidm1 = nidm1z
  endif
  
  !...Proc nid send to proc nid-1 the number of particles 
  !...that nid-1 has to receive from nid
  call MPI_Sendrecv(Nsendm1,1,MPI_INTEGER, nidm1, nidm1,&
  Nrecvp1,1,MPI_INTEGER,nidp1,rank,MPI_COMM_WORLD,STATUS,ierr)
  
  call MPI_Allreduce(Nsendm1,NsendMAX,1,MPI_INTEGER,MPI_MAX,MPI_Comm_World,ierr)
  
  
  if(NsendMAX.GT.0)then
    if(Nsendm1.GT.0)then
       !$acc parallel loop collapse(2)
       do jj=1,Nsendm1
        do ii=1,ninfop
         buffvar1(ii,jj)=part(jj,ii)
        end do
       end do
    end if
  
    !$acc host_data use_device(buffvar1,buffvar2)
    !....Proc nid send to nid-1 and receive from nid+1
    call MPI_Sendrecv(buffvar1(1,1),int(ninfop*Nsendm1),MPI_DOUBLE_PRECISION,&
    nidm1, nidm1,buffvar2(1,1),int(ninfop*Nrecvp1),MPI_DOUBLE_PRECISION,&
    nidp1, rank ,MPI_COMM_WORLD,STATUS,ierr)
    !$acc end host_data
  
    if(Nrecvp1.GT.0)then
       !....Each proc append received data
       !$acc parallel loop collapse(2)
       do jj=1,ninfop
        do ii=1,Nrecvp1
         partbuff(ii+Nrecvm1+nploc-Nsendm1-Nsendp1,jj)=buffvar2(jj,ii)
        end do
       end do
    end if
  
  end if
  
  !...Update Local Number of Paricles
  nploc=nploc+Nrecvm1+Nrecvp1-Nsendp1-Nsendm1
  !$acc parallel loop collapse(2)
  do jj=1,ninfop
   do ii=1,nploc
    part(ii,jj)=partbuff(ii,jj)
   end do
  end do
  
  
  return
end subroutine 
!***********************************************************************
subroutine ParticlesLeakage()
  use param
  !....Particles which escape from the domain along x, y and z must 
  !....be reintroduced
  implicit none
  
  !$acc kernels
  !....Along x
  part(1:nploc,Ixp)=part(1:nploc,Ixp)-lx*FLOOR(part(1:nploc,Ixp)/lx)
  !....Along y
  part(1:nploc,Iyp)=part(1:nploc,Iyp)-lx*FLOOR(part(1:nploc,Iyp)/lx)
  !....Along z
  part(1:nploc,Izp)=part(1:nploc,Izp)-lx*FLOOR(part(1:nploc,Izp)/lx)
  !$acc end kernels
    
  return
end subroutine
!***********************************************************************
end module particles