recursive subroutine flood_fill(top,s_drop,id,jd,kd)
 use param
 implicit none
 integer :: top(nx,nx,nx),s_drop(nx,nx,nx)
 integer :: id,jd,kd
 ! integer :: inew,jnew,knew
 integer :: idp,jdp,kdp
 integer :: idm,jdm,kdm


 s_drop(id,jd,kd)=1

 ! Indices increased by one (with periodicity)
 idp = mod(id, nx) + 1      ! id + 1 with periodicity
 jdp = mod(jd, nx) + 1      ! jd + 1 with periodicity
 kdp = mod(kd, nx) + 1      ! kd + 1 with periodicity

 ! Indices reduced by one (with periodicity)
 idm = mod(nx + id - 2, nx) + 1  ! id - 1 with periodicity
 jdm = mod(nx + jd - 2, nx) + 1  ! jd - 1 with periodicity
 kdm = mod(nx + kd - 2, nx) + 1  ! kd - 1 with periodicity

 ! Face-connected (search only in +-x, +-y and +-z direction 6 directions)
 !+x
 if((top(idp,jd,kd).eq.1).and.(s_drop(idp,jd,kd).eq.0)) call flood_fill(top,s_drop,idp,jd,kd)
 !-x
 if((top(idm,jd,kd).eq.1).and.(s_drop(idm,jd,kd).eq.0)) call flood_fill(top,s_drop,idm,jd,kd)

 !+y
 if((top(id,jdp,kd).eq.1).and.(s_drop(id,jdp,kd).eq.0)) call flood_fill(top,s_drop,id,jdp,kd)
 !-y
 if((top(id,jdm,kd).eq.1).and.(s_drop(id,jdm,kd).eq.0)) call flood_fill(top,s_drop,id,jdm,kd)

 !+z 
 if((top(id,jd,kdp).eq.1).and.(s_drop(id,jd,kdp).eq.0)) call flood_fill(top,s_drop,id,jd,kdp)
 !-z
 if((top(id,jd,kdm).eq.1).and.(s_drop(id,jd,kdm).eq.0)) call flood_fill(top,s_drop,id,jd,kdm)

 ! Edge-Connected Combinations (12 edges)
 !    Δx | Δy | Δz
 !1   +1 | +1 | 0
 !2   +1 | -1 | 0
 !3   -1 | +1 | 0
 !4   -1 | -1 | 0
 !5   +1 | 0 | +1
 !6   +1 | 0 | -1
 !7   -1 | 0 | +1
 !8   -1 | 0 | -1
 !9    0 | +1 | +1
 !10   0 | +1 | -1
 !11   0 | -1 | +1
 !12   0 | -1 | -1

 !1   +1 | +1 | 0
 if((top(idp,jdp,kd).eq.1) .and. (s_drop(idp,jdp,kd).eq.0)) call flood_fill(top,s_drop,idp,jdp,kd)

 !2   +1 | -1 | 0
 if((top(idp,jdm,kd).eq.1) .and. (s_drop(idp,jdm,kd).eq.0)) call flood_fill(top,s_drop,idp,jdm,kd)

 !3   -1 | +1 | 0
 if((top(idm,jdp,kd).eq.1) .and. (s_drop(idm,jdp,kd).eq.0)) call flood_fill(top,s_drop,idm,jdp,kd)

 !4   -1 | -1 | 0
 if((top(idm,jdm,kd).eq.1) .and. (s_drop(idm,jdm,kd).eq.0)) call flood_fill(top,s_drop,idm,jdm,kd)
 
 !5   +1 | 0 | +1
 if((top(idp,jd,kdp).eq.1) .and. (s_drop(idp,jd,kdp).eq.0)) call flood_fill(top,s_drop,idp,jd,kdp)
 
 !6   +1 | 0 | -1
 if((top(idp,jd,kdm).eq.1) .and. (s_drop(idp,jd,kdm).eq.0)) call flood_fill(top,s_drop,idp,jd,kdm)
 
 !7   -1 | 0 | +1
 if((top(idm,jd,kdp).eq.1) .and. (s_drop(idm,jd,kdp).eq.0)) call flood_fill(top,s_drop,idm,jd,kdp)
 
 !8   -1 | 0 | -1
 if((top(idm,jd,kdm).eq.1) .and. (s_drop(idm,jd,kdm).eq.0)) call flood_fill(top,s_drop,idm,jd,kdm)
 
 !9    0 | +1 | +1
 if((top(id,jdp,kdp).eq.1) .and. (s_drop(id,jdp,kdp).eq.0)) call flood_fill(top,s_drop,id,jdp,kdp)
 
 !10   0 | +1 | -1
 if((top(id,jdp,kdm).eq.1) .and. (s_drop(id,jdp,kdm).eq.0)) call flood_fill(top,s_drop,id,jdp,kdm)
 
 !11   0 | -1 | +1
 if((top(id,jdm,kdp).eq.1) .and. (s_drop(id,jdm,kdp).eq.0)) call flood_fill(top,s_drop,id,jdm,kdp)
 
 !12   0 | -1 | -1
 if((top(id,jdm,kdm).eq.1) .and. (s_drop(id,jdm,kdm).eq.0)) call flood_fill(top,s_drop,id,jdm,kdm)

 ! Corner-Connected Combinations (8 corners)
 !  | Δx | Δy | Δz
 !1 | +1 | +1 | +1
 !2 | +1 | +1 | -1
 !3 | +1 | -1 | +1
 !4 | +1 | -1 | -1
 !5 | -1 | +1 | +1
 !6 | -1 | +1 | -1
 !7 | -1 | -1 | +1
 !8 | -1 | -1 | -1

 !1  +1 | +1 | +1
 if ((top(idp,jdp,kdp).eq.1) .and. (s_drop(idp,jdp,kdp).eq.0)) call flood_fill(top,s_drop,idp,jdp,kdp)
 
 !2  +1 | +1 | -1
 if ((top(idp,jdp,kdm).eq.1) .and. (s_drop(idp,jdp,kdm).eq.0)) call flood_fill(top,s_drop,idp,jdp,kdm)
 
 !3  +1 | -1 | +1
 if ((top(idp,jdm,kdp).eq.1) .and. (s_drop(idp,jdm,kdp).eq.0)) call flood_fill(top,s_drop,idp,jdm,kdp)
 
 !4  +1 | -1 | -1
 if ((top(idp,jdm,kdm).eq.1) .and. (s_drop(idp,jdm,kdm).eq.0)) call flood_fill(top,s_drop,idp,jdm,kdm)
 
 !5  -1 | +1 | +1
 if ((top(idm,jdp,kdp).eq.1) .and. (s_drop(idm,jdp,kdp).eq.0)) call flood_fill(top,s_drop,idm,jdp,kdp)
 
 !6  -1 | +1 | -1
 if ((top(idm,jdp,kdm).eq.1) .and. (s_drop(idm,jdp,kdm).eq.0)) call flood_fill(top,s_drop,idm,jdp,kdm)
 
 !7  -1 | -1 | +1
 if ((top(idm,jdm,kdp).eq.1) .and. (s_drop(idm,jdm,kdp).eq.0)) call flood_fill(top,s_drop,idm,jdm,kdp)
 
 !8  -1 | -1 | -1
 if ((top(idm,jdm,kdm).eq.1) .and. (s_drop(idm,jdm,kdm).eq.0)) call flood_fill(top,s_drop,idm,jdm,kdm)

 
  return
end subroutine
!*********************************************************
!OLD
  ! ! search on diagonals (20 directions)
  ! ! ! search in x+1,y+1,[z-1,z,z+1]  --> 3 cases
  ! ! inew=mod(id+1,nx)+1
  ! ! jnew=mod(jd+1,nx)+1
  ! ! if((top(inew,jnew,kd).eq.1).and.(s_drop(inew,jnew,kd).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
  ! ! if(kd+1.le.nx)then
  ! !  if((top(inew,jnew,kd+1).eq.1).and.(s_drop(inew,jnew,kd+1).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
  ! ! endif
  ! ! if(kd-1.ge.1)then
  ! !  if((top(inew,jnew,kd-1).eq.1).and.(s_drop(inew,jnew,kd-1).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
  ! ! endif

  ! !....Search on diagonals (20 directions)
  ! ! search in x+1, y+1, [z-1, z, z+1]  --> 3 directions
  ! inew = mod(id, nx) + 1      ! id + 1 with periodicity
  ! jnew = mod(jd, nx) + 1      ! jd + 1 with periodicity
  
  ! ! z (centered)
  ! if ((top(inew, jnew, kd) == 1) .and. (s_drop(inew, jnew, kd) == 0)) then
  !   call flood_fill(top, s_drop, inew, jnew, kd)
  ! endif
  
  ! ! z+1 with periodicity
  ! knew = mod(kd, nx) + 1
  ! if ((top(inew, jnew, knew) == 1) .and. (s_drop(inew, jnew, knew) == 0)) then
  !   call flood_fill(top, s_drop, inew, jnew, knew)
  ! endif
  
  ! ! z-1 with periodicity
  ! knew = mod(nx + kd - 2, nx) + 1
  ! if ((top(inew, jnew, knew) == 1) .and. (s_drop(inew, jnew, knew) == 0)) then
  !   call flood_fill(top, s_drop, inew, jnew, knew)
  ! endif

    
  ! ! search in x+1,y-1,[z-1,z,z+1]  --> 3 cases and in x+1,y,[z-1,z+1]  --> 2 cases
  ! inew = mod(id, nx) + 1      ! id + 1 with periodicity
  ! jnew = mod(nx+jd-2,nx) + 1    ! jd - 1 with periodicity
  ! ! z (centered)
  ! if ((top(inew, jnew, kd) == 1) .and. (s_drop(inew, jnew, kd) == 0)) call flood_fill(top, s_drop, inew, jnew, kd)
  ! endif

  ! ! z+1 with periodicity
  ! knew = mod(kd, nx) + 1
  ! if ((top(inew, jnew, knew) == 1) .and. (s_drop(inew, jnew, knew) == 0)) call flood_fill(top, s_drop, inew, jnew, knew)!y-1
  ! if ((top(inew, jd  , knew) == 1) .and. (s_drop(inew, jd  , knew) == 0)) call flood_fill(top, s_drop, inew, jd  , knew)!y

  ! ! z-1 with periodicity
  ! knew = mod(nx + kd - 2, nx) + 1
  ! if ((top(inew, jnew, knew) == 1) .and. (s_drop(inew, jnew, knew) == 0)) call flood_fill(top, s_drop, inew, jnew, knew)!y-1
  ! if ((top(inew, jd  , knew) == 1) .and. (s_drop(inew, jd  , knew) == 0)) call flood_fill(top, s_drop, inew, jd  , knew)!y


  ! ! ! if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
  ! ! if(kd+1.le.nz)then
  ! !  if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
  ! !  if((top(inew,kd+1,jd).eq.1).and.(s_drop(inew,kd+1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd+1)
  ! ! endif
  ! ! if(kd-1.ge.1)then
  ! !  if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
  ! !  if((top(inew,kd-1,jd).eq.1).and.(s_drop(inew,kd-1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd-1)
  ! ! endif

  ! ! search in x-1,y+1,[z-1,z,z+1]  --> 3 cases
  ! inew = mod(nx+id-2,nx) + 1    ! id - 1 with periodicity
  ! jnew = mod(jd, nx) + 1      ! jd + 1 with periodicity


  ! ! search in x-1,y+1,[z-1,z,z+1]  --> 3 cases
  ! inew=mod(nx+id-1,nx)+1
  ! jnew=mod(jd+1,ny)+1
  ! if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
  ! if(kd+1.le.nz)then
  !  if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
  ! endif
  ! if(kd-1.ge.1)then
  !  if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
  ! endif
  
  ! ! search in x-1,y-1,[z-1,z,z+1]  --> 3 cases and in x-1,y,[z-1,z+1]  --> 2 cases
  ! jnew=mod(ny+jd-1,ny)+1
  ! if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
  ! if(kd+1.le.nz)then
  !  if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
  !  if((top(inew,kd+1,jd).eq.1).and.(s_drop(inew,kd+1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd+1)
  ! endif
  ! if(kd-1.ge.1)then
  !  if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
  !  if((top(inew,kd-1,jd).eq.1).and.(s_drop(inew,kd-1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd-1)
  ! endif
  
  ! ! search in x,[y-1,y+1],[z-1,z+1]  --> 4 cases
  ! jnew=mod(jd+1,ny)+1
  ! if(kd+1.le.nz)then
  !  if((top(id,kd+1,jnew).eq.1).and.(s_drop(id,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd+1)
  ! endif
  ! if(kd-1.ge.1)then
  !  if((top(id,kd-1,jnew).eq.1).and.(s_drop(id,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd-1)
  ! endif
  
  ! jnew=mod(ny+jd-1,ny)+1
  ! if(kd+1.le.nz)then
  !  if((top(id,kd+1,jnew).eq.1).and.(s_drop(id,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd+1)
  ! endif
  ! if(kd-1.ge.1)then
  !  if((top(id,kd-1,jnew).eq.1).and.(s_drop(id,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd-1)
  ! endif 
!END OLD