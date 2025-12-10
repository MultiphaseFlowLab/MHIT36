subroutine read_fields(nstep)

use commondata
integer :: nstep
character(len=40) :: namedir,namefile
character(len=8) :: numfile
character(len=3) :: setnum
logical :: check



namedir='../../channel/output/'
write(numfile,'(i8.8)') nstep


allocate(u(nx,ny,nz))
allocate(v(nx,ny,nz))
allocate(w(nx,ny,nz))
allocate(phi(nx,ny,nz))
allocate(theta(nx,ny,nz))


write(*,*) 'Reading step ',nstep,' out of ',nend,' , flow'
 !reading u
 if (uflag .eq. 1) then
    namefile=trim(namedir)//'u_'//numfile//'.dat'
    open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(666) u
    close(666,status='keep')
 endif
 !reading v
 if (vflag .eq.	1) then
   namefile=trim(namedir)//'v_'//numfile//'.dat'
   open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(667) v
   close(667,status='keep')
endif
!reading w
   if (wflag .eq.	1) then
   namefile=trim(namedir)//'w_'//numfile//'.dat'
   open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(668) w
   close(668,status='keep')
endif
!reading phi
   if (phiflag .eq.	1) then
   namefile=trim(namedir)//'phi_'//numfile//'.dat'
   open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(668) phi
   close(668,status='keep')
endif
!reading theta
   if (thetaflag .eq.	1) then
   namefile=trim(namedir)//'theta_'//numfile//'.dat'
   open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(668) theta
   close(668,status='keep')
endif

write(*,*) "max phi", maxval(phi)

! generate paraview output file
call generate_output(nstep)


deallocate(u)
deallocate(v)
deallocate(w)
deallocate(phi)
deallocate(theta)

return
end
