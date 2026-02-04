subroutine read_phi(nstep)
 use param
 use flowvars
 implicit none 
 integer :: nstep
 character(len=clen) :: namefile
 ! character(len=clen) :: namedir,namefile
 ! character(len=clen) :: numfile

 write(namefile, '(A,I8.8,A)') 'phi_', nstep, '.dat'
 namefile = trim(rootpath)//trim(namefile)


 ! namedir='../channel/output/'
 ! write(numfile,'(i8.8)') nstep
  
 ! !namefile=trim(namedir)//'phi_'//numfile//'.dat'
 
 open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(3) phi
 close(3,status='keep')
  
  
 return
end subroutine
    