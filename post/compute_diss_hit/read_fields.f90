subroutine read_fields(nstep)

use commondata

integer :: nstep

character(len=40) :: namedir,namefile
character(len=8) :: numfile


namedir='../../hit/output/'
write(numfile,'(i8.8)') nstep

namefile=trim(namedir)//'u_'//numfile//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) u
close(666,status='keep')

namefile=trim(namedir)//'v_'//numfile//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) v
close(666,status='keep')

namefile=trim(namedir)//'w_'//numfile//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) w
close(666,status='keep')

return
end




