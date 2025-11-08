subroutine generate_output(nstep)

use commondata

integer :: nstep,numx,numy,numz
integer :: i,k,j
character(len=40) :: namefile
character(len=80) :: buffer
character(len=10) :: lf
character(len=8) :: str1,str2,str3
character(len=16) :: str4

! end of line character
lf=achar(10)

! fields included
nfields=uflag + vflag + wflag + phiflag

!input??
x_start=1
x_end=nx
dnx=3
y_start=1
y_end=ny
dny=3
z_start=1
z_end=nz
dnz=3

numx=0
numy=0
numz=0
do i=x_start,x_end,dnx
 numx=numx+1
enddo
do j=y_start,y_end,dny
 numy=numy+1
enddo
do k=z_start,z_end,dnz
 numz=numz+1
enddo



write(namefile,'(a,i8.8,a)') './output/OUTPAR_',nstep,'.vtk'
open(66,file=trim(namefile),status='new',form='unformatted',access='stream',convert='big_endian')
! start writing vtk file
! write header
buffer='# vtk DataFile Version 3.0'//lf
write(66) trim(buffer)
buffer='Phase field'//lf
write(66) trim(buffer)
buffer='BINARY'//lf
write(66) trim(buffer)
buffer='DATASET RECTILINEAR_GRID'//lf
 write(66) trim(buffer)

 !write grid
 write(str1(1:8),'(i8)') numx
 write(str2(1:8),'(i8)') numy
 write(str3(1:8),'(i8)') numz
 buffer='DIMENSIONS '//str1//str2//str3//lf
 write(66) trim(buffer)
 buffer='X_COORDINATES '//str1//'  float'//lf ;
 write(66) trim(buffer)
 do i=x_start,x_end,dnx
  write(66) real(x(i))
 enddo
 buffer='Y_COORDINATES '//str2//'  float'//lf ;
 write(66) trim(buffer)
 do j=y_start,y_end,dny
  write(66) real(y(j))
 enddo
 buffer='Z_COORDINATES '//str3//'  float'//lf ;
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  write(66) real(z(k))
 enddo

 ! write content (data format)
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer='POINT_DATA '//str4//lf
 write(66) trim(buffer)
 write(str1(1:8),'(i8)') nfields
 buffer='FIELD FieldData '//str1//lf
 write(66) trim(buffer)


! write u field
if (uflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'u 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(u(i,j,k))
            enddo
        enddo
    enddo
endif

! write v field
if (vflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'v 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(v(i,j,k))
            enddo
        enddo
    enddo
endif

! write w field
if (wflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'w 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(w(i,j,k))
           enddo
        enddo
    enddo
endif

! write phi field
if (phiflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'phi 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(phi(i,j,k))
           enddo
        enddo
    enddo
endif

return
end
