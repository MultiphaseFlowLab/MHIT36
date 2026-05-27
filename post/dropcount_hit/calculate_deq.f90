subroutine calculate_deq(vol_cells, nstep)
  use param
  implicit none
  integer, intent(in)  :: vol_cells, nstep
  character(len=clen)  :: time, namefile
  double precision     :: vol, deq

  vol = dble(vol_cells)

  ! Equivalent diameter: twice the radius of the equivalent sphere (V = (4/3)*pi*R^3)
  deq = dx * (n6opi * vol)**(1.0d0/3.0d0)

  write(*,'(1x,a,E16.6)') 'Diameter', deq
  write(*,*)

  write(time,'(i8.8)') nstep
  namefile = 'deq_'//trim(time)//'.dat'
  open(3, file=trim(namefile), access='append', form='formatted', status='unknown')
    write(3,'(E12.5)') deq
  close(3, status='keep')

end subroutine calculate_deq
