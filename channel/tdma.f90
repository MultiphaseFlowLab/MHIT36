! Subroutine still under development, not included in the makefile at the moment (7/8/25)
subroutine tdmau
    implicit none
    integer il, jl
    ! TDMA for u velocity, solve from 1 to nz+1 where 1 is the w collocation point on the bottom and nz+1 the upper one
    ! Load all variables via modules 
    ! move pencil information in a module
    ! Set up tridiagonal system for each i and j
    ! The system is: (A_z) * w(k-1,j,i) + (B_k) * pc(k,ky,kx) + (C_k) * pc(k+1,j,i) = rhs(k,j,i)
    ! Dirilecht BC: w = 0 at k=1 and k=nz+1
    ! a, b, c, d are complex, they can be optimized defined them as real, however this may introduce additional data movement, better to use always the same a, b, c, d? 
    
    use velocity
    implicit none
    call c_f_pointer(c_devloc(vel_d), vel3d, piZ%shape)

   !$acc parallel loop collapse(2) gang private(a, b, c, d, factor, sol) 
   do jl = 1, npy
      do il = 1, npx

         !inner nodes
         do k = 2, nz-1
            a(k) =  1.0d0*dzi*dzi
            b(k) = -2.0d0*dzi*dzi 
            c(k) =  1.0d0*dzi*dzi
            d(k) =  u(k,il,jl) ! check if correct and consisten with CN
         enddo

         ! Dirilecht BC at k=1
         a(1) =  0.0d0
         b(1) =  1.0d0
         c(1) =  0.0d0
         d(1) =  0.0d0

         ! Dirilecht BC at k=nz+1
         a(nz) =  0.0d0
         b(nz) =  1.0d0
         c(nz) =  0.0d0
         d(nz) =  0.0d0

         ! Forward elimination (Thomas)
         !$acc loop seq
         do k = 2, nz
            factor = a(k)/b(k-1)
            b(k) = b(k) - factor*c(k-1)
            d(k) = d(k) - factor*d(k-1)
         end do

         ! Back substitution
         sol(nz) = d(nz)/b(nz+1)
         ! check on pivot like flutas?
         !$acc loop seq
         do k = nz-1, 1, -1
            sol(k) = (d(k) - c(k)*sol(k+1))/b(k)
         end do

         ! Store solution in array that do the back FFT
         do k=1,nz
            vel3d(k,il,jl) = sol(k)
         enddo 
      end do
   end do
end subroutine












subroutine tdmav
    implicit none
    integer il, jl
    ! TDMA for v velocity, solve from 1 to nz where 1 is the v collocation point on the bottom and nz the upper one
    ! Load all variables via modules 
    ! move pencil information in a module
    ! Set up tridiagonal system for each i and j
    ! The system is: (A_z) * w(k-1,j,i) + (B_k) * pc(k,ky,kx) + (C_k) * pc(k+1,j,i) = rhs(k,j,i)
    ! Dirilecht BC: v = 0 at k=1 and k=nz
    ! a, b, c, d are complex, they can be optimized defined them as real, however this may introduce additional data movement, better to use always the same a, b, c, d? 
    ! Same as u but just RHS is different
    
    use velocity
    implicit none
    call c_f_pointer(c_devloc(vel_d), vel3d, piZ%shape)

   !$acc parallel loop collapse(2) gang private(a, b, c, d, factor, sol) 
   do jl = 1, npy
      do il = 1, npx

         !inner nodes
         do k = 2, nz-1
            a(k) =  1.0d0*dzi*dzi
            b(k) = -2.0d0*dzi*dzi 
            c(k) =  1.0d0*dzi*dzi
            d(k) =  v(k,il,jl) ! check if correct and consisten with CN
         enddo

         ! Dirilecht BC at k=1
         a(1) =  0.0d0
         b(1) =  1.0d0
         c(1) =  0.0d0
         d(1) =  0.0d0

         ! Dirilecht BC at k=nz
         a(nz+1) =  0.0d0
         b(nz+1) =  1.0d0
         c(nz+1) =  0.0d0
         d(nz+1) =  0.0d0

         ! Forward elimination (Thomas)
         !$acc loop seq
         do k = 2, nz
            factor = a(k)/b(k-1)
            b(k) = b(k) - factor*c(k-1)
            d(k) = d(k) - factor*d(k-1)
         end do

         ! Back substitution
         sol(nz) = d(nz)/b(nz)
         ! check on pivot like flutas?
         !$acc loop seq
         do k = nz-1, 1, -1
            sol(k) = (d(k) - c(k)*sol(k+1))/b(k)
         end do

         ! Store solution in array that do the back FFT
         do k=1,nz
            vel3d(k,il,jl) = sol(k)
         enddo 
      end do
   end do
end subroutine
















subroutine tdmaw
    implicit none
    integer il, jl
    ! TDMA for w velocity, solve from 1 to ny+1 where 1 is the w collocation point on the bottom and ny+1 the upper one
    ! Load all variables via modules 
    ! move pencil information in a module
    ! Set up tridiagonal system for each i and j
    ! The system is: (A_z) * w(k-1,j,i) + (B_k) * pc(k,ky,kx) + (C_k) * pc(k+1,j,i) = rhs(k,j,i)
    ! Dirilecht BC: w = 0 at k=1 and k=nz+1
    ! a, b, c, d are complex, they can be optimized defined them as real, however this may introduce additional data movement, better to use always the same a, b, c, d? 
    
    use velocity
    implicit none
    call c_f_pointer(c_devloc(vel_d), vel3d, piZ%shape)

   !$acc parallel loop collapse(2) gang private(a, b, c, d, factor, sol) 
   do jl = 1, npy
      do il = 1, npx

         !inner nodes
         do k = 2, nz
            a(k) =  1.0d0*dzi*dzi
            b(k) = -2.0d0*dzi*dzi 
            c(k) =  1.0d0*dzi*dzi
            d(k) =  w(k,il,jl) ! check if correct and consisten with CN
         enddo

         ! Dirilecht BC at k=1
         a(1) =  0.0d0
         b(1) =  1.0d0
         c(1) =  0.0d0
         d(1) =  0.0d0

         ! Dirilecht BC at k=nz+1
         a(nz+1) =  0.0d0
         b(nz+1) =  1.0d0
         c(nz+1) =  0.0d0
         d(nz+1) =  0.0d0

         ! Forward elimination (Thomas)
         !$acc loop seq
         do k = 2, nz+1
            factor = a(k)/b(k-1)
            b(k) = b(k) - factor*c(k-1)
            d(k) = d(k) - factor*d(k-1)
         end do

         ! Back substitution
         sol(nz+1) = d(nz+1)/b(nz+1)
         ! check on pivot like flutas?
         !$acc loop seq
         do k = nz, 1, -1
            sol(k) = (d(k) - c(k)*sol(k+1))/b(k)
         end do

         ! Store solution in array that do the back FFT
         ! only up to nz, because i am not trasposing the halo, i can eventually
         do k=1,nz
            vel3d(k,il,jl) = sol(k)
         enddo 
      end do
   end do
end subroutine








subroutine tdmap
    ! TDMA for pressure, solve from 0 to ny+1 where 0 is the first ghost node and ny+1 the upper ghost node
    ! When ready, remove it from the main code and just call this function
    ! Load all variables via modules 
    ! move pencil information in a module
    ! Set up tridiagonal system for each i and j
    ! The system is: (A_z) * pc(k-1,ky,kx) + (B_k) * pc(k,ky,kx) + (C_k) * pc(k,ky,kx) = rhs(k,ky,kx)
    ! Neumann BC: d/dz pc = 0 at w collocation points
    ! Fill diagonals and rhs for each
    ! 0 and ny+1 are the ghost nodes
    ! a, b, c and are complex (because of FFT)
    
    use velocity

    implicit none

   offsets(piZ_d2z%order(1)) = piZ_d2z%lo(1) - 1
   offsets(piZ_d2z%order(2)) = piZ_d2z%lo(2) - 1
   offsets(piZ_d2z%order(3)) = piZ_d2z%lo(3) - 1
   call c_f_pointer(c_devloc(psi_d), psi3d, piZ_d2z%shape)

   xoff = offsets(1)
   yoff = offsets(2)
   npx = np(1)
   npy = np(2)

   !$acc parallel loop collapse(2) gang private(a, b, c, d, factor, sol) 
   do jl = 1, npy
      do il = 1, npx
         ! compute index global wavenumber ig and jg
         jg = yoff + jl
         ig = xoff + il

         do k = 1, nz
            a(k) =  1.0d0*dzi*dzi
            b(k) = -2.0d0*dzi*dzi - kx_d(ig)*kx_d(ig) - ky_d(jg)*ky_d(jg)
            c(k) =  1.0d0*dzi*dzi
            d(k) =  psi3d(k,il,jl)
         enddo

         ! Neumann BC at bottom
         a(0) =  0.0d0
         b(0) = -1.0d0*dzi*dzi - kx_d(ig)*kx_d(ig) - ky_d(jg)*ky_d(jg)
         c(0) =  1.0d0*dzi*dzi
         d(0) =  0.0d0

         ! ghost node elimintaion trick
         ! Neumann BC at top
         a(nz+1) =  1.0d0*dzi*dzi
         b(nz+1) = -1.0d0*dzi*dzi - kx_d(ig)*kx_d(ig) - ky_d(jg)*ky_d(jg)
         c(nz+1) =  0.0d0
         d(nz+1) =  0.0d0

         ! Enforce pressure at one point? one interior point, avodig messing up with BC
         ! need brackets?
         if (ig == 1 .and. jg == 1) then
            a(1) = 0.d0
            b(1) = 1.d0
            c(1) = 0.d0
            d(1) = 0.d0
         end if

         ! Forward elimination (Thomas)
         !$acc loop seq
         do k = 1, nz+1
            factor = a(k)/b(k-1)
            b(k) = b(k) - factor*c(k-1)
            d(k) = d(k) - factor*d(k-1)
         end do

         ! Back substitution
         sol(nz+1) = d(nz+1)/b(nz+1)
         ! check on pivot like flutas?
         !$acc loop seq
         do k = nz, 0, -1
            sol(k) = (d(k) - c(k)*sol(k+1))/b(k)
         end do

         ! Store solution in array that do the back FFT
         do k=1,nz
            psi3d(k,il,jl) = sol(k)
         enddo 
      end do
   end do


end subroutine