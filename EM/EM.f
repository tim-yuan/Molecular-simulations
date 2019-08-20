
!----------------------------------------------------------------------------
! file: EM.f
! This is a fortran file 
! This file contains the function to perform energy minimization
! Written by Tim Yuan
! Sarupria Group
! Clemson University
! Project for Numerical Method class
! Created 21st April, 2019
! Last modified 21st April, 2019
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! There are some common variables
! n is the number of particles
! pos is a nx3 array storing the atom x, y, and z position
!
!
!
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
! Subroutine force and potential calculation
!----------------------------------------------------------------------------
! This is a brute force calculation
! It is pointless to implement a neighbour search in energy minimization
! A bin nearest neighbour search may be available later
!----------------------------------------------------------------------------
        subroutine potential(pos,x,y,z,rcut,eps,sig,
     &         Up,force,n)
        integer :: n,i,j
        real*8 ::  x, y, z,fx_i,fy_i,fz_i
        real*8 ::  rcut,eps,sig,upcut,dx,dy,dz,rcut2,Up
        real*8 pos(n, 3),force(n,3)
        real*8 distan(4,1)
        real*8 :: dis2,up_i
Cf2py intent(in, out, copy) Up, force
Cf2py intent(in) rcut, eps, sig
        rcut2 = rcut**2
        upcut = (4.d0*eps*((sig**12.d0)*(rcut2**(-6.d0))-(sig**6.d0)
     &       *(rcut2**(-3.d0))))
!        print*, upcut,rcut2
        !Since we are incrementing, it is very important to clear
        Up = 0.d0
        do i = 1,n
            force(i,1)=0.d0
            force(i,2)=0.d0
            force(i,3)=0.d0
        end do

        do i = 1,(n-1)
            do j = (i+1), n
!pair wise calculation for distance and applying PBC
                dx = (pos(i, 1) - pos(j,1))
                distan(1,1) = dx - x*ANINT(dx/x)
                dy = (pos(i, 2 ) - pos(j,2))
                distan(2,1) = dy - y*ANINT(dy/y)
                dz = (pos(i, 3) - pos(j,3))
                distan(3,1) = dz - z*ANINT(dz/z)
                distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &                          +distan(3,1)**2
                dis2 = distan(4,1)
!                print*,distan(1,1),distan(2,1),distan(3,1),distan(4,1)
!force and potential calculation
                if (dis2.lt.rcut2) then
                    fx_i=4d0*eps*(12d0*(sig**12d0)*(dis2**(-7d0)) -
     &              6d0*(sig**6d0)*((dis2)**(-4d0)))* distan(1,1)

                    fy_i=4d0*eps*(12d0*(sig**12d0)*(dis2**(-7d0)) -
     &              6d0*(sig**6d0)*((dis2)**(-4d0)))* distan(2,1)

                    fz_i=4d0*eps*(12d0*(sig**12d0)*(dis2**(-7d0)) -
     &              6d0*(sig**6d0)*((dis2)**(-4d0)))* distan(3,1)

                    force(i,1) = force(i,1) +fx_i
                    force(i,2) = force(i,2) +fy_i
                    force(i,3) = force(i,3) +fz_i
                    force(j,1) = force(j,1) -fx_i
                    force(j,2) = force(j,2) -fy_i
                    force(j,3) = force(j,3) -fz_i
!potential shift applied for potential calculation
                    up_i = 4.0 * eps *  ((sig**12d0)*(dis2**(-6d0))-
     &               (sig**6d0)*(dis2**(-3d0))) - upcut
!this is used for pressure calculation
!                    wx_i = fx_i*distan(1,1)
!                    wy_i = fy_i*distan(2,1)
!                    wz_i = fz_i*distan(3,1)
                    Up = Up + up_i
!                    wx = wx+wx_i
!                    wy = wy+wy_i
!                    wz = wz+wz_i
                endif
            end do
        end do
!This is very important because we only looped from 1-(n-1)!
!        w = (wx+wy+wz)/3d0
        end


