
!----------------------------------------------------------------------
! MD nve using Fortran
! Written by Tim Yuan
!----------------------------------------------------------------------

        program main
        implicit none
        integer ::    i,j,k,n,t
        integer ::    counti, nstep
        real*8  ::    pos(500,6)
        real*8  ::    box(3)
        real*8  ::    force(500,3)
        real*8  ::    force2(500,3)
        real*8  ::    sigma, eps, rcut,w,Up
        real*8  ::    temp, ke, ET,mass,kb,dt
        sigma = 1d0
        eps = 1d0
        force = 0d0
        pos = 0d0
        box = 0d0
        rcut = 3.5*sigma
        w = 0
        Up = 0
        nstep = 1000
        temp = 0d0
        ke = 0d0
        ET = 0d0
        mass = 1d0
        kb = 1d0
        dt = 0.001
        !--------------------------------------------------------------
        ! Create box
        !--------------------------------------------------------------
        box(1)= 8d0
        box(2)= 8d0
        box(3)= 8d0
        n = 500
        counti = 1
        do while (counti .le. n) 
          do i = 1,int(box(1))
            do j = 1,int(box(2))
              do k = 1,int(box(3))
                if (counti .gt. n) then
                  exit
                end if
                pos(counti, 1) = i
                pos(counti, 2) = j
                pos(counti, 3) = k
                counti = counti + 1
              end do
            end do
          end do
        end do
        pos(500,1) = 0.5
        pos(500,2) = 0.5
        pos(500,3) = 0.5


        !--------------------------------------------------------------
        ! run the MD code
        !--------------------------------------------------------------
        call run(pos,box,rcut,eps,sigma,dt,mass,kb,nstep,n)

        end program



!----------------------------------------------------------------------
!subroutine run and write
!----------------------------------------------------------------------
!This is the first part of the velocity verlet integerator
!----------------------------------------------------------------------
        subroutine run(pos,box,rcut,eps,sigma,dt,mass,kb,nstep,n)
        integer ::    i,j,k,n,t,nstep
        real*8  ::    pos(n,6)
        real*8  ::    box(3)
        real*8  ::    force(n,3)
        real*8  ::    rcut,eps,sigma,dt,mass,kb

        real*8  ::    distan(4,1)
        real*8  ::    w, up, temp,ke
        real*8  ::    force2(n,3)

        distan = 0d0

        print*, "MD calculation starting"
        open(1,file='lj.xyz',action='write',status='replace')
        open(11,file='lj.xvg',action='write',status='replace')
        write(11,22) "#1.step 2.time 3.Up 4.KE 5.ET 6.T"
22      format(a60)
 
        call acc(pos,distan,box,rcut, eps, sigma,force,w, Up,n)
        do t = 1,nstep
!          call kinetic(pos, mass, kb, temp, ke, n)
!          et = up+ke
          call vv1(pos,dt,force,n)
          force2 = force

          call kinetic(pos, mass, kb, temp, ke, n)
          et = up+ke
 
          call acc(pos,distan,box,rcut, eps, sigma,force,w, Up,n)

          call vv2(pos,dt,force,force2,n)
!          call kinetic(pos, mass, kb, temp, ke, n)
!          et = up+ke
          

          write(1,2) n
2         format(I15.0)
          write(1,2) 
          do i = 1,n
!            print*, i,pos(i,1), pos(i,2),pos(i,3),force(i,1),force(i,2),
!     &force(i,3)
            write(1,3) 'LJP',pos(i,1), pos(i,2),pos(i,3)
3         format(a5,f17.5,f17.5,f17.5)
          end do

          write(11,33) t,t*dt,Up,Ke,Et,temp
33        format (i10,f10.3,e18.9,e18.9,e18.9,e18.9)


        end do
        print*, "MD calculation finished, fininshing output"
        close(1)
        close(11)
        end


!----------------------------------------------------------------------
!subroutine velocity verlet Integrator1
!----------------------------------------------------------------------
!This is the first part of the velocity verlet integerator
!----------------------------------------------------------------------
        subroutine vv1(pos, timestep ,force,n)
        implicit none
        integer :: n,i
        real*8 ::  timestep
        real*8 pos(n, 6)
        real*8 force(n,3)
        do i = 1,n
!position in x direction
            pos(i, 1) = pos(i,1) + timestep * pos(i, 4) + 0.5 * 
     &  timestep**2* force(i,1)
!position in y direction
            pos(i, 2) = pos(i,2) + timestep * pos(i, 5) + 0.5 * 
     &  timestep**2* force(i,2)
!position in z direction
            pos(i, 3) = pos(i,3) + timestep * pos(i, 6) + 0.5 * 
     &  timestep**2* force(i,3)
!            print*, pos(i,1), pos(i,2),pos(i,3),force(i,1),force(i,2),
!     &force(i,3)
           
        end do
        end
!----------------------------------------------------------------------
!subroutine velocity verlet Integrator
!----------------------------------------------------------------------
!
        subroutine vv2(pos,timestep,force,force2, n)
        integer :: n,i
        real*8 ::  timestep
        real*8 pos(n, 6),force(n, 3),force2(n,3)
        do i = 1,n
!velocity in x direction, timestep
            pos(i,4)=pos(i,4)+0.5*timestep*(force(i,1)+force2(i,1))
!velocity in y direction, timestep
            pos(i,5)=pos(i,5)+0.5*timestep*(force(i,2)+force2(i,2))
!velocity in z direction, timestep
            pos(i,6)=pos(i,6)+0.5*timestep*(force(i,3)+force2(i,3))
        end do
        end

!-------------------------------------------------------------------------
!subroutine Kinetic Energy and Temperature calculation
!-------------------------------------------------------------------------
!Temperature calculation 
!-------------------------------------------------------------------------
        subroutine kinetic(pos, mass, Kb, Temp, KE, n)
        integer :: i, n
        real*8 ::  mass, Kb, KE_i, Temp, KE
        real*8 pos(n, 6)
        Temp=0
        KE=0
        do i = 1,n
            KE_i = ((pos(i,4))**2 +(pos(i,5))**2+(pos(i,6))**2)/2
            KE = KE+KE_i
        end do
        KE = KE*mass
        Temp = 2d0 *KE/(3d0*n*Kb)
        end

!-------------------------------------------------------------------------
!subroutine force +acc calculation
!-------------------------------------------------------------------------
!This is a brute force calculation for all the particles
!It is very insufficient for large number of particles, but it is better
!for small numbers comparing with neighbour list calculation
!-------------------------------------------------------------------------
        subroutine acc(pos,distan,boxdim,rcut, eps, sigma,
     &        force,w, Up,n)
        integer :: n,i,j
        real*8 ::  rcut,eps,sigma,upcut,dx,dy,dz,rcut2, w,Up
        real*8 pos(n, 6),force(n,3), boxdim(3)
        real*8 distan(4, 1)
        real*8 :: dis2,fx_i,fy_i,fz_i,up_i,wx_i,wy_i,wz_i,wx,wy,wz
        force = 0d0
        Up = 0d0
        w = 0d0
        rcut2 = rcut**2
        upcut = (4.d0*eps*((sigma**12.d0)*(rcut2**(-6.d0))-(sigma**6.d0)
     &       *(rcut2**(-3.d0))))
        wx = 0
        wy = 0
        wz = 0
        do i = 1,(n-1)
            do j = (i+1), n
!pair wise calculation for distance and applying PBC
                dx = (pos(i, 1) - pos(j,1))
                distan(1,1) = dx - boxdim(1)*ANINT(dx/boxdim(1))
                dy = (pos(i, 2) - pos(j,2))
                distan(2,1) = dy - boxdim(2)*ANINT(dy/boxdim(2))
                dz = (pos(i, 3) - pos(j,3))
                distan(3,1) = dz - boxdim(3)*ANINT(dz/boxdim(3))
                distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &                          +distan(3,1)**2
                dis2 = distan(4,1)
!acceleration calculation for parwise
                if (dis2.lt.rcut2) then
                    fx_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0)) - 
     &              6d0*(sigma**6d0)*((dis2)**(-4d0)))* distan(1,1)

                    fy_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0)) - 
     &              6d0*(sigma**6d0)*((dis2)**(-4d0)))* distan(2,1)

                    fz_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0)) - 
     &              6d0*(sigma**6d0)*((dis2)**(-4d0)))* distan(3,1)

                    force(i,1) = force(i,1) +fx_i
                    force(i,2) = force(i,2) +fy_i
                    force(i,3) = force(i,3) +fz_i
                    force(j,1) = force(j,1) -fx_i
                    force(j,2) = force(j,2) -fy_i
                    force(j,3) = force(j,3) -fz_i
!       print*, pos(i,1),pos(i,2),pos(i,3),pos(j,1),pos(j,2),pos(j,3),
!    &distan(1,1),distan(2,1) ,distan(3,1),fx_i,fy_i,fz_i

!potential shift applied for potential calculation
                    up_i = 4.0 * eps *  ((sigma**12d0)*(dis2**(-6d0))-
     &               (sigma**6d0)*(dis2**(-3d0))) - upcut
!this is used for pressure calculation
                    wx_i = fx_i*distan(1,1)
                    wy_i = fy_i*distan(2,1)
                    wz_i = fz_i*distan(3,1)
                    Up = Up + up_i
                    wx = wx+wx_i
                    wy = wy+wy_i
                    wz = wz+wz_i
                endif
            end do
        end do
!This is very important because we only looped from 1-(n-1)!
        w = -(wx+wy+wz)/2d0
        end

