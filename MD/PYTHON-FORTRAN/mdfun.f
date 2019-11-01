
!----------------------------------------------------------------------
! file: mdfun.f
! This is the Fortran function file for the Mesoscopic MD
! Written by Tim (Tianmu) Yuan
! This file is in compliance with the PYTHON simulation file
! This file contains varies functions for two simulations: Langevin and 
! DPD. 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!There are some commonly used variables which is explained as follow:
!
!n is the number of particle that we are simulating
!
!pos is a nx6 array, the first three columns describes the 
!x,y,z coordinates and the last three columns describes the 
!velocities in x,y,z direction 
!
!force is a nx3 array, since mass is assumed to be 1 in the entire 
!simulation. It is also the acceleration in x,y,z direction 
!respectively.
!
!
!
!----------------------------------------------------------------------

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
Cf2py intent(in, out, copy) pos
Cf2py intent(in) n, timestep 
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
!velocity in x direction, half timestep
!            pos(i, 4) = pos(i,4) + 0.5 * timestep*(force(i,1))
!velocity in y direction, half timestep
!            pos(i, 5) = pos(i,5) + 0.5 * timestep*(force(i,2))
!velocity in z direction, half timestep
!            pos(i, 6) = pos(i,6) + 0.5 * timestep*(force(i,3))
        
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
Cf2py intent(in, out, copy) pos
Cf2py intent(in) n, timestep 
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
Cf2py intent(in, out, copy) Temp, KE
Cf2py intent(in) mass, Kb, pos
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
!subroutine neighbour list
!-------------------------------------------------------------------------
!Constructing neighbour list based on neighbour list distance
!-------------------------------------------------------------------------
        subroutine neilist(pos,boxdim,rneigh,distan,list,point,n)
        implicit none
        integer :: n,i,j,nnlist
        real*8 :: rneigh
        real*8 pos(n,6), distan(4,1), boxdim(3,1)
        integer point(n),list(n*(n-1))
        real*8 :: dx,dy,dz,dis2
Cf2py intent(out, copy) list,point
Cf2py intent(in) pos,boxdim
        nnlist = 0
        do i = 1, (n-1)
            point(i) = nnlist
            do j = (i+1), n
!calculating distance and apply PBC
               dx = (pos(i, 1) - pos(j,1))
               distan(1,1) = dx - boxdim(1,1)*ANINT(dx/boxdim(1,1))
               dy = (pos(i, 2 ) - pos(j,2))
               distan(2,1) = dy - boxdim(2,1)*ANINT(dy/boxdim(2,1))
               dz = (pos(i, 3) - pos(j,3))
               distan(3,1) = dz - boxdim(3,1)*ANINT(dz/boxdim(3,1))
               distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &            +distan(3,1)**2
               dis2 = distan(4,1)
!if the distance is within the neighbour list distance
!we add the particle
               if (dis2 .LE. rneigh**2) then
                    nnlist = nnlist+1
                    list(nnlist) = j
                endif
            end do
        end do
        point(n) = nnlist
        end

!-------------------------------------------------------------------------
!subroutine neigh check 
!-------------------------------------------------------------------------
!this tells the program when to update the neighbour list
!-------------------------------------------------------------------------
!On the safe side, I used 1.2 times the displacement to check for update
!----------------------------------------------------------------------
        subroutine check(pos,posi,rskin,boxdim,update,n)
        real*8 pos(n,6),posi(n,6),distan(4,1),boxdim(3,1)
        real*8 :: dx, dy, dz,dis2
        real*8 :: rskin,dispmx
        integer :: n,update,i
Cf2py intent(out,copy) update
        dispmx = 0.0
        update = 0
        do i = 1,n
            dx = (pos(i, 1) - posi(i,1))
            distan(1,1) = dx - boxdim(1,1)*ANINT(dx/boxdim(1,1))
            dy = (pos(i, 2 ) - posi(i,2))
            distan(2,1) = dy - boxdim(2,1)*ANINT(dy/boxdim(2,1))
            dz = (pos(i, 3) - posi(i,3))
            distan(3,1) = dz - boxdim(3,1)*ANINT(dz/boxdim(3,1))
            distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &            +distan(3,1)**2
            dis2 = distan(4,1)

            dispmx = MAX(dis2,dispmx)
!            dispmx = MAX(ABS(pos(i,2)-posi(i,2)), dispmx)
!            dispmx = MAX(ABS(pos(i,3)-posi(i,3)), dispmx)
        end do
        dispmx = 1.2 * dispmx
        if (dispmx .gt. rskin**2) then
            update = 1
        else
            update = 0
        end if
        end

!-------------------------------------------------------------------------
!subroutine velocity autocorrelation function (in progress)
!-------------------------------------------------------------------------
!!This is the program for Velocity autocorrellation function
!        subroutine autocor(pos,posi,gam,Kb,Temp,
!     &  mass,stepnumber,timestep,vacor,alyacor,n)
!        real*8 pos(n,6),posi(n,6)
!        integer :: n,i
!        real*8 :: timestep,gam,Kb,Temp,mass,stepnumber
!        real*8 :: vv,vvnorm,vafxx,vafyy,vafzz,vacor,alyacor
!Cf2py intent(out,copy) vacor,alyacor
!        vv=0.0
!        vvnorm=0.0
!        do i = 1,n
!            vafxx=pos(i,4)*posi(i,4)
!            vafyy=pos(i,5)*posi(i,5)
!            vafzz=pos(i,6)*posi(i,6)
!            vv= vv+vafxx+vafyy+vafzz
!            vvnorm=vvnorm+posi(i,4)**2+posi(i,5)**2+posi(i,6)**2
!!I used a normalising factor which is the initial velocity squared
!        end do
!        vacor=vv/vvnorm
!        alyacor=(3d0*n*Kb*Temp/mass)*dexp(-gam*
!     &  dabs(stepnumber*timestep))/vvnorm
!!alyacor is the analitical solution, and normalized 
!!by initial velocity squred
!!I time 3 here because 1/2m*v**2 = 3/2Kb*T
!        end

!-------------------------------------------------------------------------
!subroutine velocity autocorrelation function (in progress)
!-------------------------------------------------------------------------
!This is the program for Velocity autocorrellation function
        subroutine autocor(pos,vel_rec, vacf_rec, vacf_norm
     &  ,current_step, n, nstep)
        real*8 ::  pos(n,6)
        real*8 ::  vel_rec(nstep, n, 3)
        real*8 ::  vacf_rec(nstep), vacf_norm(nstep)
        integer :: n, nstep, i, j, current_step
Cf2py intent(in,out, copy) vel_rec, vacf_rec, vacf_norm
Cf2py intent(in) pos, current_step, n, nstep
        do i = 1, n
            vel_rec(current_step,i,1) = pos(i, 4)
            vel_rec(current_step,i,2) = pos(i, 5)
            vel_rec(current_step,i,3) = pos(i, 6)
        end do
        do i = 1, current_step
            do j = 1, n
                vacf_rec(current_step) = vacf_rec(current_step) +
     &                            vel_rec(1,j,1)*vel_rec(i,j,1) +
     &                            vel_rec(1,j,2)*vel_rec(i,j,2) +
     &                            vel_rec(1,j,3)*vel_rec(i,j,3)
                vacf_norm(current_step) = vacf_norm(current_step) + 1
!                print*,vel_rec(i, j, 1)
            end do
        end do
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
        real*8 pos(n, 6),force(n,3), boxdim(4,1)
        real*8 distan(4, 1)
        real*8 :: dis2,fx_i,fy_i,fz_i,up_i,wx_i,wy_i,wz_i,wx,wy,wz
Cf2py intent(in, out, copy) force, w, Up
Cf2py intent(in) rcut, eps, sigma
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
                distan(1,1) = dx - boxdim(1,1)*ANINT(dx/boxdim(1,1))
                dy = (pos(i, 2) - pos(j,2))
                distan(2,1) = dy - boxdim(2,1)*ANINT(dy/boxdim(2,1))
                dz = (pos(i, 3) - pos(j,3))
                distan(3,1) = dz - boxdim(3,1)*ANINT(dz/boxdim(3,1))
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
        w = (wx+wy+wz)/3d0
        end
!
!-------------------------------------------------------------------------
!subroutine force +acc using neighbour list calculation
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!fortran2python allows you to call PYTHON function and nprandor supplies gaussian
!random number
       subroutine nlistacc(pos,distan,boxdim,rcut, eps, sigma, 
     &       neilist,point,force,w,Up,n)

        implicit none
        integer :: n,i,j,jbeg,jend,jnab
        real*8 ::  rcut,eps,sigma,upcut,rcut2, w,Up
        real*8 pos(n, 6),force(n,3),distan(4,1),boxdim(3,1) 
        real*8 :: dis2,fx_i,fy_i,fz_i,up_i,wx_i,wy_i,wz_i,wx,wy,wz
        real*8 ::  dx, dy,dz 
        integer :: point(n),neilist(n*(n-1))
Cf2py intent(in,out,copy) force,w,Up
Cf2py intent(in) pos,rcut, eps,sigma,  neilist,point
        Up=0
        w=0
        do i=1,n
            force(i,1)=0
            force(i,2)=0
            force(i,3)=0
        end do

        rcut2 = rcut**2
        upcut = (4.d0*eps*((sigma**12.d0)*(rcut2**(-6.d0))-(sigma**6.d0)
     &*(rcut2**(-3.d0))))
        wx = 0
        wy = 0
        wz = 0
        do i = 1,(n-1)
            jbeg = point(i)+1
            jend = point(i+1)

            if (jbeg .le. jend) then
                do jnab = jbeg, jend
                    j = neilist(jnab)

                    dx = (pos(i, 1) - pos(j,1))
                    distan(1,1) = dx - boxdim(1,1)*ANINT(dx/boxdim(1,1))
                    dy = (pos(i, 2) - pos(j,2))
                    distan(2,1) = dy - boxdim(2,1)*ANINT(dy/boxdim(2,1))
                    dz = (pos(i, 3) - pos(j,3))
                    distan(3,1) = dz - boxdim(3,1)*ANINT(dz/boxdim(3,1))
                    distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &                          +distan(3,1)**2
                    dis2 = distan(4,1)

                    if (dis2 .lt. rcut2) then
                        fx_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(1,1) 

                        fy_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(2,1)

                        fz_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(3,1)
                        force(i,1) = force(i,1) + fx_i
                        force(i,2) = force(i,2) + fy_i
                        force(i,3) = force(i,3) + fz_i
                        force(j,1) = force(j,1) - fx_i
                        force(j,2) = force(j,2) - fy_i
                        force(j,3) = force(j,3) - fz_i
                        up_i = 4.0*eps*((sigma**12d0)*(dis2**(-6d0))
     &-(sigma**6d0)*(dis2**(-3d0))) - upcut
                        wx_i = fx_i*distan(1,1)
                        wy_i = fy_i*distan(2,1)
                        wz_i = fz_i*distan(3,1)
                        Up = Up + up_i
                        wx = wx+wx_i
                        wy = wy+wy_i
                        wz = wz+wz_i
                    endif
                end do
            endif
        end do
        w = (wx+wy+wz)/3d0
        end
 


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!Langevin thermostat acceleration calculation using neighbour list
!-------------------------------------------------------------------------
!fortran2python allows you to call PYTHON function and nprandor supplies gaussian
!random number
       subroutine langevin(nprandor,pos,distan,boxdim,rcut, eps, sigma, 
     &       mass,gam,neilist,point,stdev,mean,force,w,Up,n)
        implicit none
        external nprandor
        real*8 :: nprandor,stdev,mean
        integer :: n,i,j,jbeg,jend,jnab
        real*8 ::  rcut,eps,sigma,upcut,rcut2, w,Up
        real*8 pos(n, 6),force(n,3),distan(4,1),boxdim(3,1)
        real*8 :: dis2,fx_i,fy_i,fz_i,up_i,wx_i,wy_i,wz_i,wx,wy,wz
        real*8 :: mass, gam, dx, dy,dz 
        integer :: point(n),neilist(n*(n-1))
Cf2py intent(in,out,copy) force,w,Up
Cf2py intent(in) pos,rcut, eps,sigma,mass, gam, neilist,point

        rcut2 = rcut*rcut
        upcut = (4.d0*eps*((sigma**12.d0)*(rcut2**(-6.d0))-(sigma**6.d0)
     &*(rcut2**(-3.d0))))
        wx = 0
        wy = 0
        wz = 0
        Up = 0
        w = 0

        do i=1,n
            force(i,1)=0
            force(i,2)=0
            force(i,3)=0
        end do

        do i = 1,(n-1)
            jbeg = point(i)+1
            jend = point(i+1)

            if (jbeg .le. jend) then
                do jnab = jbeg, jend
                    j = neilist(jnab)

                    dx = (pos(i, 1) - pos(j,1))
                    distan(1,1) = dx - boxdim(1,1)*ANINT(dx/boxdim(1,1))
                    dy = (pos(i, 2) - pos(j,2))
                    distan(2,1) = dy - boxdim(2,1)*ANINT(dy/boxdim(2,1))
                    dz = (pos(i, 3) - pos(j,3))
                    distan(3,1) = dz - boxdim(3,1)*ANINT(dz/boxdim(3,1))
                    distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &                          +distan(3,1)**2
                    dis2 = distan(4,1)

                    if (dis2 .lt. rcut2) then
                        fx_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(1,1) 

                        fy_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(2,1)

                        fz_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(3,1)
                        force(i,1) = force(i,1) + fx_i
                        force(i,2) = force(i,2) + fy_i
                        force(i,3) = force(i,3) + fz_i
                        force(j,1) = force(j,1) - fx_i
                        force(j,2) = force(j,2) - fy_i
                        force(j,3) = force(j,3) - fz_i
                        up_i = 4.0*eps*((sigma**12d0)*(dis2**(-6d0))
     &-(sigma**6d0)*(dis2**(-3d0))) - upcut
                        wx_i = fx_i*distan(1,1)
                        wy_i = fy_i*distan(2,1)
                        wz_i = fz_i*distan(3,1)
                        Up = Up + up_i
                        wx = wx+wx_i
                        wy = wy+wy_i
                        wz = wz+wz_i
                    endif
                end do
            endif
!add friction and noise to the acceleration calculated
            force(i,1)=force(i,1)-mass*gam*pos(i,4)+nprandor(mean,stdev)
            force(i,2)=force(i,2)-mass*gam*pos(i,5)+nprandor(mean,stdev)
            force(i,3)=force(i,3)-mass*gam*pos(i,6)+nprandor(mean,stdev)
        end do
!since we loop from 1 to n-1, we have to apply it to the last particle
        force(n,1)=force(n,1)-mass*gam*pos(n,4)+nprandor(mean,stdev)
        force(n,2)=force(n,2)-mass*gam*pos(n,5)+nprandor(mean,stdev)
        force(n,3)=force(n,3)-mass*gam*pos(n,6)+nprandor(mean,stdev)
        w = (wx+wy+wz)/3d0
        end
                            



!-------------------------------------------------------------------------
!subroutine DPD force+acc using neighbour list calculation
!This is built on velocity verlet
!-------------------------------------------------------------------------
!acceleration calculation for dissipative particle dynamics
        subroutine dpdacc(nprandor,pos,distan,boxdim,rcut, eps, sigma, 
     &       gam,Kb,Tset,timestep,stdev,mean,neilist,point,force,w,Up,n)

        implicit none
        external nprandor
        integer :: n,i,j,jbeg,jend,jnab
        real*8 ::  rcut,eps,sigma,upcut,rcut2,w,Up
        real*8 pos(n, 6),force(n,3),distan(4,1),boxdim(3,1),dis
        real*8 :: dis2,fx_i,fy_i,fz_i,up_i,wx_i,wy_i,wz_i,wx,wy,wz
        real*8 :: gam,Kb,Tset,dx, dy,dz,timestep,nprandor,mean
        integer :: point(n),neilist(n*(n-1))
        real*8 :: vijx,vijy,vijz,dotprod,fr,fRx,fRy,fRz,stdev
        real*8 :: fDx,fDy,fDz,fd
Cf2py intent(in,out,copy) force,w,Up
Cf2py intent(in) pos,rcut, eps,sigma, gam, neilist,point
        rcut2 = rcut*rcut
        upcut = (4.d0*eps*((sigma**12.d0)*(rcut2**(-6.d0))-(sigma**6.d0)
     &*(rcut2**(-3.d0))))
        wx = 0
        wy = 0
        wz = 0
        Up = 0
        w=0

        do i=1,n
            force(i,1)=0
            force(i,2)=0
            force(i,3)=0
        end do

        do i = 1,(n-1)
            jbeg = point(i)+1
            jend = point(i+1)

            if (jbeg .le. jend) then
                do jnab = jbeg, jend
                    j = neilist(jnab)

                    dx = (pos(i, 1) - pos(j,1))
                    distan(1,1) = dx - boxdim(1,1)*ANINT(dx/boxdim(1,1))
                    dy = (pos(i, 2) - pos(j,2))
                    distan(2,1) = dy - boxdim(2,1)*ANINT(dy/boxdim(2,1))
                    dz = (pos(i, 3) - pos(j,3))
                    distan(3,1) = dz - boxdim(3,1)*ANINT(dz/boxdim(3,1))
                    distan(4,1) = distan(1,1)**2+distan(2,1)**2
     &                          +distan(3,1)**2
                    dis2 = distan(4,1)
                    dis = dsqrt(dis2)
!                    wdr = 2d0*(1-dis/rcut)
!                    wrr = dsqrt(wdr)
                    if (dis2 .lt. rcut2) then
                        fx_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(1,1) 

                        fy_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(2,1)

                        fz_i=4d0*eps*(12d0*(sigma**12d0)*(dis2**(-7d0))
     &-6d0*(sigma**6d0)*((dis2)**(-4d0)))*distan(3,1)
!                        force(i,1) = force(i,1) + fx_i
!                        force(i,2) = force(i,2) + fy_i
!                        force(i,3) = force(i,3) + fz_i
!                        force(j,1) = force(j,1) - fx_i
!                        force(j,2) = force(j,2) - fy_i
!                        force(j,3) = force(j,3) - fz_i
                        up_i = 4d0*eps*((sigma**12d0)*(dis2**(-6d0))
     &-(sigma**6d0)*(dis2**(-3d0))) - upcut
                        wx_i = fx_i*distan(1,1)
                        wy_i = fy_i*distan(2,1)
                        wz_i = fz_i*distan(3,1)
                        Up = Up + up_i
                        wx = wx+wx_i
                        wy = wy+wy_i
                        wz = wz+wz_i


!!dissipative force component
                        vijx = pos(i,4)-pos(j,4)
                        vijy = pos(i,5)-pos(j,5)
                        vijz = pos(i,6)-pos(j,6)
                        dotprod=(vijx*distan(1,1)+vijy*distan(2,1)+
     &  vijz*distan(3,1))
                        fd = -gam*dotprod/dis2
                        fDx=fd*distan(1,1)
                        fDy=fd*distan(2,1)
                        fDz=fd*distan(3,1)
!                        print*,nprandor(mean,stdev) 
                        fr = dsqrt(2d0*gam*Tset*Kb)/dis/dsqrt(timestep)
!                        print*,dsqrt(2d0*gam*Tset*Kb),dsqrt(timestep)
                        fRx=fr*nprandor(mean,stdev)*distan(1,1)
                        fRy=fr*nprandor(mean,stdev)*distan(2,1)
                        fRz=fr*nprandor(mean,stdev)*distan(3,1)
!                        print*, fDx, fDy, fDz,fRx,fRy,fRz
                        force(i,1)=force(i,1)+
     &  fx_i+fDx+fRx
                        force(i,2)=force(i,2)+
     &  fy_i+fDy+fRy
                        force(i,3)=force(i,3)+
     &  fz_i+fDz+fRz
                        force(j,1)=force(j,1)-
     &  fx_i-fDx-fRx
                        force(j,2)=force(j,2)-
     &  fy_i-fDy-fRy
                        force(j,3)=force(j,3)-
     &  fz_i-fDz-fRz


                    endif
                end do
            endif
        end do
        w = (wx+wy+wz)/3d0
        end

!-------------------------------------------------------------------------
!subroutine 
!-------------------------------------------------------------------------
!This is the program for berendsen thermostat
        subroutine berendsen(pos,mass,Kb,tset,dt,tau,Temp,n)
        real*8 pos(n,6)
        integer :: n,i
        real*8 :: dt, Kb,tset,mass
        real*8 :: sumv2,sumv2i,Tn,sclfactor, Temp,tau
        real*8 :: vafxx,vafyy,vafzz
Cf2py intent(in,out,copy) pos,Temp
        sumv2=0d0
        sumv2i=0d0
!        print*,mass,Kb,tset,dt,tau,Temp
        do i = 1,n
            vafxx=pos(i,4)*pos(i,4)
            vafyy=pos(i,5)*pos(i,5)
            vafzz=pos(i,6)*pos(i,6)
            sumv2=sumv2+vafxx+vafyy+vafzz
        end do
        Tn = mass * sumv2/(3d0*n*Kb)
        sclfactor=sqrt(1.0+((dt)/tau)*((tset/Tn)-1.0))
        do i = 1,n
            pos(i,4) = pos(i,4) * sclfactor
            pos(i,5) = pos(i,5) * sclfactor
            pos(i,6) = pos(i,6) * sclfactor
            vafxx=pos(i,4)*pos(i,4)
            vafyy=pos(i,5)*pos(i,5)
            vafzz=pos(i,6)*pos(i,6)
            sumv2i=sumv2i+vafxx+vafyy+vafzz
        end do
        Temp=mass * sumv2i/(3d0*n*Kb)
        end

