# Tim Yuan
# This code performs MD simulation
# It is capable of simulating canonical and microcanonical ensemble

code ="MD.py"
modified ="16th May, 2018"

import sys
import os.path
import numpy as np
import argparse
import pylab
import math
import scipy.special
import datetime
import time
np.set_printoptions(threshold=sys.maxsize)



def main():
#-------------------------------------------------------------------------
#Define some parameters
#-------------------------------------------------------------------------
    start = time.time()
    args = parseargs()
    f_in = args.inputfile
    f_out = args.outputfile
    eps = args.eps
    sigma = args.sig
    dt = args.timestep
    nstep = args.nstep
    f_ener = args.energyfile
    Tset = args.temp

#    pos_vec=[]
    mass=1.0
    Kb=1.0
    rcut = 3*sigma
    rcut2 = rcut**2
    rskin = sigma
    rneigh = rcut+rskin
    tau = 0.1

#-------------------------------------------------------------------------
#Simulation run
#-------------------------------------------------------------------------
    #initialize the position and velocity
    pos=[]
    #Read the input file and store necessary informations
    pos,n,boxdim=readfile(pos, f_in)
    #calculate the volume
    V=boxdim[0]*boxdim[1]*boxdim[2]
    #check the atoms number
    assert (n==len(pos)),\
        "Error, number of atoms in gro file not consistent"

    #open files to write the steps and energy
    fxyz=open(f_out,'w')
    fener=open(f_ener,'w')
    fener.write("{:<25s}{:<25s}{:<25s}{:<13s}{:<25s}{:<10s}\n".\
        format("Potential E", "Kinetic E", "Total E", "Tempearture","Pressure","Nlist update"))
    if ((Tset-0)<10**(-6)):
        fener.write("#The simulation is run under microcanonical ensemble\n")
    else:
        fener.write("#The simulation is run under canonical ensemble with Berendsen thermostat\n")
        fener.write("#The target temperature is "+str(Tset)+"\n")

    #loop over the steps
    n_list = neigh(pos,n, rneigh,boxdim)
    pos_prev=np.copy(pos)
    nlist_update=0
    for i in range(nstep):
        #Neighbour list check
        dispmax=0
        for k in range(n):
            d_x, d_y, d_z, dis2 = dist(pos[k], pos_prev[k], boxdim)
            dispmax=max(dis2,dispmax)
        if dispmax>(rskin**2):
            n_list = neigh(pos,n, rneigh,boxdim)
            pos_prev=np.copy(pos)
            nlist_update+=1
        #Update the velocities and positions using leapfrog
        pos, Up, KE, Tem, w = leapfrog(n,  pos, dt, boxdim, rcut, eps, sigma,n_list,mass,Kb) 

        #-------------------------------------------------------------------------
        # if the temperature is 0, the system is run under microcanonical ensemble
        # otherwise, the temperature will be controlled using Berendsen thermostat
        #-------------------------------------------------------------------------
        if ((Tset-0)>10**(-6)) and (i%10 ==0):
            pos, Tem=berendsen(pos,mass,Kb,Tset,n,dt,tau)

        if (i%100 ==0):
            print("There are "+str(nstep-i)+" steps remaining")

        #Write data to files
        fxyz.write(str(n) + '\n' + '\n')
        for k in range(int(n)):
            pos[k][0] = pos[k][0] - math.floor(pos[k][0]/boxdim[0])*boxdim[0]
            pos[k][1] = pos[k][1] - math.floor(pos[k][1]/boxdim[1])*boxdim[1]
            pos[k][2] = pos[k][2] - math.floor(pos[k][2]/boxdim[2])*boxdim[2]
            fxyz.write('LJP' + ' ' + str(pos[k][0]) + ' ' +str(pos[k][1]) + ' ' + str(pos[k][2]) + '\n')
        #Pressure calculation
        Pcc =  16*3*math.pi*(n/V)*eps*sigma**3 *(2/3*(sigma/rcut)**9 - (sigma/rcut)**3) #pressure correction
        P = (n*Kb*Tem +w)/(V) + Pcc #pressure calculated with correction
        ET = Up + KE
        fener.write("{:<25.7f}{:<25.7f}{:<25.7f}{:<13.2f}{:<25.7f}{:<10d}\n".format(Up, KE, ET, Tem,P,nlist_update))


    #close the file handle
    fxyz.close()
    fener.close()
    end = time.time()
    print("The runing time is "+str(end - start)+" s")

#-------------------------------------------------------------------------
#Define necessary functions to run the simulation
#-------------------------------------------------------------------------
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Tim Yuan\nSarupria Research Group\n"+modified+"\n\n"
                        "Perform energy minimization"
                        )
    parser.add_argument("-i", "--inputfile", help="input gro file",
                        default="input.gro", metavar='')
    parser.add_argument("-o", "--outputfile", help="specify the name of the output file",
                        default="config.xyz", metavar='')
    parser.add_argument("-energy", "--energyfile", help="specify the name of the energy (output) file",
                        default="energy.xvg", metavar='')
    parser.add_argument("-eps", "--eps", help="epsilon value for the LJ particle",
                        default=1,type=float, metavar='')
    parser.add_argument("-sig", "--sig", help="sigma value for the LJ particle",
                        default=1,type=float, metavar='')
    parser.add_argument("-dt", "--timestep", help="Time step for the MD simulation",
                        default=0.01,type=float, metavar='')
    parser.add_argument("-temp", "--temp", help="Temperature for the MD simulation"\
                            "0 would mean the simulation is under microcanonical ensemble",
                        default=0,type=float, metavar='')
    parser.add_argument("-nstep", "--nstep", help="Number of steps for the simulation",
                        default=1000,type=int, metavar='')

    args = parser.parse_args()
    return args

#-------------------------------------------------------------------------
#Now lets get into the MD functions
#-------------------------------------------------------------------------
def readfile(pos,f_in):
    if os.path.isfile(f_in):
        f_in_hand = open(f_in, 'r')
        grofile=f_in_hand.readlines()
        n = int(grofile[1].strip())
        for line in grofile[2:(len(grofile)-1)]:
            ele = [line[0:5].strip(), line[5:10].strip(),line[10:15].strip(),\
                line[15:20].strip(),line[20:28].strip(),line[28:36].strip(),\
                line[36:44].strip(),line[44:52].strip(),line[52:60].strip(),\
                line[60:68].strip()]
            pos_vec=[float(ele[4]),float(ele[5]),float(ele[6]),float(ele[7]),float(ele[8]),float(ele[9])] 
            pos.append(pos_vec)
        x=float(grofile[-1].split()[0])
        y=float(grofile[-1].split()[1])
        z=float(grofile[-1].split()[2])
        boxdim=[x,y,z]
    else:
        print("No file found, please input the correct name of the gro file")
        exit()
    pos=np.asarray(pos)
    return pos, n,boxdim
#-------------------------------------------------------------------------
#Periodic boundary conditions
#-------------------------------------------------------------------------
def dist(x, y, boxdim):
    d_x = (x[0] - y[0])
    d_x = d_x - boxdim[0]*round(d_x/boxdim[0])
    d_y = (x[1] - y[1])
    d_y =  d_y - boxdim[1]*round(d_y/boxdim[1])
    d_z = (x[2] - y[2])
    d_z =  d_z - boxdim[2]*round(d_z/boxdim[2])
    return d_x, d_y, d_z, ((d_x)**2 + (d_y)**2 + (d_z)**2)

#-------------------------------------------------------------------------
#neighbour list
#-------------------------------------------------------------------------
def neigh(pos,n,rneigh,boxdim):
    n_list =[]
    for k in range(int(n)-1):
        n_list.append([])
    for i in range(int(n)-1):
        for j in range(i+1, int(n)):
            x = pos[i]
            y = pos[j]
            d_x, d_y, d_z, dis2 = dist(x, y, boxdim)
            if dis2 <= rneigh**2 :
                n_list[i].append(j)
    return n_list

#-------------------------------------------------------------------------
#Acceleration calculated by force using L-J potential. F = -dU/dr, Force in x direction: F_x = -dU/dr *dr/dx
#-------------------------------------------------------------------------
def acc(n, rcut, eps, sigma, pos, n_list,boxdim):
    rcut2=rcut*rcut
    Upcut = ( 4.0 * eps *  ((sigma**12)*(rcut2**(-6))- (sigma**6)*(rcut2**(-3))))
    fx = np.zeros(int(n))
    fy = np.zeros(int(n))
    fz = np.zeros(int(n))
    Up = 0.0
    wxx = 0.0
    wyy = 0.0
    wzz = 0.0
    dis3 = []
    for i in range(int(n)-1):
        for j in range(len(n_list[i])):
            d_x, d_y, d_z, dis2 = dist(pos[i], pos[n_list[i][j]], boxdim)
            if dis2< rcut2:
                fx_i =4.0 * eps * (12*(sigma**12)*(dis2**(-7)) - 6*(sigma**6)*(dis2)**(-4))* d_x
                fy_i =4.0 * eps * (12*(sigma**12)*(dis2**(-7)) - 6*(sigma**6)*(dis2)**(-4))* d_y
                fz_i =4.0 * eps * (12*(sigma**12)*(dis2**(-7)) - 6*(sigma**6)*(dis2)**(-4))* d_z
                Up_i = 4.0 * eps *  ((sigma**12)*(dis2**(-6))- (sigma**6)*(dis2**(-3))) - Upcut
                Up =Up + Up_i
                fx[i] = fx[i] + fx_i
                fy[i] = fy[i] + fy_i
                fz[i] = fz[i] + fz_i
                fx[n_list[i][j]] = fx[n_list[i][j]] - fx_i
                fy[n_list[i][j]] = fy[n_list[i][j]] - fy_i
                fz[n_list[i][j]] = fz[n_list[i][j]] - fz_i
                wx_i = fx_i *d_x            #calculating virial for pressure
                wy_i = fy_i *d_y
                wz_i = fz_i *d_z
                wxx = wx_i +wxx
                wyy = wy_i +wyy
                wzz = wz_i +wzz
    w = 1/3 *(wxx+wyy+wzz)
    return fx, fy, fz, Up, w

#-------------------------------------------------------------------------
#Leap frog integrator
#-------------------------------------------------------------------------
def leapfrog(n,  pos, dt, boxdim,rcut, eps, sigma, n_list,mass,Kb):
    ax, ay, az, Up,w = acc(n, rcut, eps, sigma, pos, n_list, boxdim)
#    dt = timestep
    dt2 = dt * 0.5
    acc0 = {}
    accn = {}
    for i in range(int(n)):
        acc0[i] = [ax[i], ay[i], az[i]]
        pos[i][3] = pos[i][3] + acc0[i][0] * dt2
        pos[i][4] = pos[i][4] + acc0[i][1] * dt2
        pos[i][5] = pos[i][5] + acc0[i][2] * dt2
        pos[i][0] = pos[i][0] + pos[i][3] *dt
        pos[i][0] = pos[i][0] - math.floor(pos[i][0]/boxdim[0])*boxdim[0]
        pos[i][1] = pos[i][1] + pos[i][4] *dt
        pos[i][1] = pos[i][1] - math.floor(pos[i][1]/boxdim[1])*boxdim[1]
        pos[i][2] = pos[i][2] + pos[i][5] *dt
        pos[i][2] = pos[i][2] - math.floor(pos[i][2]/boxdim[2])*boxdim[2]
    axx, ayy, azz, Up,w = acc(n, rcut, eps, sigma, pos, n_list, boxdim)
    KE = 0
    for i in range(int(n)):
        accn[i] =  [axx[i], ayy[i], azz[i]]
        pos[i][3] = pos[i][3] + accn[i][0]*dt2
        pos[i][4] = pos[i][4] + accn[i][1]*dt2
        pos[i][5] = pos[i][5] + accn[i][2]*dt2
        KE = KE + (pos[i][3]**2 + pos[i][4]**2 +pos[i][5]**2)
    KE = 0.5 * KE
    Tem = 2.0 * mass* KE/(3*n*Kb)
    return  pos, Up, KE, Tem, w

#-------------------------------------------------------------------------
#Thermostat, Velocity scale
#-------------------------------------------------------------------------
def thermostat(pos,mass,Kb,n):
    sumv2 = 0
    sumv22 = 0
    for i in range(int(n)):
        sumv2 = sumv2 + ((pos[i][3])**2 + (pos[i][4])**2 + (pos[i][5])**2)
    Tn = mass * sumv2/(3.0 * n * Kb)
    lamdaa = math.sqrt(Tset/Tn)
    for i in range(int(n)):
        pos[i][3] = pos[i][3] * lamdaa
        pos[i][4] = pos[i][4] * lamdaa
        pos[i][5] = pos[i][5] * lamdaa
        sumv22 = sumv22 + ((pos[i][3])**2 + (pos[i][4])**2 + (pos[i][5])**2)
    Tnn = mass * sumv22/(3.0 * n * Kb)
    return pos ,Tn

#----------------------------------------------------
#Thermostat, Berendsen
#-------------------------------------------------------------------------
def berendsen(pos,mass,Kb,Tset,n,dt,tau):
    sumv2 = 0
    sumv22 = 0
    for i in range(n):
        sumv2 = sumv2 +((pos[i][3])**2 + (pos[i][4])**2 + (pos[i][5])**2)
    Tn = mass * sumv2/(3.0 * n * Kb)
#    print(mass, sumv2, n, Kb,dt)
    sclfactor = np.sqrt(1.0+(dt/tau)*((Tset/Tn)-1.0))
    for i in range(int(n)):
        pos[i][3] = pos[i][3] * sclfactor
        pos[i][4] = pos[i][4] * sclfactor
        pos[i][5] = pos[i][5] * sclfactor
        sumv22 = sumv22 + ((pos[i][3])**2 + (pos[i][4])**2 + (pos[i][5])**2)
    Tnn = mass * sumv22/(3.0 * n * Kb)
    return pos, Tnn





if __name__ == '__main__':
    main()
