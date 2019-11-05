# Tim Yuan
# This code performs MD simulation
# It is capable of simulating canonical and microcanonical ensemble
# thermostats including Berendsen, Langevin, and dissipative particle dynamics are available

code ="MD.py"
modified ="11th July, 2018"

import sys
import os.path
import numpy as np
import argparse
import pylab
import math
import scipy.special
import datetime
import random
import mdfun
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
    tset = args.tset
    thermo=args.thermo
    vacf = args.vacffile
#    pos_vec=[]
    mass=1.0
    Kb=1.0
    rcut = 3.5*sigma
    rcut2 = rcut**2
    rskin = sigma
    rneigh = rcut+rskin
    mean = 0
    distan = np.zeros((4), order='F')
    KE = 0
    Temp = 0
    Tem = 0
    w=0
    Up=0
    #for berendsen thermostat
    tau = 0.1          #coupling constant for Berendsen thermostat
    #for langevin dynamics
    gam = 1.0           #gamma value for langevin thermostat
    stdev=np.sqrt(2.0*gam*mass*Kb*tset/(dt))
    stdevdpd=1.0
    
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
    #initialize the force array
    force = np.zeros((int(n), 3),order='F')

    #open files to write the steps and energy
    fxyz=open(f_out,'w')
    fener=open(f_ener,'w')
    fener.write("{:<25s}{:<25s}{:<25s}{:<13s}{:<25s}{:<10s}\n".format("Potential E",\
             "Kinetic E", "Total E (kJ/mol)", "Tempearture (K)","Pressure (bar)","Nlist update"))
    if vacf != None:
        vel_rec = np.zeros((int(nstep),int(n),3 ), order = 'F')
        vacf_rec = np.zeros(int(nstep), order = 'F')
        vacf_norm = np.zeros(int(nstep), order = 'F')

    nlist_update=0
    neilist, point = mdfun.neilist(pos, boxdim, rneigh,distan)
    pos_prev=np.copy(pos)

    #-------------------------------------------------------------------------
    # for NVE or NVT with Berendsen thermostat
    #-------------------------------------------------------------------------
#    if (thermo != "dpd") and (thermo !="la"):
    if (tset-0<10**(-6)):
        force,w,Up=mdfun.nlistacc(pos,distan,boxdim,rcut,eps,sigma,neilist,point,force,w,Up)
        fener.write("#The simulation is run at microcanonical ensemble\n")
    else:
        if thermo=="br":
            force,w,Up=mdfun.nlistacc(pos,distan,boxdim,rcut,eps,sigma,neilist,point,force,w,Up)
            fener.write("#The simulation is run at canonical ensemble\n")
            fener.write("#Berendsen thermostat is use, target temperature is "+str(tset)+"\n")
        elif thermo=="la":
            force,w,Up=mdfun.langevin(nprandor,pos,distan,boxdim,rcut,eps,sigma, mass, gam, \
                neilist,point,stdev,mean,force,w,Up)
            fener.write("#The simulation is run at canonical ensemble\n")
            fener.write("#Langevin thermostat is use, target temperature is "+str(tset)+"\n")
        elif thermo=="dpd":
            force,w,Up=mdfun.dpdacc(nprandor,pos,distan,boxdim,rcut, eps, sigma,\
                gam,Kb,tset,dt,stdevdpd,mean,neilist,point,force,w,Up)
            fener.write("#The simulation is run at canonical ensemble\n")
            fener.write("#DPD thermostat is use, target temperature is "+str(tset)+"\n")

        else:
            print("Thers is a set temperature but no thermostat is chosen,"\
                    "please pick an available thermostat (br,la, or dpd)")
            exit()
#        print("first check") 
        #loop over the steps
    for i in range(nstep):
        if mdfun.check(pos,pos_prev,rskin,boxdim) ==1:
            neilist, point = mdfun.neilist(pos, boxdim, rneigh,distan)
            pos_prev=np.copy(pos)
            nlist_update+=1
        Temp, KE = mdfun.kinetic(pos, mass, Kb, Temp, KE)
        ET = Up + KE
        distan = np.zeros((4), order='F')
        pos = mdfun.vv1(pos,dt,force)
        force2=np.copy(force)
        if (tset-0<10**(-6)) or thermo=="br":
            force,w,Up=mdfun.nlistacc(pos,distan,boxdim,rcut,eps,sigma,neilist,point,force,w,Up)
        else:
            if thermo=="la":
                force,w,Up=mdfun.langevin(nprandor,pos,distan,boxdim,rcut,eps,sigma, mass, gam, \
                    neilist,point,stdev,mean,force,w,Up)
            if thermo=="dpd":
                force,w,Up=mdfun.dpdacc(nprandor,pos,distan,boxdim,rcut, eps, sigma,\
                    gam,Kb,tset,dt,stdevdpd,mean,neilist,point,force,w,Up)
        pos =mdfun.vv2(pos,dt,force,force2)
        #print(pos)
        if (i%1000 ==0):
            print("There are "+str(nstep-i)+" steps remaining")
        if ((tset-0)>10**(-6)) and (i%10 ==0) and thermo=="br":
            Temp=0
            pos, Temp=mdfun.berendsen(pos,mass,Kb,tset,dt,tau,Temp)
        #Write data to files
        fxyz.write(str(n) + '\n' + '\n')
        for k in range(n):
            pos[k][0] = pos[k][0] - math.floor(pos[k][0]/boxdim[0])*boxdim[0]
            pos[k][1] = pos[k][1] - math.floor(pos[k][1]/boxdim[1])*boxdim[1]
            pos[k][2] = pos[k][2] - math.floor(pos[k][2]/boxdim[2])*boxdim[2]
            fxyz.write('LJP' + ' ' + str(pos[k][0]) + ' ' +str(pos[k][1]) + ' ' + str(pos[k][2]) + '\n')
        #velocity autocorrelation function
        if vacf != None:
            nframe=i+1
            vel_rec, vacf_rec, vacf_norm = mdfun.autocor(pos, vel_rec, vacf_rec, vacf_norm, nframe, n, nstep)

        #Pressure calculation
        Pcc =  16*3*math.pi*(n/V)*eps*sigma**3 *(2/3*(sigma/rcut)**9 - (sigma/rcut)**3) #pressure correction
#        P = (n*Kb*Temp - w*2/3)/(V) + Pcc #pressure calculated with correction
        P = (2*(KE - w)/(3*V))*(100)/6.02 #calculate the pressure and convert it to bar
#         ET = Up + KE
        fener.write("{:<25.7f}{:<25.7f}{:<25.7f}{:<13.2f}{:<25.7f}{:<10d}\n".format(Up, KE, ET, Temp,P,nlist_update))
    #close the file handle
    fxyz.close()
    fener.close()
    if vacf != None:
        fvacf=open(vacf,'w')
        for i in range(nstep):
            fvacf.write("{:13.5f}{:13.5f}\n".format(i,  vacf_rec[i]/vacf_norm[i]))
        fvacf.close()
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
                        default=0.001,type=float, metavar='')
    parser.add_argument("-tset", "--tset", help="Temperature for the MD simulation"\
                            "0 would mean the simulation is under microcanonical ensemble",
                        default=0,type=float, metavar='')
    parser.add_argument("-nstep", "--nstep", help="Number of steps for the simulation",
                        default=1000,type=int, metavar='')
    parser.add_argument("-thermo", "--thermo", help="choice of thermostat, be for Berendsen,"\
                        "la for Langevin, and dpd for dissipative particle dynamics."\
                        "The default is Berendsen",
                        default="br", metavar='')
    parser.add_argument("-vacf", "--vacffile", help="specify the name of the energy (output) file",
                        default=None, metavar='')



    args = parser.parse_args()
    return args

#-------------------------------------------------------------------------
#Now lets get into the MD functions
#-------------------------------------------------------------------------
#read the fro file 
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

#Random number generator used for Fortran function
def nprandor(mean,stdev):
    x = np.random.normal(mean,stdev)
    return x



if __name__ == '__main__':
    main()
