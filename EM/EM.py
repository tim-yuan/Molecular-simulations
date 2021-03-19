# Tim Yuan
# Sarupria Group
# Clemson University
# Project for Numerical Method
# This Code takes an input gro file and perform energy minimization for a simple LJ system
# Since this is an Energy Minimization code, no velocities will be considered


code ="EM.py"
modified ="20th April, 2019"

import sys
import os.path
import numpy as np
import argparse
import pylab
import math
import scipy.special
import EM
import datetime
np.set_printoptions(threshold=sys.maxsize)


#args=parseargs()


def main():
#-------------------------------------------------------------------------
#Define some parameters
#-------------------------------------------------------------------------

    args = parseargs()
    f_in = args.inputfile
    f_out = args.outputfile
    eps = args.eps
    sig = args.sig
    ss = args.stepsize
    beta = args.beta
    nmax = args.nmax
    f_ener = args.energyfile
    gamma = args.gamma
    tol = float(args.tol)
    alg = args.algorithm
#    pos_vec=[]
    rcut = 3*sig
    rcut2 = rcut**2


#-------------------------------------------------------------------------
#Read the input file and store necessary informations
#-------------------------------------------------------------------------

    pos=[]
    if os.path.isfile(f_in):
        f_in_hand = open(f_in, 'r')
        grofile=f_in_hand.readlines()
        n = int(grofile[1].strip())
        for line in grofile[2:(len(grofile)-1)]:
            ele = [line[0:5].strip(), line[5:10].strip(),line[10:15].strip(),\
                line[15:20].strip(),line[20:28].strip(),line[28:36].strip(),\
                line[36:44].strip(),line[44:52].strip(),line[52:60].strip(),
                line[60:68].strip()]
            pos_vec=[float(ele[4]),float(ele[5]),float(ele[6])]
            pos.append(pos_vec)
        x=float(grofile[-1].split()[0])
        y=float(grofile[-1].split()[1])
        z=float(grofile[-1].split()[2])
#        print(x,y,z)
    else:
        print("No file found, please input the correct name of the gro file")
        exit()
    pos=np.asarray(pos)
#    print(len(pos))
#    print(n)
#    print(type(n))
#    print(pos)


    #-------------------------------------------------------------------------
    #Call the energy minimization function
    #-------------------------------------------------------------------------
    if  alg == 'sd':
        pos, Up, step=em_sd(pos,x,y,z,rcut,eps,sig,nmax,beta,ss,tol,gamma)
    elif alg =='sd-1':
        pos, Up, step=em_sd1(pos,x,y,z,rcut,eps,sig,nmax,beta,ss,tol,gamma)
    elif alg =='sdrd':
        pos, Up, step=em_sd_rand(pos,x,y,z,rcut,eps,sig,nmax,beta,ss,tol,gamma)
 
    #-------------------------------------------------------------------------
    #We write the energy minimized configuration to the file
    #-------------------------------------------------------------------------
    f_o_hand = open(f_out, 'w')
    f_o_hand.write(grofile[0])
    f_o_hand.write(grofile[1])
    for i in range(n):
        #Apply pbc boundary conditions and bring the atoms back to the box
        pos[i][0]=pos[i][0]-np.floor(pos[i][0]/x)*x
        pos[i][1]=pos[i][1]-np.floor(pos[i][1]/y)*y
        pos[i][2]=pos[i][2]-np.floor(pos[i][2]/z)*z
        #Now we write our output file
        line=grofile[i+2]
        ele = [line[0:5].strip(), line[5:10].strip(),line[10:15].strip(),\
            line[15:20].strip(),line[20:28].strip(),line[28:36].strip(),\
            line[36:44].strip(),line[44:52].strip(),line[52:60].strip(),
            line[60:68].strip()]

        f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(str(ele[0]), str(ele[1]), str(ele[2]), str(ele[3]), pos[i][0],\
            pos[i][1], pos[i][2], float(ele[7]), float(ele[8]), float(ele[9]) ))
    f_o_hand.write(grofile[-1])
    f_o_hand.close()

    #-------------------------------------------------------------------------
    #We write the energy to the file
    #-------------------------------------------------------------------------
    f_en_hand = open(f_ener, 'w')
    f_en_hand.write("# Tim Yuan\n# Sarupria Group\n# Clemson University\n# "+str(datetime.datetime.now())+'\n\n\n' )
    f_en_hand.write("The energy minimization used "+str(alg)+" algorithm\n")
    f_en_hand.write("A total of "+str(step)+" steps were done\n")
    f_en_hand.write("The minimal potential is "+str(Up)+"\n")

#-------------------------------------------------------------------------
#Now lets get into the energy minimization stuff
#The em_sd is the energy minimization using steepest descent
#-------------------------------------------------------------------------
def em_sd(pos,x,y,z,rcut,eps,sig,nmax,beta,ss,tol,gamma):
    #-------------------------------------------------------------------------
    #Initial the arrays necessary for calculation
    #-------------------------------------------------------------------------
    force = np.zeros((len(pos), 3),order='F')
    Up=0
    force_new = np.zeros((len(pos), 3),order='F')
    Up_new=0
    Up,force = EM.potential(pos,x,y,z,rcut,eps,sig,Up,force)
    Up_prev=np.copy(Up)
    force_prev=np.copy(force)
    pos_i=np.copy(pos)
    pos_prev=np.copy(pos)
    for i in range(nmax):
        m=1
        Up,force = EM.potential(pos,x,y,z,rcut,eps,sig,Up,force)
        print("Step "+str(i)+", the potential energy is "+str(Up))
        force1d=np.reshape(force,(3*(len(pos)),1))
        armijo_creteria=-(np.matmul(np.ndarray.transpose(force1d),force1d)/np.linalg.norm(force))
        for j in range(100):
            pos_i=pos+((beta)**m)*ss*force
            Up_new, force_new=EM.potential(pos_i,x,y,z,rcut,eps,sig,Up_new,force_new)
            #This if statement 
            #if (Up_new-Up)<=0:
            #This if statement implements Armijo rule
#            print((gamma*((beta)**m)*ss*armijo_creteria))
            if (Up_new-Up)<=(gamma*((beta)**m)*ss*armijo_creteria):
                pos=np.copy(pos_i)
                Up=Up_new
                break
            else:
                m+=1
        if (np.linalg.norm(force)<tol):
            print("\n\nThe position is not changed")
            print("Its likely that you have found a minimum or saddle point")
            print("The energy minimization stoped after "+str(i)+" number of interations")
            print("The final potential energy obtained is "+str(Up)+"in units of kT")
            step=i
            break
        else:
            step=nmax
        pos_prev=np.copy(pos)
    if i==(nmax-1):
        print("\n\nMaximum number of interation has been reached")
        print("You can change the number of interactions using -nmax")
        print("The final potential energy obtained is "+str(Up)+"in units of kT")
    return pos, Up, step

#-------------------------------------------------------------------------
#The em_sd1 moves 1 particle at each time
#The algorithm is still steepest descent
#-------------------------------------------------------------------------
def em_sd1(pos,x,y,z,rcut,eps,sig,nmax,beta,ss,tol,gamma):
    #-------------------------------------------------------------------------
    #Initial the arrays necessary for calculation
    #-------------------------------------------------------------------------
    force = np.zeros((len(pos), 3),order='F')
    Up=0
    force_new = np.zeros((len(pos), 3),order='F')
    Up_new=0

    Up,force = EM.potential(pos,x,y,z,rcut,eps,sig,Up,force)
    Up_prev=np.copy(Up)
    force_prev=np.copy(force)
    pos_i=np.copy(pos)
    pos_prev=np.copy(pos)

    for i in range(nmax):
        m=1

        Up,force = EM.potential(pos,x,y,z,rcut,eps,sig,Up,force)
#        Up_prev=np.copy(Up)
        force_abssum=np.sum(np.abs(force), axis=1)
        V_maxforce = np.amax(force_abssum)
        index_maxforce=np.where(force_abssum ==  V_maxforce)
        print("The max force is on "+str(index_maxforce)+" atom")
        print("Step "+str(i)+", the potential energy is "+str(Up))
        force1d=np.reshape(force,(3*(len(pos)),1))
        armijo_creteria=-np.matmul(np.ndarray.transpose(force1d),force1d)/np.linalg.norm(force)
        for j in range(100):
            for k in range(len(index_maxforce)):
                pos_i[index_maxforce[k]]=pos[index_maxforce[k]]+((beta)**m)*ss*force[index_maxforce[k]]
            Up_new, force_new=EM.potential(pos_i,x,y,z,rcut,eps,sig,Up_new,force_new)
#            if (Up_new-Up)<=0:
            if (Up_new-Up)<=(gamma*((beta)**m)*ss*armijo_creteria):
                pos=np.copy(pos_i)
                Up=Up_new
                break
            else:
                m+=1
        if (np.linalg.norm(force)<tol):
            print("\n\nThe potential is not changed")
            print("Its likely that you have found a minimum or saddle point")
            print("The energy minimization stoped after "+str(i)+" number of interations")
            print("The final potential energy obtained is "+str(Up)+"in units of kT")
            step = i
            break    
        else:
            step=nmax
        pos_prev=np.copy(pos)
    if i==(nmax-1):
        print("\n\nMaximum number of interation has been reached")
        print("You can change the number of interactions using -nmax")
        print("The final potential energy obtained is "+str(Up)+"in units of kT")
    return pos, Up, step
        
#-------------------------------------------------------------------------
#The em_sd1 moves 1 particle at each time
#The algorithm is still steepest descent
#-------------------------------------------------------------------------
def em_sd_rand(pos,x,y,z,rcut,eps,sig,nmax,beta,ss,tol,gamma):
    #-------------------------------------------------------------------------
    #Initial the arrays necessary for calculation
    #-------------------------------------------------------------------------
    force = np.zeros((len(pos), 3),order='F')
    Up=0
    force_new = np.zeros((len(pos), 3),order='F')
    Up_new=0

    Up,force = EM.potential(pos,x,y,z,rcut,eps,sig,Up,force)
    Up_prev=np.copy(Up)
    force_prev=np.copy(force)
    pos_i=np.copy(pos)
    pos_prev=np.copy(pos)

    for i in range(nmax):
        m=1
        Up,force = EM.potential(pos,x,y,z,rcut,eps,sig,Up,force)
        n_rd=np.random.randint(1,len(pos))
        force1d=np.reshape(force,(3*(len(pos)),1))
        armijo_creteria=-np.matmul(np.ndarray.transpose(force1d),force1d)/np.linalg.norm(force)
        n_list=[]

        for k in range(n_rd):
            rand_select=np.random.randint(len(pos))
            if rand_select not in n_list:
                n_list.append(np.random.randint(len(pos)))

        print("Step "+str(i)+", the potential energy is "+str(Up))
        print("You have randomly pick "+str(len(n_list))+" number of particles")
        for j in range(100):
            for k in range(len(n_list)):
                pos_i[n_list[k]]=pos[n_list[k]]+((beta)**m)*ss*force[n_list[k]]
            Up_new, force_new=EM.potential(pos_i,x,y,z,rcut,eps,sig,Up_new,force_new)
#            if (Up_new-Up)<=0:
            if (Up_new-Up)<=(gamma*((beta)**m)*ss*armijo_creteria):
                pos=np.copy(pos_i)
                Up=Up_new
                break
            else:
                m+=1
        if (np.linalg.norm(force)<tol):
            print("\n\nThe position is not changed")
            print("Its likely that you have found a minimum or saddle point")
            print("The energy minimization stoped after "+str(i)+" number of interations")
            print("The final potential energy obtained is "+str(Up)+"in units of kT")
            step = i
            break    
        else:
            step=nmax
        pos_prev=np.copy(pos)
    if i==(nmax-1):
        print("\n\nMaximum number of interation has been reached")
        print("You can change the number of interactions using -nmax")
        print("The final potential energy obtained is "+str(Up)+"in units of kT")
    return pos, Up, step

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Tim Yuan\nSarupria Research Group\n"+modified+"\n\n"
                        "Perform energy minimization"
                        )
    parser.add_argument("-i", "--inputfile", help="input gro file",
                        default="input.gro", metavar='')
    parser.add_argument("-o", "--outputfile", help="specify the name of the output file",
                        default="config.gro", metavar='')
    parser.add_argument("-energy", "--energyfile", help="specify the name of the energy (output) file",
                        default="energy.xvg", metavar='')
    parser.add_argument("-eps", "--eps", help="epsilon value for the LJ particle",
                        default=1,type=float, metavar='')
    parser.add_argument("-sig", "--sig", help="sigma value for the LJ particle",
                        default=1,type=float, metavar='')
    parser.add_argument("-ss", "--stepsize", help="Stepsize for Armijo rule",
                        default=0.05,type=float, metavar='')
    parser.add_argument("-beta", "--beta", help="beta value for Armijo rule (0-1), default 0.5",
                        default=0.5,type=float, metavar='')
    parser.add_argument("-gamma", "--gamma", help="gamma value for Armijo rule (0-1), default 0.5",
                        default=0.5,type=float, metavar='')
    parser.add_argument("-nmax", "--nmax", help="Maximimum number of steps for the calculation, default 5,000",
                        default=2000,type=int, metavar='')
    parser.add_argument("-tol", "--tol", help="Tolerance, default 500",
                        default=500,type=float, metavar='')
    parser.add_argument("-alg", "--algorithm", help="pick 'sd' for steepest descent, or 'sd-1'"\
                            " for 1 atom steepest descent, or 'sdrd' for randomly picking"\
                            " x number of particles and perform energy minimization. default sd",
                        default='sd', metavar='')
 
    args = parser.parse_args()
    return args




if __name__ == '__main__':
    main()
