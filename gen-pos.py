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
np.set_printoptions(threshold=sys.maxsize)


#args=parseargs()


def main():
    args = parseargs()
    f_out = args.outputfile
    n=args.num
    x=float(args.xdim)
    y=float(args.ydim)
    z=float(args.zdim)
    vmax=args.vmax
    pos = np.zeros((int(n), 6),order='F')
    sumvx = 0
    sumvy = 0
    sumvz = 0
    sumv2 = 0
    for i in range (n):
        pos[i][0] = np.random.uniform(0,x)
        pos[i][1] = np.random.uniform(0,y)
        pos[i][2] = np.random.uniform(0,z)
        pos[i][3] = np.random.normal(0,vmax)
        pos[i][4] = np.random.normal(0,vmax)
        pos[i][5] = np.random.normal(0,vmax)
        #remove center of mass motions
        sumvx +=  pos[i][3]   
        sumvy +=  pos[i][4]   
        sumvz +=  pos[i][5]
    avevx = sumvx/n   
    avevy = sumvy/n   
    avevz = sumvz/n
    #To check we indeed removed the center of mass
    sumx = 0
    sumy = 0
    sumz = 0

    for i in range (n):
        pos[i][3] = (pos[i][3] - avevx)   
        pos[i][4] = (pos[i][4] - avevy)   
        pos[i][5] = (pos[i][5] - avevz)   
        sumx += pos[i][3]
        sumy += pos[i][4]
        sumz += pos[i][5]

    if (sumx<10**-3 and sumy<10**-3 and sumz<10**-3):
        f_o_hand = open(f_out, 'w')
        f_o_hand.write("Initial structure from random placement\n")
        f_o_hand.write(str(n)+"\n")
        for i in range (n):
            f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".format(str(i+1), "LJs", "LJ",str(i+1),(pos[i][0]), (pos[i][1]), (pos[i][2]), (pos[i][3]), (pos[i][4]), (pos[i][5]) ))
        boxlength  = str("        "+str(x)+"       "+str(y)+"       "+str(z)+"     \n")
        f_o_hand.write(boxlength)
        f_o_hand.close()
    else:
        print("Unable to create the inital configuration. Please check the code")



def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Tim Yuan\nSarupria Research Group\n"+modified+"\n\n"
                        "Perform energy minimization"
                        )
    parser.add_argument("-o", "--outputfile", help="specify the name of the output file",
                        default="initial.gro", metavar='')
    parser.add_argument("-n", "--num", help="Specificy the number of particles",
                        default=500, type=int, metavar='')
    parser.add_argument("-x", "--xdim", help="Specificy the x dimension of the box",
                        default=8, type=float, metavar='')
    parser.add_argument("-y", "--ydim", help="Specificy the y dimension of the box",
                        default=8, type=float, metavar='')
    parser.add_argument("-z", "--zdim", help="Specificy the z dimension of the box",
                        default=8, type=float, metavar='')
    parser.add_argument("-vmax", "--vmax", help="Specificy maximum velocity for x, y, and z component",
                        default=1, type=float, metavar='')
  
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    main()

