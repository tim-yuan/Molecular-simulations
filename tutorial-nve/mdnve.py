#This code is written by Tim Yuan
#Created on the 14th July, 2021
#For the python fortran tutorial

import argparse
import numpy as np
import mdnve 
import sys
import os.path

code="mdnve.py"
modified="14th July, 2022"
def main():
    args = parseargs()
    f_in = args.inputfile
    eps = args.eps
    sigma = args.sig
    dt = args.timestep
    nstep = args.nstep
    #--------------------------------------------------------------
    #Input data
    #--------------------------------------------------------------
    rcut = 3*sigma
    mass = 1
    kb = 1
    #--------------------------------------------------------------
    #Read input file
    #--------------------------------------------------------------
    readfile(f_in)




    #--------------------------------------------------------------
    #Call fortran function to run MD
    #--------------------------------------------------------------

    ###############################################################
    ###############################################################
    ###############################################################
    #           Exercise for you:                                 #
    #           Call the run function                             #
    #           from fortran, and perform                         #
    #           MD calculation using                              #
    #           config.gro as a starting poin                     #
    ###############################################################
    ###############################################################
    ###############################################################







#---------------------------------------------------------------------
#read the fro file 
#---------------------------------------------------------------------
def readfile(f_in):
    pm.pos=[]
    if os.path.isfile(f_in):
        f_in_hand = open(f_in, 'r')
        grofile=f_in_hand.readlines()
        pm.n = int(grofile[1].strip())
        for line in grofile[2:(len(grofile)-1)]:
            ele = [line[0:5].strip(), line[5:10].strip(),line[10:15].strip(),\
                line[15:20].strip(),line[20:28].strip(),line[28:36].strip(),\
                line[36:44].strip(),line[44:52].strip(),line[52:60].strip(),\
                line[60:68].strip()]
            pm.pos_vec=[float(ele[4]),float(ele[5]),float(ele[6]),float(ele[7]),float(ele[8]),float(ele[9])]
            pm.pos.append(pm.pos_vec)
        x=float(grofile[-1].split()[0])
        y=float(grofile[-1].split()[1])
        z=float(grofile[-1].split()[2])
        pm.box=[x,y,z]
    else:
        print("No file found, please input the correct name of the gro file")
        exit()
    pm.box=np.asarray(pm.box, dtype=np.float64)
    pm.pos=np.asarray(pm.pos, dtype=np.float64)



###########################################################################
# Class for the parameters                                                #
###########################################################################
class pm:
    '''This class includes the parameters from input file'''

#---------------------------------------------------------------------
#Input argument
#---------------------------------------------------------------------
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Tim Yuan\nSarupria Research Group\n"+modified+"\n\n"
                        "Perform energy minimization"
                        )
    parser.add_argument("-i", "--inputfile", help="input gro file",
                        default="input.gro", metavar='')
    parser.add_argument("-eps", "--eps", help="epsilon value for the LJ particle",
                        default=1,type=float, metavar='')
    parser.add_argument("-sig", "--sig", help="sigma value for the LJ particle",
                        default=1,type=float, metavar='')
    parser.add_argument("-dt", "--timestep", help="Time step for the MD simulation",
                        default=0.001,type=float, metavar='')
    parser.add_argument("-nstep", "--nstep", help="Number of steps for the simulation",
                        default=1000,type=int, metavar='')



    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()

