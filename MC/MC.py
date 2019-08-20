#Monte Carlo for NVT system
#MSM 2018 Dr. Sarupria
#Code written by Tim Yuan
import numpy, math, pylab, datetime, random
#-------------------------------------------------------------------
#Setting parameters
position = input("Random generation of lattice initial position (random/lattice): ")
pos = []
PE = []
sigma = 1.0                     #sigma value dimensionless
eps = 1.0                       #eps dimensionless value
n_step = int(input("Specify the Number of steps you would like to run: "))                   #number of steps
delta = sigma * 0.05                    #delta value represents the step change
Kb = 1.0                        #dimensionless Kb
T = 1.0                         #dimensionless Temperature
beta = 1.0/(Kb*T)               #beta value
rcut = 3*sigma                  #cutoff distance
rcut2 = rcut**2                 
density = float(input("Specify the number density you would like to simulate: "))                          #density that we are interested
#-------------------------------------------------------------------
if position == "random":
    numpar =float(input("number of particles if you want to generate from random position(eg: 500): "))                          #number of particles
    V = numpar/density
    boxlength = numpy.cbrt(V)
elif position =="lattice":
    grid = float(input("cubic root number of particles if you want to generate from lattice config(input 5 will give 125 particles): "))                      #1D grid number
    gridL = numpy.cbrt(1/float(density))                     #Grid length
    V = (grid*gridL)**3
    numpar = grid**3                       #1 atom is assigned in every grid
    boxlength =  gridL * grid                #total length 
#-------------------------------------------------------------------
##generating initial positions with lattice config 
def inipos(grid, gridL):
    for x in range(int(grid)):
        for y in range(int(grid)):
            for z in range(int(grid)):
                xx = x * (gridL) + gridL/2
                yy = y * (gridL) + gridL/2
                zz = z * (gridL) + gridL/2
                boxgrid = [xx, yy, zz]
                pos.append(boxgrid)
    return pos
#-------------------------------------------------------------------
#for random initial position
##generating initial positions with random choice
def randompos(numpar, boxlength):
    for i in range(int(numpar)):
        randompos = [random.uniform(0.0, boxlength), random.uniform(0.0, boxlength), random.uniform(0.0, boxlength)]
        pos.append(randompos)
    return pos
##-------------------------------------------------------------------
#Periodic boundary condition
def dist(x, y, L):
    d_x = (x[0] - y[0]) 
    d_x = d_x - L*round(d_x/L)
    d_y = (x[1] - y[1])  
    d_y =  d_y - L*round(d_y/L)
    d_z = (x[2] - y[2])
    d_z =  d_z - L*round(d_z/L)
    return d_x, d_y, d_z, ((d_x)**2 + (d_y)**2 + (d_z)**2)
#-------------------------------------------------------------------
#Calculating potential from LJ
def potential(pos, numpar):
    Up = 0
    for i in range(int(numpar)-1):
        for j in range(i+1, int(numpar)):
            x = pos[i]
            y = pos[j]
            d_x, d_y, d_z, dis2 = dist(x , y , boxlength )
            if dis2 == 0:
                Up_i = 0
            elif dis2 < rcut2:
                Up_i = 4.0 * eps *  ((sigma**12)*(dis2**(-6))- (sigma**6)*(dis2**(-3))) - ( 4.0 * eps *  ((sigma**12)*(rcut2**(-6))- (sigma**6)*(rcut2**(-3))))
            else:
                Up_i = 0
            Up = Up + Up_i
    return Up
#
#-------------------------------------------------------------------
#Implementing MC
if position =="random":
    pos = randompos(numpar, boxlength)
elif position =="lattice":
    pos = inipos(grid,gridL)
f = open('MCrcut' + str(datetime.datetime.now().strftime("%Y-%m-%d%H-%M-%S")), 'a')
ff = open('MCposition' + str(datetime.datetime.now().strftime("%Y-%m-%d%H-%M-%S")) + ".xyz", 'a')
Up = potential(pos, numpar)
b = [0.0, 0.0, 0.0]
naccept = 0
aveUp = 0
upt = 0
for nsp in range(n_step):
    Upt = 0
    Uptn = 0
    a = random.choice(pos)
    for d in pos:
        if d != a:
            distx, disty, distz, dist2 = dist( a, d, boxlength)
            if dist2 ==0:
                Up_i = 0
            elif dist2 < rcut2:
                Up_i = 4.0 * eps *  ((sigma**12)*(dist2**(-6))- (sigma**6)*(dist2**(-3)))
            else: 
                Up_i = 0
            Upt = Upt + Up_i
    b[0] = (a[0] + random.uniform(-delta, delta)) 
    b[0] = b[0] - math.floor(b[0]/boxlength)*boxlength
    b[1] = (a[1] + random.uniform(-delta, delta))
    b[1] = b[1] - math.floor(b[1]/boxlength)*boxlength
    b[2] = (a[2] + random.uniform(-delta, delta))
    b[2] = b[2] - math.floor(b[2]/boxlength)*boxlength
#    print(a, b)
    for c in pos:
        if c != a:
            newdistx, newdisty, newdistz, newdist2 = dist( b, c, boxlength)
            if newdist2 < rcut2:
                Up_n_i = 4.0 * eps *  ((sigma**12)*(newdist2**(-6))- (sigma**6)*(newdist2**(-3)))
            else:
                Up_n_i = 0
            Uptn = Uptn + Up_n_i
        deltaUp = Uptn - Upt
    if deltaUp < 0.0:
        a[:] = b
        Up = Up + deltaUp 
        naccept = naccept + 1
    elif random.uniform(0.0, 1.0) < numpy.exp(-beta*deltaUp):
        a[:] = b
        Up = Up + deltaUp 
        naccept = naccept + 1
    else: 
        a[:] = a
        Up = Up
        naccept = naccept + 0 
    upt = upt + Up
    aveUp = upt/(nsp+1.0)
    acceptratio = naccept/(nsp+1.0)
    print(Up, aveUp, deltaUp, naccept, acceptratio)
    f.write(str(Up) + ' ' + str(aveUp) + ' ' + str(naccept)+ ' ' + str(acceptratio) + '\n')
    ff.write(str(numpar) + '\n' + '\n')
    for j in range(int(numpar)):
        ff.write('LJP' + ' ' + str(pos[j][0]) + ' ' +str(pos[j][1]) + ' ' + str(pos[j][2]) + '\n')
