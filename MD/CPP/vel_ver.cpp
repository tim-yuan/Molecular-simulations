#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;




int vel_ver1(double par_fun[][6],double acc[][3],double dt, int n)
{
    for (int i=0; i<n; i++)
    {
        par_fun[i][0] = par_fun[i][0] +dt*par_fun[i][3]+ 0.5*dt*dt*acc[i][0];
        par_fun[i][1] = par_fun[i][1] +dt*par_fun[i][4]+ 0.5*dt*dt*acc[i][1];
        par_fun[i][2] = par_fun[i][2] +dt*par_fun[i][5]+ 0.5*dt*dt*acc[i][2];
    }


    return 0;
}


int vel_ver2(double par_fun[][6],double acc[][3],double acc2[][3], double dt, int n)
{
    for (int i=0; i<n; i++)
    {
        par_fun[i][3] = par_fun[i][3] + 0.5*dt*(acc[i][0]+acc2[i][0]);
        par_fun[i][4] = par_fun[i][4] + 0.5*dt*(acc[i][1]+acc2[i][1]);
        par_fun[i][5] = par_fun[i][5] + 0.5*dt*(acc[i][2]+acc2[i][2]);
    }


    return 0;
}



