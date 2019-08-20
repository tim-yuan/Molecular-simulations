#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;




int pbc(int n, double par_fun[][6],double xx, double yy, double zz )
{
    for (int i=0; i<n; i++)
    {
        par_fun[i][0]=par_fun[i][0]-xx*floor(par_fun[i][0]/xx);
        par_fun[i][1]=par_fun[i][1]-yy*floor(par_fun[i][1]/yy);
        par_fun[i][2]=par_fun[i][2]-zz*floor(par_fun[i][2]/zz);

    }
    return 0;
}





