#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;




int n_list(int n, double par_fun[][6], int list[], int point[],double x, double y, double z, double rneigh)
{
    double dx;
    double dy;
    double dz;
    double distan[4];
    int nnlist;
    nnlist=0;
//    cout<<x<<" xdis "<<y<<" ydis "<<z<<" zdis "<<endl;

    for(int i=0; i<(n-1);i++)
    {
        point[i]=nnlist;
        for(int j=(i+1); j<n; j++)
        {
//            cout<<"j= "<<j<<endl;
            dx = par_fun[i][0]-par_fun[j][0];
//            cout<<"lround x "<<lround(dx/x)<<endl;
            distan[0]=dx-x*lround(dx/x);
            dy = par_fun[i][1]-par_fun[j][1];
            distan[1]=dy-y*lround(dy/y);
            dz = par_fun[i][2]-par_fun[j][2];
            distan[2]=dz-z*lround(dz/z);
            distan[3]=distan[0]*distan[0]+distan[1]*distan[1]+distan[2]*distan[2];
//            cout<<"dx= "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;
//            cout<<"rdx= "<<distan[0]<<" rdy= "<<distan[1]<<" rdz= "<<distan[2]<<endl;
//            cout<<"Distance between "<<i<<" and "<<j<<" is "<< distan[3]<<endl;
            if ( distan[3] < rneigh*rneigh )
            {
                list[nnlist] = j;
                nnlist = nnlist+1;
            }
        }
//        cout<<point[i]<<endl;
    }
//    cout<<nnlist<<endl;
    point[n-1]=nnlist;
//    for(int i=0; i<nnlist;i++)
//    {
//        cout<<"list for "<<i<<" "<<list[i]<<endl;
//    }
//    for(int i=0; i<n;i++)
//    {
//        cout<<"point "<<i<<" "<<point[i]<<endl;
//    }



    return 0;
}





