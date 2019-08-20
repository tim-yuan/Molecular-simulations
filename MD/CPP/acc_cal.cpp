#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <iomanip>

using namespace std;




int acc_cal(int n, double par_fun[][6],int list[],int point[],double rcut, double eps, double sigma,double acc[][3], double xx,double yy, double zz,  double Up_pres[2])
{
    double rcut2=rcut*rcut;
    double upcut=4.*eps*((pow(sigma,12.)*pow(rcut2,-6.))-(pow(sigma,6.)*pow(rcut2,-3.)));
    double wx =0.;
    double wy =0.;
    double wz =0.;
    double distan[4];
    double dx;
    double dy;
    double dz;
    double fx_i;
    double fy_i;
    double fz_i;
    double up_i;
    double wx_i;
    double wy_i;
    double wz_i;
    int j;

    Up_pres[0] = 0.;
    Up_pres[1] = 0.;
    wx = 0.;
    wy = 0.;
    wz = 0.;
    fx_i=0.; 
    fy_i=0.;
    fz_i=0.;
    up_i=0.;
    wx_i=0.;
    wy_i=0.;
    wz_i=0.;
    for (int i=0; i<n; i++)
    {
        acc[i][0]=0.0;
        acc[i][1]=0.0;
        acc[i][2]=0.0;
    }
        
    for (int i=0; i<(n-1); i++)
    {
        int jbeg = point[i];
        int jend = point[i+1];
        if (jbeg <= jend)
        {
            for (int jab=jbeg ;jab<jend; jab++  )
            {
                j = list[jab];
                dx = par_fun[i][0]-par_fun[j][0];
                distan[0]=dx-xx*lround(dx/xx);
                dy = par_fun[i][1]-par_fun[j][1];
                distan[1]=dy-yy*lround(dy/yy);
                dz = par_fun[i][2]-par_fun[j][2];
                distan[2]=dz-zz*lround(dz/zz);
//                cout<<"dx "<<distan[0]<<" dy "<<distan[1]<<" dz "<<distan[2]<<endl;
                distan[3]=distan[0]*distan[0]+distan[1]*distan[1]+distan[2]*distan[2];
//                cout<<distan[3];
                if (distan[3]==0)
                {
                    cout<<"Error: calculating mirror imagine"<<endl;
                }
                if (distan[3] < rcut2)
                {
                    fx_i = 4.*eps*(12*(pow(sigma,12)*pow(distan[3],-7))-6*(pow(sigma,6)*pow(distan[3],-4)))*distan[0];    
                    fy_i = 4.*eps*(12*(pow(sigma,12)*pow(distan[3],-7))-6*(pow(sigma,6)*pow(distan[3],-4)))*distan[1];    
                    fz_i = 4.*eps*(12*(pow(sigma,12)*pow(distan[3],-7))-6*(pow(sigma,6)*pow(distan[3],-4)))*distan[2];    
                    acc[i][0] +=fx_i;
                    acc[i][1] +=fy_i;
                    acc[i][2] +=fz_i;
//                    acc[j][0] +=(-fx_i);
//                    acc[j][1] +=(-fy_i);
//                    acc[j][2] +=(-fz_i);
                    acc[j][0] =acc[j][0]-fx_i;
                    acc[j][1] =acc[j][1]-fy_i;
                    acc[j][2] =acc[j][2]-fz_i;
                    up_i = 4*eps*((pow(sigma,12)*pow(distan[3],-6))-(pow(sigma,6)*pow(distan[3],-3)))-upcut;
                    wx_i = fx_i*distan[0];
                    wy_i = fy_i*distan[1];
                    wz_i = fz_i*distan[2];
                    Up_pres[0] += up_i;
                    wx += wx_i;
                    wy += wy_i;
                    wz += wz_i;

//                    cout<<setprecision(16)<<"potential_i"<<up_i<<endl;
                }
            }
        }
//        cout<<setprecision(16)<<acc[i][0]<<" accy "<<acc[i][1]<<" accz "<<acc[i][2]<<endl;
    }
//    cout<<setprecision(16)<<acc[n-1][0]<<" accy last one "<<acc[n-1][1]<<" accz "<<acc[n-1][2]<<endl;
//    for (int i=0; i<n; i++)
//    {
//        cout<<"acc x for "<<i<<" "<< acc[i][0]<< endl;
//        cout<<"acc y for "<<i<<" "<< acc[i][1]<< endl;
//        cout<<"acc z for "<<i<<" "<< acc[i][2]<< endl;
//    }

    Up_pres[1] = (wx+wy+wz)/3.;
//    cout<<setprecision(16)<<"potential energy"<<Up_pres[0]<<endl;
//    cout<<"Up "<<Up_pres[0]<<"\t"<<"w "<<Up_pres[1]<<endl;
    return 0;
}





