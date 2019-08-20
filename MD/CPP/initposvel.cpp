#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <stdio.h>
using namespace std;


int initpos(int n, double  pos_vel[][6], double rho)
{
    double lattice = cbrt(rho);
    int a = round(cbrt(n));

    int num=0;
    for (int i =0; i<a; i++)
    {
        for (int j =0; j<a; j++)
        {
            for (int k =0; k<a; k++)
            {
                pos_vel[num][0]=i*lattice+lattice/2.;
                pos_vel[num][1]=j*lattice+lattice/2.;
                pos_vel[num][2]=k*lattice+lattice/2.;
                num+=1;
            }
        }
    }        
    return 0;
}


int initvel(int n, double pos_vel[][6], double v_std)
{
    int a = 0;
    random_device rd;
    mt19937 mt(rd());
    normal_distribution<double> dist(a, v_std);
    for (int i=0; i<n; i++)
    {
        pos_vel[i][3]=dist(mt);
        pos_vel[i][4]=dist(mt);
        pos_vel[i][5]=dist(mt);
    }
    double sumvx=0;
    double sumvy=0;
    double sumvz=0;
    for (int i=0; i<n; i++)
    {
        sumvx += pos_vel[i][3];
        sumvy += pos_vel[i][4];
        sumvz += pos_vel[i][5];
    }
    double avevx=sumvx/n;
    double avevy=sumvy/n;
    double avevz=sumvz/n;
//    double fsumvx=0.;
//    double fsumvy=0.;
//    double fsumvz=0.;
    for (int i=0; i<n; i++)
    {
        pos_vel[i][3]=pos_vel[i][3]-avevx;
        pos_vel[i][4]=pos_vel[i][4]-avevy;
        pos_vel[i][5]=pos_vel[i][5]-avevz;
//        fsumvx+= pos_vel[i][3];
//        fsumvy+= pos_vel[i][4];
//        fsumvz+= pos_vel[i][5];

    }
//    cout<<fsumvx<<"\t"<<fsumvy<<"\t"<<fsumvz<<endl;
    return 0;
}

int output(int n, double pos_vel[][6],double x,double y ,double z)
{
//    ofstream outfile ("initial.gro");
//    outfile << "initial"<< "\n";
//    outfile << n << "\n";  
    FILE * outfile;
    outfile = fopen ("initial.gro", "w");
    fprintf (outfile, "%s\n", "initial.gro");
    fprintf (outfile, "%d\n", n);
    for (int i=0;i<n;i++)
    {
//        outfile <<fixed<<setprecision(8)<<"\t"<<i<<"LJP\t"<<"LJP\t"<<i+1<<"\t"<<pos_vel[i][0]<<"\t"<< pos_vel[i][1]<<"\t"<<pos_vel[i][2]<<"\t"<<pos_vel[i][3]<<"\t"<<pos_vel[i][4]<<"\t"<<pos_vel[i][5]<<"\n";
//        printf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f", i, "LJP", "LJP",i,pos_vel[i][0],pos_vel[i][1],pos_vel[i][2],pos_vel[i][3],pos_vel[i][4],pos_vel[i][5]);
       fprintf  (outfile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", (i+1), "LJP", "LJP",(i+1),pos_vel[i][0],pos_vel[i][1],pos_vel[i][2],pos_vel[i][3],pos_vel[i][4],pos_vel[i][5]);
    }    
//    outfile << x <<"\t"<< y <<"\t"<< z <<"\t" << "\n";
//    outfile.close();
    fprintf (outfile, "%f\t %f\t %f\n", x,y,z);
    fclose (outfile);
    return 0;
} 
