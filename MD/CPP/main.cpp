#include <iostream>
#include "mdcpp.h"
#include <cmath>
#include <fstream>
#include <random>
#include <time.h>
#include <iomanip>
using namespace std;




//int main( int argc, char *argv[] )
int main()
{
    clock_t t1,t2;
    string initial= "yes";
    int latsize;
//    cout<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<endl;
//    if (argv[1] !=NULL)
//    {
//        initial == argv[1];
//    }
//    else
//    {
//        initial = "yes";
//    }
//
//    if (argv[2] != NULL)
//    {
//        latsize = stoi(argv[2]);
//    }
//    else
//    {
//        latsize = 5;
//    }

    int n;
    double dim[3];
    double xx;
    double yy;
    double zz;
    string inputfile="initial.gro";
    
//    string inputfile=argv[1];
//------------------------------------------------------------------------------------
    if (initial == "yes")
    {
        latsize = 8;
        n = latsize*latsize*latsize;
    }
    else 
    {
        cout<<"reading from initial.gro"<<endl;
//        ifstream input(inputfile);
//        if(input.is_open())
//        {
//            int linenumber =0;
//            while(!input.eof())
//            {
//                string number_1;
//                getline(input,number_1);
//                linenumber+=1;
//            }
//            int linenum = 0;
//            string length[linenumber];
//
//        ifstream infile(inputfile);
//        if(infile.is_open())
//            while(!infile.eof())
//            {
//                string number;
//                string aaa;
//                getline(infile,number);
//                aaa=number.c_str();
//                length[linenum]=aaa;
//                linenum+=1;
//            }
//            n = atoi(length[1].c_str());
//
//            size_t pos=0;
//            string delimiter="\t";
//            string token;
//            double xyz[3];
//            int coord = 0;
//            linenum=linenum-2;
//            while((pos=length[linenum].find(delimiter)) != string::npos)
//            {
//                token=length[linenum].substr(0,pos);
//                length[linenum].erase(0,pos+delimiter.length());
//                xyz[coord]=stod(token);
//                coord+=1;
//            }
//
//            xx = (xyz[0]);
//            yy = (xyz[1]);
//            zz = (xyz[2]);
//        }

        FILE * input;
        input = fopen (inputfile.c_str(), "r");
        char mystring [100];
        fgets (mystring,100,input);        
        fgets (mystring,100,input);        
        n=stoi(mystring);
        cout<<n<<endl;
        fclose(input);
    }
//-------------------------------------------------------------------------------------
    double pos_vel[n][6];    
    double par_fun[n][6];
    double par_fun_nei[n][6];
    double acc[n][3]={0};
    double acc2[n][3]={0};
    int list[n*(n-1)]={0};
    int point[n];
    double sigma = 1.;
    double eps = 1.;
    double rneigh=3.5*sigma;
    double rskin=sigma;
    double rcut=2.5*sigma;
    double KE;
    double Temp;
    double Up_pres[2]={0,0};
    double dt=0.0005;
    int nsteps=1000 ;
    double mass = 1.0;
    double Kb = 1.0;
//    if (argv[3] !=NULL)
//    {
//        nsteps = stoi(argv[3]);
//    }
//    else
//    {
//        nsteps = 10000;
//    }




//-------------------------------------------------------------------------------------
    if (initial == "yes")
    {
        double rho = 1.;
        double v_std=1.;
        cout<<"Creating from a lattice position with "<<n<<" particles, with density "<<rho<<", v_std "<<v_std<<endl;
        xx = cbrt(n/rho);
        yy = cbrt(n/rho);
        zz = cbrt(n/rho);
        for (int i=0;i<n;i++)
        {
            initpos(n,pos_vel,rho);
            initvel(n,pos_vel,v_std);
            output(n,pos_vel, xx, yy, zz);
        }
    }    
//---------------------------------------------------------------------------------------
    readfile(inputfile,n,par_fun,dim);
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------
    xx = dim[0];
    yy = dim[1];
    zz = dim[2];

    n_list(n,par_fun,list,point,xx,yy,zz,rneigh);


    ofstream xyzfile ("md.xyz");
    ofstream enefile ("energy.xvg");
    xyzfile <<n<<"\n"<<"\n";    
    enefile <<"steps"<<"\t"<<"Potential energy"<<"\t"<<"pressure"<<"\t"<<"Kinetic Energy"<<"\t"<<"Total Energy"<<"\t"<<"Temperature"<<endl;
    for (int i=0; i<n; i++)        
    {
        xyzfile <<fixed<<setprecision(8)<<"LJP"<<"\t"<<par_fun[i][0]<<"\t"<<par_fun[i][1]<<"\t"<<par_fun[i][2]<<"\n";
        par_fun_nei[i][0] = par_fun[i][0];
        par_fun_nei[i][1] = par_fun[i][1];
        par_fun_nei[i][2] = par_fun[i][2];
//        par_fun_nei[i][3] = par_fun[i][3];
//        par_fun_nei[i][4] = par_fun[i][4];
//        par_fun_nei[i][5] = par_fun[i][5];
    }

    t1=clock();
    for (int j=0; j<nsteps; j++)
    {
        if (j%2000==0)
        {
            cout<<"Number of steps to finish is: "<<(nsteps-j)<<endl;
        }

//----------------------------------------------------------------------------------------
//Neighbour list update
//----------------------------------------------------------------------------------------
        double dispmx;
        double dispmy;
        double dispmz;
        double maxdisp;
        dispmx = 0.;
        dispmy = 0.;
        dispmz = 0.;
        maxdisp = 0.;
        for (int i=0; i<n; i++)
        {
            dispmx = (par_fun[i][0]-par_fun_nei[i][0])*(par_fun[i][0]-par_fun_nei[i][0]);
            dispmy = (par_fun[i][1]-par_fun_nei[i][1])*(par_fun[i][1]-par_fun_nei[i][1]);
            dispmz = (par_fun[i][2]-par_fun_nei[i][2])*(par_fun[i][2]-par_fun_nei[i][2]);
            maxdisp = max((dispmx+dispmy+dispmz),maxdisp);
        }

        if ((1.2*maxdisp) > (rskin*rskin))
        {
            n_list(n,par_fun,list,point,xx,yy,zz,rneigh);
            for (int i=0; i<n; i++)        
            {
                par_fun_nei[i][0] = par_fun[i][0];
                par_fun_nei[i][1] = par_fun[i][1];
                par_fun_nei[i][2] = par_fun[i][2];
            }
        }

//-----------------------------------------------------------------------------------------
        acc_cal(n,par_fun,list,point,rcut,eps,sigma,acc,xx,yy,zz,Up_pres);
        vel_ver1(par_fun, acc, dt, n);
        KE=0;
        Temp=0;
        double KE_i=0;
        for (int i=0; i<n; i++)
        {
            acc2[i][0]=acc[i][0];
            acc2[i][1]=acc[i][1];
            acc2[i][2]=acc[i][2];
//            printf("%f\t %f \t %f \n", acc2[i][0],acc2[i][1],acc2[i][2]);
        }    
        KE=KE*mass;
        Temp = 2.0*KE/(3.0*n*Kb);
        acc_cal(n,par_fun,list,point,rcut,eps,sigma,acc,xx,yy,zz,Up_pres);
        vel_ver2(par_fun, acc, acc2, dt, n);
//        printf("%f\t %f \t %f \n", Up_pres[0], KE, KE+Up_pres[0]);

        for (int i=0; i<n; i++)
        {
            KE_i= (par_fun[i][3]*par_fun[i][3]+par_fun[i][4]*par_fun[i][4]+par_fun[i][5]*par_fun[i][5])/2.0;
            KE += KE_i;
        }
        enefile<<fixed<<setprecision(8)<<j<<"\t"<<Up_pres[0]<<"\t"<<Up_pres[1]<<"\t"<<KE<<"\t"<<KE+Up_pres[0]<<"\t"<<Temp<<endl;


        xyzfile <<n<<"\n"<<"\n";
        for (int i=0; i<n; i++)
        {
            pbc(n,par_fun,xx,yy,zz);
            xyzfile <<fixed<<setprecision(8)<<"LJP"<<"\t"<<par_fun[i][0]<<"\t"<<par_fun[i][1]<<"\t"<<par_fun[i][2]<<"\n";
        }
    }

    xyzfile.close();
    enefile.close();

    t2=clock();
    double diff =( ((double)t2 - (double)t1)/ (double) CLOCKS_PER_SEC);
    cout<<diff<<endl;

    return 0;
}





