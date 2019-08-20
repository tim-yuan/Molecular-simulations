#ifndef MDCPP_H
#define MDCPP_H

void initpos(int n, double pos_vel[][6], double rho);
void initvel(int n, double pos_vel[][6], double v_std);
void output(int n, double pos_vel[][6],double x,double y,double z);
void readfile(std::string inputfile, int n, double par_fun[][6],double dim[3]);
void n_list(int n, double par_fun[][6], int list[], int point[],double xx, double yy, double zz, double rneigh);
void acc_cal(int n,double par_fun[][6],int list[],int point[],double rcut, double eps, double sigma,double acc[][3], double xx,double yy, double zz,  double Up_pres[2]);
void vel_ver1(double par_fun[][6],double acc[][3],double dt, int n);
void vel_ver2(double par_fun[][6],double acc[][3],double acc2[][3], double dt, int n);
void pbc(int n, double par_fun[][6],double xx, double yy, double zz );






#endif
