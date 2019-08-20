#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <stdio.h>
using namespace std;

//int readfile(std::string inputfile, int n, double par_fun[][6])
//{
//    ifstream input(inputfile);
//    if(input.is_open())
//    {
////        input.seekg(23);
////        string s;
////        s.resize(68 - 23);
////        input.read(&s[0], 68 - 23);
////        cout<<s<<endl;
////    }
////    return 0;
////}
////
//
//        int count=n+10;
//
//        string length[count];
//        int linenumber =0;
//        while(!input.eof())
//        {
//            string number;
//            string aaa;
//            getline(input,number);
//
//            aaa=number.c_str();
//            cout<<aaa<<endl;
//            length[linenumber]=aaa;
//            linenumber+=1;
//        }
//
//        int n_mem = atof(length[1].c_str());
//
//        size_t pos=0;
//        string delimiter=" ";
//        string token;
////        double par_fun[n][6];
//        for (int i=2; i<(n+2);i++)
//        {
////            cout<<"!!!"<<endl;
//            int num=0;
//            int j = (i-2);
//            int counting=0;
//            while((pos=length[i].find(delimiter)) != string::npos)
//            {
//                token=length[i].substr(0,pos);
//                cout<<"!!!"<<endl;
//                cout<<token<<endl;
//                length[i].erase(0,pos+delimiter.length());
////                if (num>3)
////                {
////                    par_fun[j][counting]=stod(token);
////                    counting +=1;
////                }
//                num+=1;
//            }
//            par_fun[j][5]=stod(length[i]);
//        }
//        for (int i=0; i<n;i++)
//        {
////            cout<<par_fun[i][0]<<"\t"<< par_fun[i][1]<<"\t"<<par_fun[i][2]<<"\t"<<par_fun[i][3]<<"\t"<<par_fun[i][4]<<"\t"<<par_fun[i][5]<<"\t"<<endl;
//        }
//    }
//}

int readfile(std::string inputfile, int n, double par_fun[][6], double dim[3])
{                                                              
//    ifstream input(inputfile);
//    if(input.is_open())
    FILE * input;
    input = fopen (inputfile.c_str(), "r");

    char mystring [100];
//    int count = 0;
    int aa, bb;
    char cc, dd;
    int * aaa = &aa;
    int * bbb = &bb;
    char * ccc = &cc;
    char * ddd = &dd;
    float i_mem, j_mem, k_mem, l_mem, m_mem, n_mem;
    float * i_poi =&i_mem; 
    float * j_poi =&j_mem;
    float * k_poi =&k_mem;
    float * l_poi =&l_mem;
    float * m_poi =&m_mem;
    float * n_poi =&n_mem;
    float x_mem, y_mem, z_mem;
    float * x_poi=&x_mem;    
    float * y_poi=&y_mem;    
    float * z_poi=&z_mem;
    fgets (mystring, 100, input);
    fgets (mystring, 100, input);
    for (int i=0; i<(n); i++)
    {
        fscanf (input, "%5d%5s%5s%5d%8f%8f%8f%8f%8f%8f", aaa, ccc, ddd, bbb,i_poi, j_poi, k_poi, l_poi,m_poi, n_poi);
//        cout<<i_mem<<"\t"<<j_mem<<"\t"<<k_mem<<"\t"<<l_mem<<"\t"<<m_mem<<"\t"<<n_mem<<endl;
        par_fun[i][0]=i_mem;
        par_fun[i][1]=j_mem;
        par_fun[i][2]=k_mem;
        par_fun[i][3]=l_mem;
        par_fun[i][4]=m_mem;
        par_fun[i][5]=n_mem;
    }
    fscanf (input, "%f", x_poi);
    fscanf (input, "%f", y_poi);
    fscanf (input, "%f", z_poi);
    dim[0]=x_mem;
    dim[1]=y_mem;
    dim[2]=z_mem;
//    cout<<x_mem<<" "<<y_mem<<" "<<z_mem<<endl;

    fclose(input);
    return 0;
}




