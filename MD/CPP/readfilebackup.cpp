#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;

int readfile(std::string inputfile, int n, double par_fun[][6])
{
    ifstream input(inputfile);
    if(input.is_open())
    {
        int count=n+10;

        string length[count];
        int linenumber =0;
        while(!input.eof())
        {
            string number;
            string aaa;
            getline(input,number);

            aaa=number.c_str();
//            cout<<aaa<<endl;
            length[linenumber]=aaa;
            linenumber+=1;
        }

        int nn = atof(length[1].c_str());

        size_t pos=0;
        string delimiter=" ";
        string token;
//        double par_fun[n][6];
        for (int i=2; i<(n+2);i++)
        {
//            cout<<"!!!"<<endl;
            int num=0;
            int j = (i-2);
            int counting=0;
            while((pos=length[i].find(delimiter)) != string::npos)
            {
                token=length[i].substr(0,pos);
                cout<<"!!!"<<endl;
                length[i].erase(0,pos+delimiter.length());
                if (num>3)
                {
                    par_fun[j][counting]=stod(token);
                    counting +=1;
                }
                num+=1;
            }
            par_fun[j][5]=stod(length[i]);
        }
        for (int i=0; i<n;i++)
        {
//            cout<<par_fun[i][0]<<"\t"<< par_fun[i][1]<<"\t"<<par_fun[i][2]<<"\t"<<par_fun[i][3]<<"\t"<<par_fun[i][4]<<"\t"<<par_fun[i][5]<<"\t"<<endl;
        }
    }



}
