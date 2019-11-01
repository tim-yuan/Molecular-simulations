#!/bin/bash
module load gcc/5.4.0
#g++ -std=c++11 -g -Wall -o main main.cpp readfile.cpp n_list.cpp initposvel.cpp acc_cal.cpp vel_ver.cpp pbc.cpp
g++ -std=c++11 -Wall -o main main.cpp readfile.cpp n_list.cpp initposvel.cpp acc_cal.cpp vel_ver.cpp pbc.cpp
