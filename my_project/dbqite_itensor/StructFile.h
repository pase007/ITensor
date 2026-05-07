//
// Created by Pascal Knoll on 15.03.26.
//

#ifndef MY_PROJECT_STRUCTFILE_H
#define MY_PROJECT_STRUCTFILE_H
#include "itensor/all.h"
#include <iostream>
#include <iomanip>
using namespace std;
using cplx = complex<double>;

// Algorithm Parameters
struct AlgoLoopParams{
    double s_step;
    string s_string;
    int K;
    double infid_target;
    string str_infid_target;
};

struct AlgoOptParams{
    double s_min;
    double s_bin;
    double s_N;
    int K_max;
    double infid_target;
    string str_infid_target;
};

struct AlgoOptParamsVec{
    double s_min;
    double s_bin;
    double s_N;
    int K_max;
    vector <double> infid_target;
    vector <string> str_infid_target;
};

struct AlgoInfidParams{
    double s_step;
    string s_string;
    int K_max;
    double infid_min;
    int infid_N;
};


//Helper for Hamiltonians
struct PauliTerm{
    cplx coeff;
    char left;
    char right;
};


// Data Containers
struct Row{
    int k;
    double energy;
    double fidelity;
};

struct RowOpt {
    double s_step;
    int k;
    double infidelity;
};

struct RowInfid {
    double infid_step;
    int k;
    double infidelity;
};



#endif //MY_PROJECT_STRUCTFILE_H