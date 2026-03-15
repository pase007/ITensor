//
// Created by Pascal Knoll on 14.03.26.
//

#ifndef MY_PROJECT_TESTMODEL_H
#define MY_PROJECT_TESTMODEL_H
#include "itensor/all.h"
#include "functions.h"
#include "data.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;
using itensor::Index;
using itensor::ITensor;
using cplx = complex<double>;

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

class TestModel {
private:
    Index s_;

    ITensor H_;
    ITensor H2_;
    ITensor H3_;

    ITensor ket0_;
    ITensor phi_;

public:
    TestModel();

    // Getters
    ITensor const& H() const { return H_; }
    ITensor const& H2() const { return H2_; }
    ITensor const& H3() const { return H3_; }

    Index const& s() const { return s_; }
    ITensor const& phi() const { return phi_; }
    ITensor const& ket0() const { return ket0_; }

    // Optional helper checks
    void checkHamiltonian() const;
    void printSummary() const;
    void basicModelLoop(const AlgoLoopParams& params) const;
    void optimizeStepsLoop(const AlgoOptParams& params) const;
};
#endif //MY_PROJECT_TESTMODEL_H