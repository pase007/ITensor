//
// Created by Pascal Knoll on 15.03.26.
//

#ifndef MY_PROJECT_TWOQUBITTESTMODEL_H
#define MY_PROJECT_TWOQUBITTESTMODEL_H
#include "itensor/all.h"
#include "functions.h"
#include "data.h"
#include "StructFile.h"
#include "BuildingHamiltonians.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;
using itensor::Index;
using itensor::ITensor;
using cplx = complex<double>;


class TwoQubitTestModel {
private:
    Index s1_;
    Index s2_;

    ITensor H_;

    ITensor ket0_;
    ITensor phi_;

public:
    TwoQubitTestModel(OneQubitOp op1, OneQubitOp op2);

    // Getters
    const Index& s1() const { return s1_; }
    const Index& s2() const { return s2_; }

    const ITensor& H() const { return H_; }

    ITensor const& phi() const { return phi_; }
    ITensor const& ket0() const { return ket0_; }

    // Optional helper checks
    void buildHamiltonian(ITensor const& H1, ITensor const& H2);
    void checkHamiltonian() const;
    void printSummary() const;
    void basicModelLoop(const AlgoLoopParams& params) const;
    void optimizeStepsLoop(const AlgoOptParams& params) const;
};
#endif //MY_PROJECT_TWOQUBITTESTMODEL_H