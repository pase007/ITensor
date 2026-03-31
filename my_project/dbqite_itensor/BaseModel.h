//
// Created by Pascal Knoll on 25.03.26.
//

#ifndef MY_PROJECT_BASEMODEL_H
#define MY_PROJECT_BASEMODEL_H
#include "itensor/all.h"
#include "functions.h"
#include "data.h"
#include "StructFile.h"
#include "BuildingHamiltonians.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
using namespace std;
using itensor::Index;
using itensor::ITensor;
using cplx = complex<double>;

class BaseModel {
protected:
    Index s_;
    ITensor Id_s_;

    ITensor H_;
    ITensor H0_;
    ITensor j1_;
    ITensor j2_;
    ITensor j3_;

    ITensor ket0_;
    ITensor phi_;

    // Aditional Hamiltonian options
    ITensor H2_;
    ITensor H3_;

    // Protected constructor: only derived classes can construct
    BaseModel(int dim, const string& name)
        : s_(dim, name){}

public:
    virtual ~BaseModel() = default;

    // Getters
    ITensor const& H() const { return H_; }
    ITensor const& H0() const { return H0_; }
    ITensor const& j1() const { return j1_; }
    ITensor const& j2() const { return j2_; }
    ITensor const& j3() const { return j3_; }

    Index const& s() const { return s_; }
    ITensor const& Id_s() const { return Id_s_; }
    ITensor const& phi() const { return phi_; }
    ITensor const& ket0() const { return ket0_; }

    // Shared methods
    void checkHamiltonian() const;
    void printSummary() const;
    void basicModelLoop(const AlgoLoopParams& params) const;
    void optimizeStepsLoop(const AlgoOptParams& params) const;
    void degradingInfidLoop(const AlgoInfidParams& params) const;
};

#endif //MY_PROJECT_BASEMODEL_H