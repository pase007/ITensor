//
// Created by Pascal Knoll on 26.03.26.
//

#ifndef MY_PROJECT_BASEMODEL2_H
#define MY_PROJECT_BASEMODEL2_H
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
#include <limits>
#include <utility>

using namespace std;
using itensor::Index;
using itensor::ITensor;
using itensor::combiner;
using itensor::replaceTags;
using cplx = complex<double>;


class BaseModel2 {
protected:
    Index s1_;
    Index s2_;
    IndexSet s12_ = {s1_, s2_};

    ITensor H_;
    ITensor U0_;

    ITensor psi1_;
    ITensor psi2_;
    ITensor psi0_;
    double E0_;

    ITensor ket01_;
    ITensor ket02_;
    ITensor ket0_;
    vector<ITensor> groundspace_;

    ITensor p0_;

    // Protected constructor: only derived classes can construct
    BaseModel2(Index s1, Index s2)
        : s1_(s1), s2_(s2) {}

public:
    virtual ~BaseModel2() = default;

    // Getters
    const Index& s1() const { return s1_; }
    const Index& s2() const { return s2_; }
    const IndexSet& s12() const { return s12_; }

    const ITensor& H() const { return H_; }
    const ITensor& U0() const { return U0_; }

    ITensor const& psi1() const { return psi1_; }
    ITensor const& psi2() const { return psi2_; }
    ITensor const& psi0() const { return psi0_; }

    ITensor const& ket01() const { return ket01_; }
    ITensor const& ket02() const { return ket02_; }
    ITensor const& ket0() const { return ket0_; }

    // Optional helper checks
    void buildHamiltonian(ITensor const& H1, ITensor const& H2);
    void checkHamiltonian() const;

    pair<double, vector<ITensor>> groundSpace2(double tol) const;
    void checkGroundStateDegeneracy(double tol) const;

    void printSummary() const;
    void print_state() const;

    ITensor transformToMatrix(ITensor const& A) const;
    ITensor transformToVector(ITensor const& psi) const;

    void basicModelLoop(const AlgoLoopParams& params) const;
    void optimizeStepsLoop(const AlgoOptParams& params) const;
    void degradingInfidLoop(const AlgoInfidParams& params) const;

};
#endif //MY_PROJECT_BASEMODEL2_H