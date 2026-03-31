//
// Created by Pascal Knoll on 30.03.26.
//

#ifndef MY_PROJECT_CHAINMODEL_H
#define MY_PROJECT_CHAINMODEL_H
#include "itensor/all.h"
#include "functions.h"
#include "data.h"
#include "StructFile.h"
#include "BuildingHamiltonians.h"
#include "BaseModel.h"
#include "OneSigmaSite.h"
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
using cplx = complex<double>;


class ChainModel {
protected:
    IndexSet Iset_;
    vector<OneSigmaSite> models_;
    int NrSites_;
    double g_;

    //Full Hamiltonian
    bool PBC_;
    ITensor H_;

    //Start State
    ITensor phi_;

    // States of the system
    vector<ITensor> groundspace_;
    vector<double> Energies_;
    double E0_;

    //Algorithm definitions
    ITensor p0_;
    ITensor U0_;




public:
    ChainModel(double g, int N, bool PBC);

    // Getters
    IndexSet const& Iset() const { return Iset_; }
    ITensor const& H() const { return H_; }
    int const& NrSites() const { return NrSites_; }
    double const& g() const { return g_; }

    ITensor const& phi() const { return phi_; }

    // Shared methods
    void checkHamiltonian() const;
    void printSummary() const;

    void buildHamiltonian();
    pair<double, vector<ITensor>> groundSpace2(double tol) const;
    void checkGroundStateDegeneracy(double tol) const;

    void basicModelLoop(const AlgoLoopParams& params) const;
    void optimizeStepsLoop(const AlgoOptParams& params) const;
    void degradingInfidLoop(const AlgoInfidParams& params) const;
};


#endif //MY_PROJECT_CHAINMODEL_H