//
// Created by Pascal Knoll on 21.02.26.
//

#ifndef MY_PROJECT_FUNCTIONS_H
#define MY_PROJECT_FUNCTIONS_H
#include "itensor/all_mps.h"
#include <complex>
#include <array>
#include <cmath>
using itensor::prime;
using itensor::noPrime;
using itensor::dag;
using itensor::eltC;
using itensor::ITensor;
using itensor::Index;
using itensor::replaceInds;
using itensor::IndexSet;
using itensor::swapTags;
using namespace std;
using cplx = complex<double>;


// --- Build Itensor from Matrix ---
ITensor opFromArray(Index const& s, cplx a00, cplx a01, cplx a10, cplx a11);
//Itensor opFromDims();

// ---- Build ITensor for Gates ----
ITensor idGate(Index const& s);
ITensor phaseOnZero(Index const& s, double theta);
ITensor phaseOnState(ITensor const& psi_ref, Index const& s, double theta);
ITensor exp_i_theta_H(ITensor const& H, double theta);

ITensor removeTrivial(ITensor T);
// ---  Build Itensor Operations ---
ITensor applyGate(ITensor const& G, ITensor psi);
ITensor composeGate(ITensor V, ITensor U, Index const& s);
ITensor adjointGate(ITensor const& U);
double unitary_Defect(ITensor const& U, Index const& s);
ITensor projectorOnState(ITensor const& psi, Index const& s);
ITensor canonGate(ITensor G, Index const& s);
ITensor reunitarize_polar_gate(ITensor const& U, Index const& s, double eps = 1E-10);
ITensor reunitarize_polar_svd(ITensor const& U, Index const& s);

// --- Build functions for measuring Observables ---
double expectation(ITensor const& psi, ITensor const& Op_sp_s, Index const& s);
double fidelity(ITensor const& psi, ITensor const& phi);





#endif //MY_PROJECT_FUNCTIONS_H