//
// Created by Pascal Knoll on 05.05.26.
//

#ifndef MY_PROJECT_FUNCTIONSMPS_H
#define MY_PROJECT_FUNCTIONSMPS_H
#include "itensor/all.h"
#include "OneSigmaSite.h"
#include "BuildingHamiltonians.h"
#include <vector>
#include <string>
using namespace itensor;

// Calculations with MPOs
MPO composeMPO(MPO const& A, MPO const& B);
MPO dagMPO(MPO const& U);
MPO buildIdentityMPO(SiteSet const& sitesVec, int NrSites);
MPO addMPO(MPO const& A, MPO const& B, Cplx a, Cplx b);

// Unitarization checks and fix
double unitarityDefectOnState(MPO const& U, MPS const& psi);
MPO reunitarizeMPO(MPO const& U, SiteSet const& sitesVec, int NrSites);
void testTwoSiteGate(SiteSet const& sitesVec, vector<OneSigmaSite> const& models, int i, double theta, double g);

// ----
void normalizeMPOPrimeLevels(MPO& K);
int maxBondDim(MPO const& M);

// Helper for product State
MPS buildProductStateHelper(SiteSet const& sitesVec, int NrSites, IndexSet const& links);

// Helper for WarmUpState
MPO buildWarmUpHelper(SiteSet const& sitesVec, int NrSites, std::vector<Index> const& links);

// Helper for R0 Gate
// ---

// Helper functions for building trotterization of A
pair<ITensor,ITensor> gateToMPOBond(ITensor G, int i, string tag);
ITensor bondHamiltonian(vector<OneSigmaSite> const& models, int i, double g);
ITensor twoSiteGate(vector<OneSigmaSite> const& models, int i, double theta, double g);
MPS applyTrotterLayerToMPS(MPS psi, vector<OneSigmaSite> const& models, double theta, bool oddLayer, int NrSites, double g, Args const& args);

// Helper functions to build MPO
ITensor makeBulkMPO_(Index const& wL, Index const& wR, OneSigmaSite const& site, double coupling);
ITensor makeFirstMPO_(Index const& wR, OneSigmaSite const& site, double coupling);
ITensor makeLastMPO_(Index const& wL, OneSigmaSite const& site);

#endif //MY_PROJECT_FUNCTIONSMPS_H