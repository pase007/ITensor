//
// Created by Pascal Knoll on 07.05.26.
//
#include "SigmaChainAnalysis.h"

using namespace itensor;

SigmaChainAnalysis::SigmaChainAnalysis(double g, int N, bool PBC, bool quiet)
    : NrSites_(N), g_(g), PBC_(PBC), p0_(N), H_(N), U0_(N), groundstate_(N){

    // Build chain of sites and models
    vector<Index> temp_sites;
    for(int i = 0; i < NrSites_; ++i){
        models_.emplace_back(g_);
        temp_sites.push_back(models_.back().s());
    }
    sitesVec_ = SiteSet(temp_sites);

    // Build MPS and H-MPO
    buildProductStateMPS();
    buildHamiltonianMPO();

    // Groundstate search with DMRG
    E0_ = groundStateDMRG(8, 200, 1E-10, quiet);
    checkGroundstateMPS();
}



// ------------------------------------ Groundstate/E0 function --------------------------------------
void SigmaChainAnalysis::checkGroundstateMPS() const {
    double Echeck = real(innerC(groundstate_, H_, groundstate_));
    double normcheck = real(innerC(groundstate_, groundstate_));

    cout << "DMRG energy H  = " << E0_ << "\n";
    cout << "Energy check = " << Echeck << "\n";
    cout << "Norm = " << normcheck << "\n";
}



// ---------------------------- Build MPS |000...> as start state -----------------------------------
void SigmaChainAnalysis::buildProductStateMPS(){
    vector<Index> templinks;
    for(int i = 1; i < NrSites_; ++i) {
        templinks.push_back(Index(1, "Link,l=" + to_string(i)));
    }
    IndexSet links(templinks);
    p0_ = buildProductStateHelper(sitesVec_, NrSites_, links);
    p0_.orthogonalize();
    p0_.normalize();
}



// ----------------------------- Build MPO version of the Hamiltonian ----------------------------
// --- (DMRGtensor.pdf site 60) ---
void SigmaChainAnalysis::buildHamiltonianMPO(){
    double coupling = -3.0 / (4.0 * g_ * g_);

    // Link indicies vector, size 5
    std::vector<Index> links;
    links.reserve(NrSites_ - 1);
    for(int i = 1; i <= NrSites_ - 1; ++i){
        links.emplace_back(5, "Link,w=" + std::to_string(i));
    }

    // Build Hamiltonian H as an MPO - here closed boundary conditions
    H_.setA(1, makeFirstMPO_(links.at(0), models_.at(0), coupling));
    for(int i = 2; i <= NrSites_ - 1; ++i){
        H_.setA(i, makeBulkMPO_(links.at(i-2), links.at(i-1), models_.at(i-1), coupling));
    }
    H_.setA(NrSites_, makeLastMPO_(links.back(), models_.back()));
}



// ----------------------------------- DMRG for Groundstate ---------------------------------------
double SigmaChainAnalysis::groundStateDMRG(int nsweeps, int maxdim_last, double cutoff, bool quiet){
    auto sweeps = Sweeps(nsweeps);
    sweeps.maxdim() = 20, 40, 80, maxdim_last, maxdim_last, maxdim_last;
    sweeps.cutoff() = cutoff;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7, 1E-8, 0.0;

    auto psi = p0_;
    auto [energy, psiGS] = dmrg(H_, psi, sweeps, {"Silent", quiet});

    E0_ = energy;
    groundstate_ = psiGS;
    groundstate_.normalize();

    return E0_;
}



// ----------------------------------- Add Groundstate to H ---------------------------------------





// ----------------------------------- DMRG for Excited State ---------------------------------------
