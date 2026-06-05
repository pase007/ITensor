//
// Created by Pascal Knoll on 07.05.26.
//
#include "SigmaChainAnalysis.h"
#include <cmath> // ADDED correlation-length helpers
#include <utility>

using namespace itensor;

SigmaChainAnalysis::SigmaChainAnalysis(double g, int N, bool PBC, bool quiet)
    : NrSites_(N), g_(g), PBC_(PBC), p0_(N), H_(N), U0_(N), groundstate_(N), excState_(N),
      E0_(0.0), E1_(0.0), gap_(0.0), var0_(0.0), var1_(0.0){

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
    E0_ = groundStateDMRG(14, 250, 1E-10, quiet);
    checkGroundstateMPS();
    E1_ = excStateDMRG(27, 250, 1E-10, quiet, 7.0);
    gap_ = E1_ - E0_;
    checkExcitedState();
}



// ------------------------------------ Groundstate/E0 function --------------------------------------
void SigmaChainAnalysis::checkGroundstateMPS() const {
    double Echeck = real(innerC(groundstate_, H_, groundstate_));
    double normcheck = real(innerC(groundstate_, groundstate_));

    cout << "DMRG energy H  = " << E0_ << "\n";
    cout << "Energy check = " << Echeck << "\n";
    cout << "Energy variance = " << var0_ << "\n";
    cout << "Norm = " << normcheck << "\n\n" ;
}

void SigmaChainAnalysis::checkExcitedState() const {
    double Echeck = real(innerC(excState_, H_, excState_));
    double normcheck = real(innerC(excState_, excState_));
    double overlap = std::abs(innerC(excState_, groundstate_));

    cout << "DMRG exc energy H  = " << E1_ << "\n";
    cout << "Energy check = " << Echeck << "\n";
    cout << "Energy variance = " << var1_ << "\n";
    cout << "Norm = " << normcheck << "\n";
    cout << "GS overlap = " << overlap << "\n";
    cout << "Mass gap E1-E0 = " << gap_ << "\n\n";
}

double SigmaChainAnalysis::energyVariance(MPS const& psi, double energy) const {
    double H2 = real(innerC(H_, psi, H_, psi));
    double variance = H2 - energy * energy;
    return std::max(variance, 0.0);
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
    groundstate_ = std::move(psiGS);
    groundstate_.normalize();
    var0_ = energyVariance(groundstate_, E0_);

    return E0_;
}



// ----------------------------------- DMRG for Excited State ---------------------------------------
double SigmaChainAnalysis::excStateDMRG(int nsweeps, int maxdim_last, double cutoff, bool quiet, double c){
    auto sweeps = Sweeps(nsweeps);
    sweeps.maxdim() = 20, 40, 80, maxdim_last, maxdim_last, maxdim_last;
    sweeps.cutoff() = cutoff;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7, 1E-8, 0.0;

    auto psi = p0_;
    vector<MPS> penaltyStates = {groundstate_};
    auto [energy, psiES] = dmrg(H_, penaltyStates, psi, sweeps, {"Weight", c, "Silent", quiet});

    E1_ = energy;
    excState_ = std::move(psiES);
    excState_.normalize();
    var1_ = energyVariance(excState_, E1_);

    return E1_;
}


// --- ADDED correlation-length helpers: direct MPS contractions without custom SiteSet ops ---
ITensor SigmaChainAnalysis::localOp(int site, int comp) const {
    if(site < 1 || site > NrSites_) {
        itensor::error("localOp: site out of range");
    }
    auto const& m = models_.at(site - 1);
    if(comp == 0) return m.H0(); // ADDED correlation-length helpers: comp 0 is local Laplace-Beltrami kinetic term
    if(comp == 1) return m.j1();
    if(comp == 2) return m.j2();
    if(comp == 3) return m.j3();

    itensor::error("localOp: comp must be 0, 1, 2, or 3");
    return ITensor();
}

Cplx SigmaChainAnalysis::onePointC(int site, int comp) const {
    if(site < 1 || site > NrSites_) {
        itensor::error("onePointC: site out of range");
    }

    auto psi = groundstate_;
    psi.position(site);

    auto val = psi(site)
             * localOp(site, comp)
             * dag(prime(psi(site), sitesVec_(site)));

    return eltC(val);
}

Cplx SigmaChainAnalysis::twoPointC(int i, int j, int comp_i, int comp_j) const {
    if(i < 1 || i > NrSites_ || j < 1 || j > NrSites_) {
        itensor::error("twoPointC: site out of range");
    }
    if(i == j) {
        itensor::error("twoPointC: same-site operator products are not implemented");
    }
    if(i > j) {
        std::swap(i, j);
        std::swap(comp_i, comp_j);
    }

    auto psi = groundstate_;
    psi.position(i);

    ITensor leftEnv(1.0);
    if(i > 1) {
        auto leftLink = commonIndex(psi(i), psi(i - 1));
        leftEnv = delta(dag(leftLink), prime(leftLink));
    }

    auto partial = (leftEnv * psi(i) * localOp(i, comp_i)) * dag(prime(psi(i)));

    for(int site = i + 1; site < j; ++site) {
        partial *= psi(site);
        partial *= dag(prime(psi(site), "Link"));
    }

    auto leftLinkAtJ = commonIndex(psi(j), partial);
    partial *= psi(j);

    auto val = (partial * localOp(j, comp_j))
             * dag(prime(prime(psi(j), "Site"), leftLinkAtJ));

    return eltC(val);
}

double SigmaChainAnalysis::onePoint(int site, int comp) const {
    return real(onePointC(site, comp));
}

double SigmaChainAnalysis::twoPoint(int i, int j, int comp_i, int comp_j) const {
    return real(twoPointC(i, j, comp_i, comp_j));
}

double SigmaChainAnalysis::connectedCorr(int i, int j, int comp_i, int comp_j) const {
    auto cij = twoPointC(i, j, comp_i, comp_j);
    auto ci = onePointC(i, comp_i);
    auto cj = onePointC(j, comp_j);
    return real(cij - ci * cj);
}

vector<std::pair<int,double>>
SigmaChainAnalysis::correlationProfile(int i0, int comp_i, int comp_j, int rMin, int rMax) const {
    if(i0 < 1 || i0 > NrSites_) {
        itensor::error("correlationProfile: i0 out of range");
    }
    if(rMax < 0) {
        rMax = NrSites_ - i0;
    }

    vector<std::pair<int,double>> profile;
    for(int r = rMin; r <= rMax; ++r) {
        int j = i0 + r;
        if(j > NrSites_) break;
        profile.push_back({r, connectedCorr(i0, j, comp_i, comp_j)});
    }
    return profile;
}

double SigmaChainAnalysis::fitCorrelationLength(int i0, int comp_i, int comp_j,
                                                int rMin, int rMax) const {
    auto profile = correlationProfile(i0, comp_i, comp_j, rMin, rMax);

    double sx = 0.0;
    double sy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;
    int n = 0;

    for(auto const& [r, corr] : profile) {
        double absCorr = std::abs(corr);
        if(absCorr <= 1E-14 || !std::isfinite(absCorr)) continue;

        double x = static_cast<double>(r);
        double y = std::log(absCorr);

        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        ++n;
    }

    if(n < 2) {
        itensor::error("fitCorrelationLength: not enough nonzero data points");
    }

    double denom = n * sxx - sx * sx;
    if(std::abs(denom) <= 1E-14) {
        itensor::error("fitCorrelationLength: singular linear fit");
    }

    double slope = (n * sxy - sx * sy) / denom;
    if(slope >= 0.0) {
        itensor::error("fitCorrelationLength: fitted slope is non-negative");
    }

    return -1.0 / slope;
}
// --- END ADDED correlation-length helpers ---



// ----------------------------------- Energy Density convergence ---------------------------------------
double SigmaChainAnalysis::energyDensityConvergenceTest(double g, int N) const {
    vector<double> nrmE0;
    return 1.0;
}
