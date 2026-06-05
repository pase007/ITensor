//
// Created by Pascal Knoll on 07.05.26.
//

#ifndef MY_PROJECT_SIGMACHAINANALYSIS_H
#define MY_PROJECT_SIGMACHAINANALYSIS_H
#include "itensor/all.h"
#include "OneSigmaSite.h"
#include "functionsMPS.h"
#include <vector>
#include <string>
#include <utility> // ADDED correlation-length helpers

using namespace itensor;
using std::string;
using std::vector;

class SigmaChainAnalysis {
protected:
    int NrSites_;
    double g_;
    bool PBC_;

    vector<OneSigmaSite> models_;
    SiteSet sitesVec_;

    // MPS / MPO objects
    MPS p0_;
    MPO H_;
    MPO U0_;

    MPS groundstate_;
    MPS excState_;
    double E0_;
    double E1_;
    double gap_;
    double var0_;
    double var1_;
    vector<double> energiesVec_;

public:
    SigmaChainAnalysis(double g, int N, bool PBC = false, bool quiet = true);

    // Getters
    int NrSites() const { return NrSites_; }
    double g() const { return g_; }
    bool PBC() const { return PBC_; }

    SiteSet const& sitesVec() const { return sitesVec_; }
    MPS const& p0() const { return p0_; }
    MPO const& H() const { return H_; }
    MPS const& groundState() const { return groundstate_; }
    MPS const& excState() const { return excState_; }
    Real groundEnergy() const { return E0_; }
    Real excEnergy() const { return E1_; }
    Real gap() const { return gap_; }
    Real groundVariance() const { return var0_; }
    Real excVariance() const { return var1_; }
    vector<double> const& energiesVec() const { return energiesVec_; }

    // --- ADDED correlation-length helpers: public analysis interface ---
    ITensor localOp(int site, int comp) const;
    Cplx onePointC(int site, int comp) const;
    Cplx twoPointC(int i, int j, int comp_i, int comp_j) const;
    double onePoint(int site, int comp) const;
    double twoPoint(int i, int j, int comp_i, int comp_j) const;
    double connectedCorr(int i, int j, int comp_i, int comp_j) const;
    vector<std::pair<int,double>> correlationProfile(int i0, int comp_i, int comp_j,
                                                     int rMin = 1, int rMax = -1) const;
    double fitCorrelationLength(int i0, int comp_i, int comp_j,
                                int rMin = 2, int rMax = -1) const;
    // --- END ADDED correlation-length helpers ---

    // Experiment algorithm loops

private:
    // Build MPS, MPO
    void buildProductStateMPS();
    void buildHamiltonianMPO();

    // State Checks
    void checkGroundstateMPS() const;
    void checkExcitedState() const;
    double energyVariance(MPS const& psi, double energy) const;

    // Sweep Functions
    double groundStateDMRG(int nsweeps = 8, int maxdim_last = 200, double cutoff = 1E-10, bool quiet = true);
    double excStateDMRG(int nsweeps = 8, int maxdim_last = 200, double cutoff = 1E-10, bool quiet = true, double c = 5.0);

    // Energy analysis
    double energyDensityConvergenceTest(double g, int N) const;
    double massGapScan(double g, int N);


    // Gate applications







};

#endif //MY_PROJECT_SIGMACHAINANALYSIS_H
