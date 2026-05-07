#ifndef MY_PROJECT_CHAINMODELMPS_H
#define MY_PROJECT_CHAINMODELMPS_H

#include "itensor/all.h"
#include "OneSigmaSite.h"
#include "functionsMPS.h"
#include <vector>
#include <string>

using namespace itensor;
using std::string;
using std::vector;

class ChainModelMPS {
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
    double E0_;

public:
    ChainModelMPS(double g, int N, bool PBC = false, bool quiet = true);

    // Getters
    int NrSites() const { return NrSites_; }
    double g() const { return g_; }
    bool PBC() const { return PBC_; }

    SiteSet const& sitesVec() const { return sitesVec_; }
    MPS const& p0() const { return p0_; }
    MPO const& H() const { return H_; }
    MPS const& groundState() const { return groundstate_; }
    Real groundEnergy() const { return E0_; }

    // Experiment algorithm loops
    void basicModelLoop(const AlgoLoopParams& params) const;

private:
    // Build MPS, MPO and find groundstate(energy)
    void buildProductStateMPS();
    void buildHamiltonianMPO();
    void checkGroundstateMPS() const;
    double groundStateDMRG(int nsweeps = 8, int maxdim_last = 200, double cutoff = 1E-10, bool quiet = true);

    // Gate applications
    MPS applyA(MPS psi, double theta) const;
    MPS applyAd(MPS psi, double theta) const;

    MPS applyR0(MPS const& psi, double theta) const;
    MPS applyR0dag(MPS const& psi, double theta) const;

    MPS applyU0(MPS psi) const;
    MPS applyU0dag(MPS psi) const;

    MPS applyU(int k, MPS psi, vector<double> const& theta_history) const;
    MPS applyUdag(int k, MPS psi, vector<double> const& theta_history) const;




};

#endif
