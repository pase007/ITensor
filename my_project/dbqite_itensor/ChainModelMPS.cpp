#include "ChainModelMPS.h"


using namespace itensor;

ChainModelMPS::ChainModelMPS(double g, int N, bool PBC, bool quiet)
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
void ChainModelMPS::checkGroundstateMPS() const {
    double Echeck = real(innerC(groundstate_, H_, groundstate_));
    double normcheck = real(innerC(groundstate_, groundstate_));

    cout << "DMRG energy H  = " << E0_ << "\n";
    cout << "Energy check = " << Echeck << "\n";
    cout << "Norm = " << normcheck << "\n";
}



// ---------------------------- Build MPS |000...> as start state -----------------------------------
void ChainModelMPS::buildProductStateMPS(){
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
void ChainModelMPS::buildHamiltonianMPO(){
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



// -------------------------------- Exp Hamiltonian Trotterization --------------------------------
MPS ChainModelMPS::applyA(MPS psi, double theta) const{
    auto args = Args("Cutoff",1E-12, "MaxDim",200);

    psi = applyTrotterLayerToMPS(psi, models_, theta/2.0, true, NrSites_, g_, args);  // odd half step
    psi.orthogonalize();
    psi = applyTrotterLayerToMPS(psi, models_, theta,    false, NrSites_, g_, args);  // even full step
    psi.orthogonalize();
    psi = applyTrotterLayerToMPS(psi, models_, theta/2.0, true, NrSites_, g_, args);  // odd half step

    psi.orthogonalize();

    return psi;
}

MPS ChainModelMPS::applyAd(MPS psi, double theta) const{
    return applyA(psi, -theta);
}



// ---------------------------------- MPO version of evolution U0 ----------------------------------
MPS ChainModelMPS::applyU0(MPS psi) const{
    return psi;
}

MPS ChainModelMPS::applyU0dag(MPS psi) const{
    return psi;
}



// ------------------------------------ Reflection Gate R0 -------------------------------------------
MPS ChainModelMPS::applyR0(MPS const& psi, double theta) const{
    auto alpha = std::exp(Cplx(0.0, theta)) - Cplx(1.0, 0.0);
    auto overlap = innerC(p0_, psi);
    auto correction = p0_;
    correction *= alpha * overlap;

    auto out = sum(psi, correction, {"Cutoff", 1E-12, "MaxDim", 200});
    return out;
}

MPS ChainModelMPS::applyR0dag(MPS const& psi, double theta) const{
    return applyR0(psi, -theta);
}



// -------------------------- Apply U/Udag string as functions recursively -------------------------------
MPS ChainModelMPS::applyU(int k, MPS psi, vector<double> const& theta_history) const{
    if(k == 0)
        return applyU0(psi);

    double theta = theta_history.at(k-1);

    // U_k = A U_{k-1} R0 U_{k-1}† A† U_{k-1}
    // Rightmost operator acts first.
    psi = applyU(k-1, psi, theta_history);
    psi = applyAd(psi, theta);
    psi = applyUdag(k-1, psi, theta_history);
    psi = applyR0(psi, theta);
    psi = applyU(k-1, psi, theta_history);
    psi = applyA(psi, theta);
    return psi;
}

MPS ChainModelMPS::applyUdag(int k, MPS psi, vector<double> const& theta_history) const{
    if(k == 0)
        return applyU0dag(psi);

    double theta = theta_history.at(k-1);

    // U_k† = U_{k-1}† A U_{k-1} R0† U_{k-1}† A†
    // Rightmost operator acts first.
    psi = applyAd(psi, theta);
    psi = applyUdag(k-1, psi, theta_history);
    psi = applyR0dag(psi, theta);
    psi = applyU(k-1, psi, theta_history);
    psi = applyA(psi, theta);
    psi = applyUdag(k-1, psi, theta_history);

    return psi;
}



// ----------------------------------- DMRG for Groundstate ---------------------------------------
double ChainModelMPS::groundStateDMRG(int nsweeps, int maxdim_last, double cutoff, bool quiet){
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





// ----------------------------------- Algorithm Section ---------------------------------------
void ChainModelMPS::basicModelLoop(const AlgoLoopParams& params) const {
    const double s_step = params.s_step;
    double theta = sqrt(s_step);
    const string s_string = params.s_string;
    const int K = params.K;
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;
    cout << "Step s = " << s_step << "\n";

    //Data container
    vector<Row> rows;
    rows.reserve(K+1);
    vector<double> theta_history;
    theta_history.reserve(K);
    double schedule_factor = 0.85;
    double Var, Ek, Fk, IFk;
    double nrm;

    // Start psi as U0 * p0_
    MPS psi = applyU0(p0_);

    // Norm checks
    cout << "ground norm = " << real(innerC(groundstate_, groundstate_)) << "\n";

    cout << "k\tEnergy(<X>)\t\tVariance(H)\t\tInfidelity(|->)\t\tNorm(psi)\t\tmaxBondDim(psi)\n";
    for (int k = 1; k <= K; ++k) {
        // First do Measurements
        Ek = real(innerC(psi, H_, psi));
        Var = real(innerC(psi, H_, H_, psi)) - Ek*Ek;
        Fk = norm(innerC(groundstate_, psi));
        IFk = 1 - Fk;
        nrm = real(innerC(psi,psi));
        if (nrm > 1.5) {
            cout << "Norm exploded! nrm(psi) = " << nrm << "\n";
            break;
        }
        rows.push_back(Row{k-1, Ek, Fk});
        cout << k-1 << "\t" << Ek << "\t" << Var << "\t" << IFk << "\t" << nrm << "\t" << maxLinkDim(psi) << "\n";

        if (IFk < infid_target) {
            cout << "\nInfidelity of " << str_infid_target << " after " << k << " steps achieved!" << endl;
            break;
        }
        if (k==K) break;

        // Then start doing the evolution
        double shed = theta * pow(sqrt(schedule_factor), k-1);
        theta_history.push_back(shed);
        MPS psiNext = psi;
        psiNext = applyAd(psiNext, shed);
        psiNext = applyUdag(k-1, psiNext, theta_history);
        psiNext = applyR0(psiNext, shed);
        psiNext = applyU(k-1, psiNext, theta_history);
        psiNext = applyA(psiNext, shed);
        psi = psiNext;
    }
    // write data
    write_csv("data" + s_string + ".csv", rows);
}
