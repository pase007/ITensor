#include "ChainModelMPS.h"
#include <future>
#include <limits>
#include <utility>


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

MPS ChainModelMPS::evolveOneStep(MPS psi, int current_level, double theta_step, vector<double> const& theta_history) const{
    psi = applyAd(psi, theta_step);
    psi = applyUdag(current_level, psi, theta_history);
    psi = applyR0(psi, theta_step);
    psi = applyU(current_level, psi, theta_history);
    psi = applyA(psi, theta_step);
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
    groundstate_ = std::move(psiGS);
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
    double schedule_factor = params.schedule_factor;
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
        psi = evolveOneStep(psi, k-1, shed, theta_history);
    }
    // write data
    write_csv("data" + s_string + ".csv", rows);
}

void ChainModelMPS::adaptiveSModelLoop(const AlgoLoopParams& params) const {
    const double s_step = params.s_step;
    const string s_string = params.s_string;
    const int K = params.K;
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;
    cout << "Adaptive initial s = " << s_step << "\n";

    vector<Row> rows;
    rows.reserve(K+1);
    vector<double> theta_history;
    theta_history.reserve(K);
    double current_s = s_step;

    struct CandidateResult {
        bool stable;
        double s;
        double theta;
        double energy;
        double norm;
        MPS psi;
    };
    auto evaluateCandidates = [&](vector<double> const& candidates,
                                  MPS const& current_psi,
                                  int current_level,
                                  vector<double> const& current_history,
                                  double max_s) {
        vector<std::future<CandidateResult>> trials;
        trials.reserve(candidates.size());
        for (double candidate_s : candidates) {
            if (candidate_s <= 0.0) continue;
            if (candidate_s > max_s) continue;

            trials.push_back(std::async(std::launch::async, [this, candidate_s, current_history, current_psi, current_level]() {
                double candidate_theta = sqrt(candidate_s);
                vector<double> trial_history = current_history;
                trial_history.push_back(candidate_theta);

                MPS trial_psi = evolveOneStep(current_psi, current_level, candidate_theta, trial_history);
                double trial_norm = real(innerC(trial_psi, trial_psi));
                if (!std::isfinite(trial_norm) || trial_norm > 1.5) {
                    return CandidateResult{false, candidate_s, candidate_theta,
                                           std::numeric_limits<double>::infinity(),
                                           trial_norm, std::move(trial_psi)};
                }

                double trial_energy = real(innerC(trial_psi, H_, trial_psi));
                return CandidateResult{true, candidate_s, candidate_theta,
                                       trial_energy, trial_norm, std::move(trial_psi)};
            }));
        }

        bool found_candidate = false;
        double best_energy = std::numeric_limits<double>::infinity();
        CandidateResult best_result{false, 0.0, 0.0, best_energy, 0.0, MPS()};
        for (auto& trial : trials) {
            CandidateResult result = trial.get();
            if (!result.stable) continue;
            if (result.energy < best_energy) {
                found_candidate = true;
                best_energy = result.energy;
                best_result = std::move(result);
            }
        }

        if (!found_candidate) {
            return CandidateResult{false, 0.0, 0.0,
                                   std::numeric_limits<double>::infinity(),
                                   0.0, MPS()};
        }
        return best_result;
    };

    MPS psi = applyU0(p0_);

    cout << "ground norm = " << real(innerC(groundstate_, groundstate_)) << "\n";
    cout << "k\tEnergy(<X>)\t\tVariance(H)\t\tInfidelity(|->)\t\tNorm(psi)\t\tmaxBondDim(psi)\n";
    for (int k = 1; k <= K; ++k) {
        double Ek = real(innerC(psi, H_, psi));
        double Var = real(innerC(psi, H_, H_, psi)) - Ek*Ek;
        double Fk = norm(innerC(groundstate_, psi));
        double IFk = 1 - Fk;
        double nrm = real(innerC(psi,psi));
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

        if (k == 1) {
            double initial_theta = sqrt(s_step);
            theta_history.push_back(initial_theta);
            cout << "initial s_0 = " << s_step
                 << " (theta = " << initial_theta << ")\n";
            psi = evolveOneStep(psi, k-1, initial_theta, theta_history);
            continue;
        }

        vector<double> candidates = params.s_candidates;
        if (candidates.empty()) {
            candidates = {0.5*s_step, 0.6*s_step, 0.7*current_s, 0.8*current_s, 0.9*current_s, current_s};
        }

        CandidateResult coarse_result = evaluateCandidates(candidates, psi, k-1, theta_history, current_s);
        if (!coarse_result.stable) {
            cout << "No stable adaptive s candidate found at k = " << k << "\n";
            break;
        }

        double coarse_s = coarse_result.s;
        CandidateResult best_result;
        if (params.refine_s) {
            vector<double> refine_candidates = {0.8*coarse_s, 0.9*coarse_s, coarse_s,
                                                1.1*coarse_s, 1.2*coarse_s};
            CandidateResult refine_result = evaluateCandidates(refine_candidates, psi, k-1, theta_history, current_s);
            if (refine_result.stable) {
                best_result = std::move(refine_result);
            } else {
                best_result = std::move(coarse_result);
            }
        } else {
            best_result = std::move(coarse_result);
        }

        cout << "adaptive s_" << k-1 << " = " << best_result.s
             << " (theta = " << best_result.theta
             << ", trial E = " << best_result.energy
             << ", dE = " << best_result.energy - Ek
             << ", coarse s = " << coarse_s
             << ", refined = " << (params.refine_s ? "yes" : "no") << ")\n";

        theta_history.push_back(best_result.theta);
        current_s = best_result.s;
        psi = std::move(best_result.psi);
    }

    write_csv("data" + s_string + ".csv", rows);
}
