//
// Created by Pascal Knoll on 30.03.26.
//

#include "ChainModel.h"

ChainModel::ChainModel(double g, int N, bool PBC) {
    PBC_ = PBC;
    NrSites_ = N;
    g_ = g;

    // Init of ITensors
    H_ = ITensor(0.0);
    U0_ = ITensor(1.0);
    p0_ = ITensor(1.0);

    // Build vector for models and indicies
    vector<Index> Ivec;
    for (int i=0; i<NrSites_; i++) {
        // Models and corresponding IndexSet
        OneSigmaSite sigma_i = OneSigmaSite(g_);
        Index s_i = sigma_i.s();
        Ivec.push_back(s_i);
        models_.push_back(sigma_i);

        // Warm up Gate U_0 and trivial state |00...>
        U0_ *= hadamard4(s_i);
        ITensor p_i = ITensor(s_i);
        p_i.set(s_i(1), 1.0);
        p_i.set(s_i(2), 0.0);
        p_i.set(s_i(3), 0.0);
        p_i.set(s_i(4), 0.0);
        p0_ *= p_i;
    }
    Iset_= IndexSet(Ivec);

    // Build he system's Hamiltonian
    buildHamiltonian();

    // taget: ground state
    pair<double, vector<ITensor>> gs = groundSpace2(1E-14);
    E0_   = gs.first;
    groundspace_ = gs.second;
    checkGroundStateDegeneracy(1E-14);

}





//----------------------------------------- Build Hamiltonian of the system --------------------------------------------
void ChainModel::buildHamiltonian(){
    for (int i=0; i<NrSites_; i++) {
        ITensor H0_local = ITensor(1.0);
        for (int j=0; j<NrSites_; j++) {
            const auto& sigma_j = models_[j];
            if (i==j) {
                H0_local *= sigma_j.H0();
            } else {
                H0_local *= sigma_j.Id_s();
            }
        }
        if(i==0){
            H_ = H0_local;
        }else{
            H_ += H0_local;
        }
    }

    for(int i = 0; i < NrSites_ - 1; ++i){
        ITensor link1 = ITensor(1.0);
        ITensor link2 = ITensor(1.0);
        ITensor link3 = ITensor(1.0);

        for(int j = 0; j < NrSites_; ++j){
            const auto& sigma_j = models_[j];

            if(j == i || j == i+1){
                link1 *= sigma_j.j1();
                link2 *= sigma_j.j2();
                link3 *= sigma_j.j3();
            }else{
                link1 *= sigma_j.Id_s();
                link2 *= sigma_j.Id_s();
                link3 *= sigma_j.Id_s();
            }
        }
        H_ += -3.0/(4.0*pow(g_,2)) * (link1 + link2 + link3);
    }
    if (PBC_ == true) {
        ITensor link1 = ITensor(1.0);
        ITensor link2 = ITensor(1.0);
        ITensor link3 = ITensor(1.0);

        for(int j = 0; j < NrSites_; ++j){
            const auto& sigma_j = models_[j];

            if(j == 0 || j == NrSites_-1){
                link1 *= sigma_j.j1();
                link2 *= sigma_j.j2();
                link3 *= sigma_j.j3();
            }else{
                link1 *= sigma_j.Id_s();
                link2 *= sigma_j.Id_s();
                link3 *= sigma_j.Id_s();
            }
        }
        H_ += -3.0/(4.0*pow(g_,2)) * (link1 + link2 + link3);
    }
}


void ChainModel::printSummary() const{
    cout << "Site indicies: ";
    for(auto const& Inds : inds(H_)){
        cout << Inds <<  "\n";
    }
    //cout << "H = " << H_ << "\n";
    cout << "GS: E = " << E0_ << " and Xi = " << groundspace_[0] << "\n";
    //cout << "start-state: psi0 = " << psi0_ << "\n";
}


//------------------------------------ Diagonalization: Groundstate and Energies ---------------------------------------
pair<double, vector<ITensor>> ChainModel::groundSpace2(double tol) const{
    ITensor D, V;
    diagHermitian(H_, V, D);

    auto dinds = inds(D);
    if(length(dinds) != 2){
        itensor::error("groundSpace2: D does not have exactly 2 indices");
    }

    Index d  = dinds[0];
    Index dp = dinds[1];

    if(dim(d) != dim(dp)){
        itensor::error("groundSpace2: eigenvalue indices of D have different dimensions");
    }

    double Emin = std::numeric_limits<double>::infinity();
    for(int n = 1; n <= dim(d); ++n){
        double lam = real(eltC(D, d(n), dp(n)));
        if(lam < Emin) Emin = lam;
    }

    auto vinds = inds(V);
    Index u;
    bool found_u = false;

    for(auto const& I : vinds){
        if(!hasIndex(Iset_, I)){
            u = I;
            found_u = true;
            break;
        }
    }

    if(!found_u){
        itensor::error("groundSpace2: could not identify eigenvector index in V");
    }

    vector<ITensor> gs_vecs;

    for(int n = 1; n <= dim(d); ++n){
        double lam = real(eltC(D, d(n), dp(n)));

        if(std::abs(lam - Emin) < tol){
            ITensor psi = V * setElt(u(n));
            auto nrm = norm(psi);
            if(nrm > 0.0) psi /= nrm;
            gs_vecs.push_back(psi);
        }
    }

    return {Emin, gs_vecs};
}

void ChainModel::checkGroundStateDegeneracy(double tol) const{
        ITensor D, V;
        diagHermitian(H_, V, D);

        auto dinds = inds(D);
        if(length(dinds) != 2)
            itensor::error("checkSpectrum: D does not have exactly 2 indices");

        Index d  = dinds[0];
        Index dp = dinds[1];

        // Finde minimalen Eigenwert
        double Emin = std::numeric_limits<double>::infinity();
        for(int n = 1; n <= dim(d); ++n){
            double lam = real(eltC(D, d(n), dp(n)));
            if(lam < Emin) Emin = lam;
        }

        std::cout << "================ Spectrum =================\n";

        // Finde Eigenvektor-Index
        auto vinds = inds(V);
        Index u;
        bool found_u = false;

        for(auto const& I : vinds){
           if(!hasIndex(Iset_, I)){
               u = I;
               found_u = true;
               break;
            }
        }
        if(!found_u)
            itensor::error("checkSpectrum: could not find eigenvector index");

        int degeneracy = 0;

        for(int n = 1; n <= dim(d); ++n){
            double lam = real(eltC(D, d(n), dp(n)));

            bool isGS = std::abs(lam - Emin) < tol;
            if(isGS) degeneracy++;

            std::cout << "Eigenvalue " << n << " = " << lam;
            if(isGS) std::cout << "   <-- GS";
            std::cout << "\n";

            // Für die ersten paar oder für GS-Zustände drucken
            if( isGS){
                ITensor psi = V * setElt(u(n));
                auto nrm = norm(psi);
                if(nrm > 0.0) psi /= nrm;

                std::cout << "State " << n << ":\n";
                std::cout << psi << "\n";
            }

            std::cout << "------------------------------------------\n";
        }

        std::cout << "Ground-state energy: " << Emin << "\n";
        std::cout << "Degeneracy: " << degeneracy << "\n";
        std::cout << "==========================================\n";
}








// ------------------------------------------ Algorithm Section ---------------------------------------------------
void ChainModel::basicModelLoop(const AlgoLoopParams& params) const {
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
    //cout << "norm psi\t\tk\tEnergy(<X>)\t\tInfidelity(|->)\n"; //"\t\tunitarity defect\n";
    double Ek, Fk, IFk;

    // Start U0 = Warm up
    ITensor U = U0_;
    ITensor psi0 = applyGate(U, p0_);

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H_, theta);
    ITensor Ad = adjointGateN(A);
    ITensor R0 = phaseOnStateN(p0_, Iset_, theta);

    for(int k = 0; k <= K; k++) {
        // Build ITensor gate from U and apply to phi(k-1)
        ITensor psi = applyGate(U, p0_);
        //cout << "Psi = " << psi << "\n";
        //cout << norm(psi) << "\t";

        // Calculate Observables and save them
        Ek = expectation(psi, H_);
        Fk = fidelityToSubspace(psi, groundspace_);
        IFk = 1 - Fk;
        rows.push_back(Row{k, Ek, Fk});

        //cout << k << "\t" << Ek << "\t" << IFk << "\n";
        if (IFk < infid_target) {
            //cout << "\nInfidelity of " << str_infid_target << " after " << k << " steps achieved!" << endl;
            break;
        }
        //cout << "--------------------------------- I Am Here ----------------------------------" << endl;
        // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
        ITensor Ud = adjointGateN(U);
        ITensor Ui = U;

        Ui = composeGateN(A, U, Iset_);
        Ui = composeGateN(Ui, R0, Iset_);
        Ui = composeGateN(Ui, Ud, Iset_);
        Ui = composeGateN(Ui, Ad, Iset_);
        //cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
        Ui = composeGateN(Ui, U, Iset_);

        //Next step
        if (unitary_DefectN(Ui, Iset_) > 1E-12) {
            Ui = reunitarize_polar_gateN(Ui, Iset_);
            U = Ui;
        }else{
            U = Ui;
        }
    }
    //write_csv("data4.csv", rows);
    write_csv("data" + s_string + ".csv", rows);

}

void ChainModel::degradingInfidLoop(const AlgoInfidParams& params) const {
    const double s_step = params.s_step;
    const string s_string = params.s_string;
    double theta = sqrt(s_step);
    const int K_max = params.K_max;
    const double infid_min = params.infid_min;
    const int infid_N = params.infid_N;
    cout << "Step s = " << s_step << "\n";

    // Data Container
    vector<RowInfid> rows;
    rows.reserve(infid_N+1);
    double Fk, IFk;
    //cout << "current infid\tsteps k" << endl;

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H_, theta);
    ITensor Ad = adjointGateN(A);
    ITensor R0 = phaseOnStateN(p0_, Iset_, theta);

    for (int ek = 0; ek < infid_N; ek++) {
        double infid_step;
        if (ek % 2 == 0) {
            infid_step = infid_min * pow(10.0, -ek);
        } else {
            infid_step = infid_min * pow(10.0, -ek);
        }
        //cout << infid_step << "\n";

        // Start U0 = Warm up
        ITensor U = U0_;

        for(int k = 0; k <= K_max; k++) {
            // Build ITensor gate from U and apply to phi(k-1)
            ITensor psi = applyGate(U, p0_);

            // Calculate Observables and save them
            Fk = fidelityToSubspace(psi, groundspace_);
            IFk = 1 - Fk;

            if (IFk > infid_step) {
                if (k<K_max) {
                    // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
                    ITensor Ud = adjointGateN(U);
                    ITensor Ui = U;

                    Ui = composeGateN(A, U, Iset_);
                    Ui = composeGateN(Ui, R0, Iset_);
                    Ui = composeGateN(Ui, Ud, Iset_);
                    Ui = composeGateN(Ui, Ad, Iset_);
                    //cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
                    Ui = composeGateN(Ui, U, Iset_);

                    //Next step
                    if (unitary_DefectN(Ui, Iset_) > 1E-12) {
                        //cout << "U: " << unitary_Defect(Ui, s_comb);
                        //cout << " -- start Unitarization -- ";
                        Ui = reunitarize_polar_gateN(Ui, Iset_);
                        //cout << "Unitarization complete -- ";
                        //cout << "U: " << unitary_Defect(Ui, s_comb) << "\n";
                        U = Ui;
                    }else{
                        U = Ui;
                    }

                }else if (k == K_max) {
                    //cout << k << endl;
                    rows.push_back(RowInfid{infid_step, k, IFk});
                }
            }else{
                //cout << k << endl;
                rows.push_back(RowInfid{infid_step, k, IFk});
                break;
            }
        }
        write_csv_eps("data" + s_string + "eps.csv", rows);
    }
}


void ChainModel::optimizeStepsLoop(const AlgoOptParams& params) const {
    const double s_min = params.s_min;
    const double s_bin = params.s_bin;
    const int s_N = params.s_N;
    const int K_max = params.K_max;
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;
    cout << "Infid Target = " << infid_target << "\n";

    // Data Container
    vector<RowOpt> rows;
    rows.reserve(s_bin+1);
    double Fk, IFk;
    //cout << "current s\tsteps k" << endl;

    for (int j = 0; j <= s_N; j++) {
        double s_step = s_min + s_bin*j;
        double theta = sqrt(s_step);
        //cout << s_step << "\t";

        // Start U0 = Warm up
        ITensor U = U0_;

        // Contruct Unitaries A and R
        ITensor A  = exp_i_theta_H(H_, theta);
        ITensor Ad = adjointGateN(A);
        ITensor R0 = phaseOnStateN(p0_, Iset_, theta);

        for(int k = 0; k <= K_max; k++) {
            // Build ITensor gate from U and apply to phi(k-1)
            ITensor psi = applyGate(U, p0_);

            // Calculate Observables and save them
            Fk = fidelityToSubspace(psi, groundspace_);
            IFk = 1 - Fk;

            if (IFk > infid_target) {
                if (k<K_max) {
                    // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
                    ITensor Ud = adjointGateN(U);
                    ITensor Ui = U;

                    Ui = composeGateN(A, U, Iset_);
                    Ui = composeGateN(Ui, R0, Iset_);
                    Ui = composeGateN(Ui, Ud, Iset_);
                    Ui = composeGateN(Ui, Ad, Iset_);
                    //cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
                    Ui = composeGateN(Ui, U, Iset_);

                    //Next step
                    if (unitary_DefectN(Ui, Iset_) > 1E-10) {
                        Ui = reunitarize_polar_gateN(Ui, Iset_);
                        U = Ui;
                    }else{
                        U = Ui;
                    }

                }else if (k == K_max) {
                    //cout << k << endl;
                    rows.push_back(RowOpt{s_step, k, IFk});
                }
            }else{
                //cout << k << endl;
                rows.push_back(RowOpt{s_step, k, IFk});
                break;
            }
        }
    }
    write_csv_opt("data" + str_infid_target + "opt.csv", rows);
}