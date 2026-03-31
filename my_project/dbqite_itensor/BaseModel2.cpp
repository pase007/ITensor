//
// Created by Pascal Knoll on 26.03.26.
//

#include "BaseModel2.h"

// Printing Details
void BaseModel2::printSummary() const{
    cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    //cout << "H = " << H_ << "\n";
    cout << "GS: E = " << E0_ << " and Xi = " << groundspace_[0] << "\n";
    cout << "start-state: psi0 = " << psi0_ << "\n";
}

void BaseModel2::print_state() const {
    //cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    //cout << "state Phi = " << phi_ << "\n";
}


pair<double, vector<ITensor>> BaseModel2::groundSpace2(double tol) const{
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
        if(I != s1_ && I != s2_){
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

void BaseModel2::checkGroundStateDegeneracy(double tol) const{
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
        Index u;
        bool found_u = false;
        for(auto const& I : inds(V)){
            if(I != s1_ && I != s2_){
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



// Transformation between Representations
ITensor BaseModel2::transformToMatrix(ITensor const& A) const {
    vector<Index> ket_inds;
    vector<Index> bra_inds;

    for(auto const& i : inds(A)){
        if(primeLevel(i) == 0)
            ket_inds.push_back(i);
    }

    reverse(ket_inds.begin(), ket_inds.end());

    for(auto const& k : ket_inds){
        bra_inds.push_back(prime(k));
    }

    auto [cket, ket_combined] = combiner(ket_inds);
    auto [cbra, bra_combined] = combiner(bra_inds);

    ITensor M = cbra * (cket * A);

    // only prime the first combined index
    auto bra_p = prime(bra_combined);
    M.replaceInds({bra_combined}, {bra_p});

    //cout << "bra prime level = " << primeLevel(bra_p) << "\n";
    //cout << "ket prime level = " << primeLevel(ket_combined) << "\n";

    //cout << inds(M) << endl;
    //cout << "M = " << M << "\n";

    return M;
}

ITensor BaseModel2::transformToVector(ITensor const& psi) const {
    vector<Index> ket_inds;

    for(auto const& i : inds(psi)) {
        if(primeLevel(i) == 0)
            ket_inds.push_back(i);
    }

    reverse(ket_inds.begin(), ket_inds.end());

    auto [cket, ket_combined] = combiner(ket_inds);
    auto V = cket * psi;

    //cout << "V = " << V << "\n";
    //cout << inds(V) << endl;
    //cout << "V prime level = " << primeLevel(inds(V)[0]) << "\n";
    return V;
}



// ------------------------------------------ Algorithm Section ---------------------------------------------------
void BaseModel2::basicModelLoop(const AlgoLoopParams& params) const {
    const double s_step = params.s_step;
    double theta = sqrt(s_step);
    const string s_string = params.s_string;
    const int K = params.K;
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;

    //Data container
    vector<Row> rows;
    rows.reserve(K+1);
    cout << "norm psi\t\tk\tEnergy(<X>)\t\tInfidelity(|->)\n"; //"\t\tunitarity defect\n";
    double Ek, Fk, IFk;

    // Start U0 = Warm up
    ITensor U = U0_;
    ITensor psi0 = applyGate(U, p0_);

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H_, theta);
    ITensor Ad = adjointGateN(A);
    ITensor R0 = phaseOnStateN(p0_, s12_, theta);

    for(int k = 0; k <= K; k++) {
        // Build ITensor gate from U and apply to phi(k-1)
        ITensor psi = applyGate(U, p0_);
        //cout << "Psi = " << psi << "\n";
        cout << norm(psi) << "\t";

        // Calculate Observables and save them
        Ek = expectation(psi, H_);
        Fk = fidelityToSubspace(psi, groundspace_);
        IFk = 1 - Fk;
        rows.push_back(Row{k, Ek, Fk});

        cout << k << "\t" << Ek << "\t" << IFk << "\n";
        if (IFk < infid_target) {
            cout << "\nInfidelity of " << str_infid_target << " after " << k << " steps achieved!" << endl;
            break;
        }
        //cout << "--------------------------------- I Am Here ----------------------------------" << endl;
        // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
        ITensor Ud = adjointGateN(U);
        ITensor Ui = U;

        Ui = composeGateN(A, U, s12_);
        Ui = composeGateN(Ui, R0, s12_);
        Ui = composeGateN(Ui, Ud, s12_);
        Ui = composeGateN(Ui, Ad, s12_);
        //cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
        Ui = composeGateN(Ui, U, s12_);

        //Next step
        if (unitary_DefectN(Ui, s12_) > 1E-12) {
            Ui = reunitarize_polar_gateN(Ui, s12_);
            U = Ui;
        }else{
            U = Ui;
        }
    }
    //write_csv("data4.csv", rows);
    write_csv("data" + s_string + ".csv", rows);

}

void BaseModel2::degradingInfidLoop(const AlgoInfidParams& params) const {
    const double s_step = params.s_step;
    const string s_string = params.s_string;
    double theta = sqrt(s_step);
    const int K_max = params.K_max;
    const double infid_min = params.infid_min;
    const int infid_N = params.infid_N;

    // Data Container
    vector<RowInfid> rows;
    rows.reserve(infid_N+1);
    double Fk, IFk;
    cout << "current infid\tsteps k" << endl;

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H_, theta);
    ITensor Ad = adjointGateN(A);
    ITensor R0 = phaseOnStateN(p0_, s12_, theta);

    for (int ek = 0; ek < infid_N; ek++) {
        double infid_step;
        if (ek % 2 == 0) {
            infid_step = infid_min * pow(10.0, -ek);
        } else {
            infid_step = infid_min * pow(10.0, -ek);
        }
        cout << infid_step << "\t";

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

                    Ui = composeGateN(A, U, s12_);
                    Ui = composeGateN(Ui, R0, s12_);
                    Ui = composeGateN(Ui, Ud, s12_);
                    Ui = composeGateN(Ui, Ad, s12_);
                    //cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
                    Ui = composeGateN(Ui, U, s12_);

                    //Next step
                    if (unitary_DefectN(Ui, s12_) > 1E-12) {
                        //cout << "U: " << unitary_Defect(Ui, s_comb);
                        //cout << " -- start Unitarization -- ";
                        Ui = reunitarize_polar_gateN(Ui, s12_);
                        //cout << "Unitarization complete -- ";
                        //cout << "U: " << unitary_Defect(Ui, s_comb) << "\n";
                        U = Ui;
                    }else{
                        U = Ui;
                    }

                }else if (k == K_max) {
                    cout << k << endl;
                    rows.push_back(RowInfid{infid_step, k, IFk});
                }
            }else{
                cout << k << endl;
                rows.push_back(RowInfid{infid_step, k, IFk});
                break;
            }
        }
        write_csv_eps("data" + s_string + "eps.csv", rows);
    }
}


void BaseModel2::optimizeStepsLoop(const AlgoOptParams& params) const {
    const double s_min = params.s_min;
    const double s_bin = params.s_bin;
    const int s_N = params.s_N;
    const int K_max = params.K_max;
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;

    // Data Container
    vector<RowOpt> rows;
    rows.reserve(s_bin+1);
    double Fk, IFk;
    cout << "current s\tsteps k" << endl;

    for (int j = 0; j <= s_N; j++) {
        double s_step = s_min + s_bin*j;
        double theta = sqrt(s_step);
        cout << s_step << "\t";

        // Start U0 = Warm up
        ITensor U = U0_;

        // Contruct Unitaries A and R
        ITensor A  = exp_i_theta_H(H_, theta);
        ITensor Ad = adjointGateN(A);
        ITensor R0 = phaseOnStateN(p0_, s12_, theta);

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

                    Ui = composeGateN(A, U, s12_);
                    Ui = composeGateN(Ui, R0, s12_);
                    Ui = composeGateN(Ui, Ud, s12_);
                    Ui = composeGateN(Ui, Ad, s12_);
                    //cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
                    Ui = composeGateN(Ui, U, s12_);

                    //Next step
                    if (unitary_DefectN(Ui, s12_) > 1E-10) {
                        Ui = reunitarize_polar_gateN(Ui, s12_);
                        U = Ui;
                    }else{
                        U = Ui;
                    }

                }else if (k == K_max) {
                    cout << k << endl;
                    rows.push_back(RowOpt{s_step, k, IFk});
                }
            }else{
                cout << k << endl;
                rows.push_back(RowOpt{s_step, k, IFk});
                break;
            }
        }
    }
    write_csv_opt("data" + str_infid_target + "opt.csv", rows);
}