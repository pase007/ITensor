//
// Created by Pascal Knoll on 15.03.26.
//
#include "TwoQubitTestModel.h"

TwoQubitTestModel::TwoQubitTestModel(OneQubitOp op1, OneQubitOp op2)
    : s1_(2,"Qubit1"), s2_(2,"Qubit2"){
    ITensor H1 = BuildingHamiltonians::make(s1_, op1);
    ITensor H2 = BuildingHamiltonians::make(s2_, op2);
    ITensor I1 = BuildingHamiltonians::make(s1_, op2);
    ITensor I2 = BuildingHamiltonians::make(s2_, op1);

    H_ = ((H1 * H2) + (I1 * I2))/2;
    //cout << "H = " << H_ << endl;

    //H_ = transformToMatrix(H_);
    //cout << "transformed H" << endl;

    // IndexSet
    s12_ = IndexSet(s1_, s2_);


    // taget: ground state |-->
    ket01_ = ITensor(s1_);
    ket01_.set(s1_(1), 1.0/sqrt(2.0));
    ket01_.set(s1_(2), -1.0/sqrt(2.0));
    ket02_ = ITensor(s2_);
    ket02_.set(s2_(1), 1.0/sqrt(2.0));
    ket02_.set(s2_(2), -1.0/sqrt(2.0));

    ket0_ = ket01_ * ket02_;
    //phi_ = transformToVector(phi_);
    //cout << "transformed phi" << endl;

    // start state |psi0> = |00>
    psi1_ = ITensor(s1_);
    psi1_.set(s1_(1), 1.0);
    psi1_.set(s1_(2), 0.0);
    psi2_ = ITensor(s2_);
    psi2_.set(s2_(1), 1.0);
    psi2_.set(s2_(2), 0.0);

    psi0_ = psi1_ * psi2_;
}


void TwoQubitTestModel::printSummary() const{
    cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    cout << "H = " << H_ << "\n";
    //cout << "phi = " << phi_ << "\n";
    //cout << "ket0 = " << ket0_ << "\n";
}

void TwoQubitTestModel::print_state() const {
    //cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    //cout << "state Phi = " << phi_ << "\n";
}

ITensor TwoQubitTestModel::transformToMatrix(ITensor const& A) const {
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

ITensor TwoQubitTestModel::transformToVector(ITensor const& psi) const {
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


void TwoQubitTestModel::basicModelLoop(const AlgoLoopParams& params) const {
    const double s_step = params.s_step;
    double theta = sqrt(s_step);
    const string s_string = params.s_string;
    const int K = params.K;
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;

    //Data container
    vector<Row> rows;
    rows.reserve(K+1);
    cout << "norm psi\t\tk\tEnergy(<X>)\t\tInfidelity(|->)\t\tunitarity defect\n";
    double Ek, Fk, IFk;

    // Start U0 = I
    ITensor U = idGateN(s12_);
    ITensor I = idGateN(s12_);

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H_, theta);
    ITensor Ad = adjointGateN(A);
    ITensor R0 = phaseOnStateN(psi0_, s12_, theta);

    for(int k = 0; k <= K; k++) {
        // Build ITensor gate from U and apply to phi(k-1)
        ITensor psi = applyGate(U, psi0_);
        cout << norm(psi) << "\t";

        // Calculate Observables and save them
        Ek = expectation(psi, H_);
        Fk = fidelity(psi, ket0_);
        IFk = 1 - Fk;
        rows.push_back(Row{k, Ek, Fk});

        cout << k << "\t" << Ek << "\t" << IFk << "\t";
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
        cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
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
    write_csv("data" + s_string + "D.csv", rows);

}

void TwoQubitTestModel::degradingInfidLoop(const AlgoInfidParams& params) const {
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
    ITensor R0 = phaseOnStateN(psi0_, s12_, theta);

    for (int ek = 0; ek < infid_N; ek++) {
        double infid_step;
        if (ek % 2 == 0) {
            infid_step = infid_min * pow(10.0, -ek);
        } else {
            infid_step = infid_min * pow(10.0, -ek);
        }
        cout << infid_step << "\t";

        // Start U0 = I
        ITensor U = idGateN(s12_);
        ITensor I = idGateN(s12_);

        for(int k = 0; k <= K_max; k++) {
            // Build ITensor gate from U and apply to phi(k-1)
            ITensor psi = applyGate(U, psi0_);

            // Calculate Observables and save them
            Fk = fidelity(psi, ket0_);
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
                    cout << "-U: " << unitary_DefectN(Ui, s12_) << "\n";
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
        write_csv_eps("data" + s_string + "epsD2.csv", rows);
    }
}


void TwoQubitTestModel::optimizeStepsLoop(const AlgoOptParams& params) const {
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

        // Start U0 = I
        ITensor U = idGateN(s12_);
        ITensor I = idGateN(s12_);

        // Contruct Unitaries A and R
        ITensor A  = exp_i_theta_H(H_, theta);
        ITensor Ad = adjointGateN(A);
        ITensor R0 = phaseOnStateN(psi0_, s12_, theta);

        for(int k = 0; k <= K_max; k++) {
            // Build ITensor gate from U and apply to phi(k-1)
            ITensor psi = applyGate(U, psi0_);

            // Calculate Observables and save them
            Fk = fidelity(psi, ket0_);
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
    write_csv_opt("data" + str_infid_target + "optD2.csv", rows);
}





















