//
// Created by Pascal Knoll on 25.03.26.
//
#include "BaseModel.h"
void BaseModel::checkHamiltonian() const{
    if(!(hasInds(H0_, prime(s_)) && hasInds(H0_, s_) && order(H0_) == 2)){
        cout << inds(H0_) << std::endl;
        itensor::error("H does not have indices (prime(s), s)");
    }
}

void BaseModel::printSummary() const{
    cout << "Site index: " << s_ << "\n";
    cout << "H = " << H0_ << "\n";
    cout << "phi = " << phi_ << "\n";
    cout << "ket0 = " << ket0_ << "\n";
}
void BaseModel::basicModelLoop(const AlgoLoopParams& params) const {
    // Parameters
    const double s_step = params.s_step;
    const string s_string = params.s_string;
    const int K = params.K;
    const double theta = sqrt(s_step);
    const double infid_target = params.infid_target;
    const string str_infid_target = params.str_infid_target;

    //Data container
    vector<Row> rows;
    rows.reserve(K+1);
    cout << "k\tEnergy(<X>)\t\tInfidelity(|->)\t\tunitarity defect\n";
    double Ek, Fk, IFk;

    // Start U0 = I
    ITensor U = idGate(s_);
    ITensor I = idGate(s_);

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H_, theta);
    ITensor Ad = adjointGate(A);
    ITensor R0 = phaseOnState(ket0_, s_, theta);

    for(int k = 0; k <= K; k++){
        // Build ITensor gate from U and apply to |0>
        ITensor psi = applyGate(U, ket0_);
        //cout << norm(psi) << "\t";

        // Calculate Observables and save them
        Ek = expectation(psi, H_);
        Fk = fidelity(psi, phi_);
        IFk = 1 - Fk;
        rows.push_back(Row{k, Ek, Fk});

        cout << k << "\t" << Ek << "\t" << IFk << "\t";
        if (IFk < infid_target) {
            cout << "\nInfidelity of " << str_infid_target << " after " << k << " steps achieved!" << endl;
            break;
        }

        // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
        ITensor Ud = adjointGate(U);
        ITensor Ui = U;

        Ui = composeGate(A, U, s_);
        Ui = composeGate(Ui, R0, s_);
        Ui = composeGate(Ui, Ud,s_);
        Ui = composeGate(Ui, Ad,s_);
        //cout << "-U: " << unitary_Defect(Ui,s_) << "\t";
        Ui = composeGate(Ui, U,s_);

        //Next step
        if (unitary_Defect(Ui, s_) > 1E-12) {
            cout << endl;
            //cout << "U: " << unitary_Defect(Ui,s_);
            //cout << " -- start Unitarization -- ";
            Ui = reunitarize_polar_gate(Ui, s_);
            //cout << "Unitarization complete -- ";
            //cout << "U: " << unitary_Defect(Ui,s_) << "\n";
            U = Ui;
        }else{
            U = Ui;
            cout << endl;
        }
    }
    //write_csv("data.csv", rows);
    write_csv("data" + s_string + ".csv", rows);
}


void BaseModel::degradingInfidLoop(const AlgoInfidParams& params) const {
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
    ITensor Ad = adjointGate(A);
    ITensor R0 = phaseOnState(ket0_, s_, theta);

    for (int ek = 0; ek < infid_N; ek++) {
        double infid_step;
        if (ek % 2 == 0) {
            infid_step = infid_min * pow(10.0, -ek);
        } else {
            infid_step = infid_min * pow(10.0, -ek);
        }
        cout << infid_step << "\t";

        // Start U0 = I
        ITensor U = idGate(s_);
        ITensor I = idGate(s_);

        for(int k = 0; k <= K_max; k++) {
            ITensor psi = applyGate(U, ket0_);

            // Calculate Observables and save them
            Fk = fidelity(psi, phi_);
            IFk = 1 - Fk;

            if (IFk > infid_step) {
                if (k<K_max) {
                    // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
                    ITensor Ud = adjointGate(U);
                    ITensor Ui = U;

                    Ui = composeGate(A, U, s_);
                    Ui = composeGate(Ui, R0, s_);
                    Ui = composeGate(Ui, Ud,s_);
                    Ui = composeGate(Ui, Ad,s_);
                    //cout << "-U: " << unitary_Defect(Ui,s_) << "\t";
                    Ui = composeGate(Ui, U,s_);

                    //Next step
                    if (unitary_Defect(Ui, s_) > 1E-12) {
                        //cout << endl;
                        //cout << "U: " << unitary_Defect(Ui,s_);
                        //cout << " -- start Unitarization -- ";
                        Ui = reunitarize_polar_gate(Ui, s_);
                        //cout << k << "Unitarization complete -- ";
                        //cout << "U: " << unitary_Defect(Ui,s_) << "\n";
                        U = Ui;
                    }else{
                        U = Ui;
                        //cout << endl;
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
    }
    write_csv_eps("data" + s_string + "eps.csv", rows);
}


void BaseModel::optimizeStepsLoop(const AlgoOptParams& params) const {
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
        ITensor U = idGate(s_);
        ITensor I = idGate(s_);

        // Contruct Unitaries A and R
        ITensor A  = exp_i_theta_H(H_, theta);
        ITensor Ad = adjointGate(A);
        ITensor R0 = phaseOnState(ket0_, s_, theta);

        for(int k = 0; k <= K_max; k++) {
            ITensor psi = applyGate(U, ket0_);

            // Calculate Observables and save them
            Fk = fidelity(psi, phi_);
            IFk = 1 - Fk;

            if (IFk > infid_target) {
                if (k<K_max) {
                    // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
                    ITensor Ud = adjointGate(U);
                    ITensor Ui = U;

                    Ui = composeGate(A, U, s_);
                    Ui = composeGate(Ui, R0, s_);
                    Ui = composeGate(Ui, Ud,s_);
                    Ui = composeGate(Ui, Ad,s_);
                    //cout << "-U: " << unitary_Defect(Ui,s_) << "\t";
                    Ui = composeGate(Ui, U,s_);

                    //Next step
                    if (unitary_Defect(Ui, s_) > 1E-12) {
                        //cout << endl;
                        //cout << "U: " << unitary_Defect(Ui,s_);
                        //cout << " -- start Unitarization -- ";
                        Ui = reunitarize_polar_gate(Ui, s_);
                        //cout << "Unitarization complete -- ";
                        //cout << "U: " << unitary_Defect(Ui,s_) << "\n";
                        U = Ui;
                    }else{
                        U = Ui;
                        //cout << endl;
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