//
// Created by Pascal Knoll on 21.02.26.
//
#include "functions.h"
#include "data.h"
#include "TestModel.h"
#include "StructFile.h"
#include "TwoQubitTestModel.h"
#include "BuildingHamiltonians.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

///////////////////////// Main ////////////////////////////
int main(){
    cout << fixed << setprecision(14);
    TestModel model;
    TwoQubitTestModel model3(OneQubitOp::I, OneQubitOp::X);
    //model3.printSummary();

    // Model
    double theta = sqrt(0.5);
    double infid_target = 1E-10;
    string str_infid_target = "1E-10";
    int K = 30;
    TwoQubitTestModel model2(OneQubitOp::I, OneQubitOp::X);
    //model2.printSummary();
    ITensor ket0 = model2.ket0();
    ITensor psi0 = model2.psi();
    ITensor H = model2.H();
    cout << "H = " << H << endl;
    IndexSet s12 = IndexSet(model2.s1(), model2.s2());

    //Data container
    vector<Row> rows;
    rows.reserve(K+1);
    cout << "norm psi\tk\tEnergy(<X>)\t\tInfidelity(|->)\t\tunitarity defect\n";
    double Ek, Fk, IFk;

    // Start U0 = I
    ITensor U = idGateN(s12);
    ITensor I = idGateN(s12);
    //cout << "U0 = " << U << endl;
    //cout << "ket0 = " << ket0 << endl;

    // Contruct Unitaries A and R
    ITensor A  = exp_i_theta_H(H, theta);
    //cout << "A = " << A << endl;
    ITensor Ad = adjointGateN(A);
    //cout << "Ad = " << Ad << endl;
    ITensor R0 = phaseOnStateN(psi0, s12, theta);
    //cout << "R0 = " << R0 << endl;

    for(int k = 0; k <= K; k++){
        // Build ITensor gate from U and apply to phi(k-1)
        ITensor psi = applyGate(U, psi0);
        cout << norm(psi) << "\t";

        // Calculate Observables and save them
        Ek = expectation(psi, H);
        Fk = fidelity(psi, ket0);
        IFk = 1 - Fk;
        rows.push_back(Row{k, Ek, Fk});

        cout << k << "\t" << Ek << "\t" << IFk << "\t";
        if (IFk < infid_target) {
            cout << "\nInfidelity of " << str_infid_target << " after " << k << " steps achieved!" << endl;
            break;
        }

        // DB-QITE recursion: U_{k+1} = A * U * Rk * U' * A' * U
        ITensor Ud = adjointGateN(U);
        //cout << "U = " << U << endl;
        //cout << "Ud = " << Ud << endl;
        //cout << "--------------------------------- I Am Here ----------------------------------" << endl;
        ITensor Ui = U;

        Ui = composeGateN(A, U, s12);
        Ui = composeGateN(Ui, R0, s12);
        Ui = composeGateN(Ui, Ud, s12);
        Ui = composeGateN(Ui, Ad, s12);
        cout << "-U: " << unitary_DefectN(Ui, s12) << "\n";
        Ui = composeGateN(Ui, U, s12);
        //cout << "Ui = " << Ui << "\n";

        //Next step
        if (unitary_DefectN(Ui, s12) > 1E-12) {
            //cout << "U: " << unitary_Defect(Ui, s_comb);
            //cout << " -- start Unitarization -- ";
            //Ui = reunitarize_polar_gate(Ui, s_comb);
            //cout << "Unitarization complete -- ";
            //cout << "U: " << unitary_Defect(Ui, s_comb) << "\n";
            U = Ui;
        }else{
            U = Ui;
        }
    }
    write_csv("data4.csv", rows);








    // ------------ Generate new Data or Plot existing Data -------------
    bool Data = false;
    bool Plot = true;

    if (Data == true){
        // Params for Fidelity in time Algo
        AlgoLoopParams loop_params;
        loop_params.s_step = 0.6;
        loop_params.s_string = "0.6";
        loop_params.K = 100;
        loop_params.infid_target = 1E-8;
        loop_params.str_infid_target = "1E-8";

        // Params for optimization for given epsilon
        AlgoOptParams opt_params;
        opt_params.s_min = 0.05;
        opt_params.s_bin = 0.05;
        opt_params.s_N = 20;
        opt_params.K_max = 100;
        opt_params.infid_target = 1E-13;
        opt_params.str_infid_target = "1E-13";

        // Params for reacing target infidelities
        AlgoInfidParams eps_params;
        eps_params.s_step = 0.6;
        eps_params.s_string = "0.6";
        eps_params.K_max = 150;
        eps_params.infid_min = 1E-2;
        eps_params.infid_N = 15;

        //Run algorithm
        //model.basicModelLoop(loop_params);
        //model.optimizeStepsLoop(opt_params);
        model.degradingInfidLoop(eps_params);
    }

    if (Plot == true) {
        plot_with_python_S("data4.csv", "plot.png"); // creates plot_energy.png + plot_fidelity.png
        //plot_with_python({"data0.2.csv", "data0.4.csv", "data0.5.csv", "data0.6.csv", "data0.8.csv", "data1.0.csv"}, "plot.png", "python3", "plot.py");

        //plot_with_python({"data1E-3opt.csv", "data1E-5opt.csv", "data1E-7opt.csv", "data1E-9opt.csv", "data1E-11opt.csv", "data1E-13opt.csv"}, "plotOpti.png", "python3", "plot.py", "opti");

        //plot_with_python({"data0.6eps.csv","data0.5eps.csv", "data0.4eps.csv"}, "plot_eps.png", "python3", "plot.py", "infid_trace");
    }


    return 0;
}

