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

    // Model
    TestModel model;
    TwoQubitTestModel model2(OneQubitOp::I, OneQubitOp::Z);

    if (true) {
        model2.printSummary();
        model2.printAsMatrix(model2.H());
        model2.printAsVector(model2.phi());
    }





    // ------------ Generate new Data or Plot existing Data -------------
    bool Data = false;
    bool Plot = false;

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
        //plot_with_python_S("data.csv", "plot.png"); // creates plot_energy.png + plot_fidelity.png
        //plot_with_python({"data0.2.csv", "data0.4.csv", "data0.5.csv", "data0.6.csv", "data0.8.csv", "data1.0.csv"}, "plot.png", "python3", "plot.py");

        //plot_with_python({"data1E-3opt.csv", "data1E-5opt.csv", "data1E-7opt.csv", "data1E-9opt.csv", "data1E-11opt.csv", "data1E-13opt.csv"}, "plotOpti.png", "python3", "plot.py", "opti");

        plot_with_python({"data0.6eps.csv","data0.5eps.csv", "data0.4eps.csv"}, "plot_eps.png", "python3", "plot.py", "infid_trace");
    }


    return 0;
}

