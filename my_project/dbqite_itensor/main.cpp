//
// Created by Pascal Knoll on 21.02.26.
//
#include "functions.h"
#include "data.h"
#include "TestModel.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

///////////////////////// Main ////////////////////////////
int main(){
    cout << fixed << setprecision(14);

    // Model
    TestModel model;

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


    //Main Loop
    bool Data = false;
    bool Plot = true;
    if (Data == true){
        //model.basicModelLoop(loop_params);
        model.optimizeStepsLoop(opt_params);
    }

    if (Plot == true) {
        //plot_with_python("data.csv", "plot.png"); // creates plot_energy.png + plot_fidelity.png
        //plot_with_python({"data0.2.csv", "data0.4.csv", "data0.5.csv", "data0.6.csv", "data0.8.csv", "data1.0.csv"}, "plot.png", "python3", "plot.py");
        //plot_with_python({"data0.01IF.csv","data0.05IF.csv","data0.1IF.csv","data0.2IF.csv", "data0.4IF.csv"}, "plotIF.png", "python3", "plot.py");
        //plot_with_python({"dataLt1.0.csv","dataLt2.0.csv","dataLt3.0.csv", "dataLt4.0.csv", "dataLt6.0.csv", "dataLt9.0.csv"}, "plotLt.png", "python3", "plot.py");


        plot_with_python({"data1E-3opt.csv", "data1E-5opt.csv", "data1E-7opt.csv", "data1E-9opt.csv", "data1E-11opt.csv", "data1E-13opt.csv"}, "plotOpti.png", "python3", "plot.py", "opti");
    }


    return 0;
}

