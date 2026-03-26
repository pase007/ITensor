//
// Created by Pascal Knoll on 21.02.26.
//
#include "functions.h"
#include "data.h"
#include "TestModel.h"
#include "StructFile.h"
#include "TwoQubitTestModel.h"
#include "BuildingHamiltonians.h"
#include "OneSigmaSite.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void getData1(TestModel& model){
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
    opt_params.infid_target = 1E-3;
    opt_params.str_infid_target = "1E-3";

    // Params for reacing target infidelities
    AlgoInfidParams eps_params;
    eps_params.s_step = 0.2;
    eps_params.s_string = "0.2";
    eps_params.K_max = 150;
    eps_params.infid_min = 1E-2;
    eps_params.infid_N = 15;

    //Run algorithm
    //model.basicModelLoop(loop_params);
    //model.optimizeStepsLoop(opt_params);
    model.degradingInfidLoop(eps_params);
}

void getData2(TwoQubitTestModel& model2){
	// Params for Fidelity in time Algo
    vector<double> s_step_vec = {0.1, 0.6, 0.5, 0.4, 0.2};
    vector<string> s_step_vec_str = {"0.1", "0.6", "0.5", "0.4", "0.2"};

    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoLoopParams loop_params2;
        loop_params2.s_step = s_step_vec[i];
        loop_params2.s_string = s_step_vec_str[i];
        loop_params2.K = 100;
        loop_params2.infid_target = 1E-14;
        loop_params2.str_infid_target = "1E-14";
        //model2.basicModelLoop(loop_params2);
    }


    // Params for optimization for given epsilon
    vector<double> e_step_vec = {1E-3, 1E-5, 1E-7, 1E-9, 1E-11, 1E-13};
    vector<string> e_step_vec_str = {"1E-3", "1E-5", "1E-7", "1E-9", "1E-11", "1E-13"};

    for (int i = 0; i < e_step_vec.size(); i++) {
        AlgoOptParams opt_params2;
        opt_params2.s_min = 0.3;
        opt_params2.s_bin = 0.05;
        opt_params2.s_N = 12;
        opt_params2.K_max = 100;
        opt_params2.infid_target = e_step_vec[i];
        opt_params2.str_infid_target = e_step_vec_str[i];
        //model2.optimizeStepsLoop(opt_params2);
    }

    // Params for reacing target infidelities
    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoInfidParams eps_params2;
        eps_params2.s_step = s_step_vec[i];
        eps_params2.s_string = s_step_vec_str[i];
        eps_params2.K_max = 100;
        eps_params2.infid_min = 1E-2;
        eps_params2.infid_N = 14;
        //model2.degradingInfidLoop(eps_params2);
    }
}

void getDataSigma1(OneSigmaSite& sigma1){
	vector<double> s_step_vec = {0.1, 0.6, 0.5, 0.4, 0.2, 0.8, 1.0};
    vector<string> s_step_vec_str = {"0.1S", "0.6S", "0.5S", "0.4S", "0.2S", "0.8S", "1.0S"};

	// Params for Fidelity in time Algo
	for (int i = 0; i < s_step_vec.size(); i++) {
	    AlgoLoopParams loop_params;
    	loop_params.s_step = s_step_vec[i];
    	loop_params.s_string = s_step_vec_str[i];
    	loop_params.K = 100;
    	loop_params.infid_target = 1E-8;
    	loop_params.str_infid_target = "1E-8S";
		sigma1.basicModelLoop(loop_params);
	}

    // Params for optimization for given epsilon
	vector<double> e_step_vec = {1E-3, 1E-5, 1E-7, 1E-9, 1E-11, 1E-13};
    vector<string> e_step_vec_str = {"1E-3S", "1E-5S", "1E-7S", "1E-9S", "1E-11S", "1E-13S"};

	for (int i = 0; i < e_step_vec.size(); i++) {
    	AlgoOptParams opt_params;
	    opt_params.s_min = 0.05;
    	opt_params.s_bin = 0.05;
    	opt_params.s_N = 20;
	    opt_params.K_max = 100;
    	opt_params.infid_target = e_step_vec[i];
    	opt_params.str_infid_target = e_step_vec_str[i];
		sigma1.optimizeStepsLoop(opt_params);
	}

    // Params for reacing target infidelities
    for (int i = 0; i < s_step_vec.size(); i++) {
    	AlgoInfidParams eps_params;
    	eps_params.s_step = s_step_vec[i];
    	eps_params.s_string = s_step_vec_str[i];
    	eps_params.K_max = 150;
    	eps_params.infid_min = 1E-2;
    	eps_params.infid_N = 15;
		sigma1.degradingInfidLoop(eps_params);
	}
}

///////////////////////// Main ////////////////////////////
int main(){
    cout << fixed << setprecision(14);

    // ------------ Generate new Data or Plot existing Data -------------
    bool Data1 = false;
    bool Data2 = false;
	bool DataSig1 = true;
    bool Plot = true;

	//Test of Principle with 2 Hamiltonians
	if (Data1 == true){
		TestModel model;
        getData1(model);
    }
    if (Data2 == true) {
		TwoQubitTestModel model2(OneQubitOp::I, OneQubitOp::X);
        getData2(model2);
    }

    //Sigma Model Algos
	if (DataSig1 == true){
		double g = 1.0;
		OneSigmaSite sigma1(g);
		sigma1.printSummary();
		getDataSigma1(sigma1);
	}

    if (Plot == true) {
        //plot_with_python_S("data0.5D.csv", "plot.png", "python3", "plot.py");

        //plot_with_python({"data0.2.csv", "data0.4.csv", "data0.5.csv", "data0.6.csv", "data0.8.csv", "data1.0.csv"}, "plot.png", "python3", "plot.py");
        //plot_with_python({ "data0.1D.csv", "data0.6D.csv", "data0.5D.csv", "data0.2D.csv"}, "plotN2.png", "python3", "plot.py"); // creates plot_energy.png + plot_fidelity.png
        plot_with_python({"data0.2S.csv", "data0.4S.csv", "data0.5S.csv", "data0.6S.csv", "data0.8S.csv", "data1.0S.csv"}, "plotSig.png", "python3", "plot.py");


        //plot_with_python({"data1E-3opt.csv", "data1E-5opt.csv", "data1E-7opt.csv", "data1E-9opt.csv", "data1E-11opt.csv", "data1E-13opt.csv"}, "plotOpti.png", "python3", "plot.py", "opti");
        //plot_with_python({"data1E-3optD2.csv", "data1E-5optD2.csv", "data1E-7optD2.csv", "data1E-9optD2.csv", "data1E-11optD2.csv", "data1E-13optD2.csv"}, "plotOptiD2.png", "python3", "plot.py", "opti");
        plot_with_python({"data1E-3Sopt.csv", "data1E-5Sopt.csv", "data1E-7Sopt.csv", "data1E-9Sopt.csv", "data1E-11Sopt.csv", "data1E-13Sopt.csv"}, "plotOptiSig.png", "python3", "plot.py", "opti");


        //plot_with_python({"data0.8epsD2.csv", "data0.6epsD2.csv","data0.5epsD2.csv", "data0.4epsD2.csv"}, "plot_epsD2_fit.png", "python3", "plot.py", "infid_fit_trace");
        //plot_with_python({"data0.8eps.csv", "data0.6eps.csv","data0.5eps.csv", "data0.4eps.csv", "data0.2eps.csv"}, "plot_eps_fit.png", "python3", "plot.py", "infid_fit_trace");
        plot_with_python({"data0.8Seps.csv", "data0.6Seps.csv","data0.5Seps.csv", "data0.4Seps.csv", "data0.2Seps.csv"}, "plot_eps_sig.png", "python3", "plot.py", "infid_trace");

    }


    return 0;
}



