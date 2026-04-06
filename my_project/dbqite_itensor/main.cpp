//
// Created by Pascal Knoll on 21.02.26.
//
#include "functions.h"
#include "data.h"
#include "StructFile.h"
#include "BuildingHamiltonians.h"
#include "TestModel.h"
#include "TwoQubitTestModel.h"
#include "OneSigmaSite.h"
#include "TwoSigmaSite.h"
#include "ChainModel.h"
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
    vector<double> s_step_vec = {0.1, 0.6, 0.5, 0.4, 0.2, 0.8};
    vector<string> s_step_vec_str = {"0.1D2", "0.6D2", "0.5D2", "0.4D2", "0.2D2", "0.8D2"};

    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoLoopParams loop_params2;
        loop_params2.s_step = s_step_vec[i];
        loop_params2.s_string = s_step_vec_str[i];
        loop_params2.K = 100;
        loop_params2.infid_target = 1E-14;
        loop_params2.str_infid_target = "1E-14";
        model2.basicModelLoop(loop_params2);
    }

    // Params for optimization for given epsilon
    vector<double> e_step_vec = {1E-3, 1E-5, 1E-7, 1E-9, 1E-11, 1E-13};
    vector<string> e_step_vec_str = {"1E-3D2", "1E-5D2", "1E-7D2", "1E-9D2", "1E-11D2", "1E-13D2"};

    for (int i = 0; i < e_step_vec.size(); i++) {
        AlgoOptParams opt_params2;
        opt_params2.s_min = 0.3;
        opt_params2.s_bin = 0.05;
        opt_params2.s_N = 12;
        opt_params2.K_max = 100;
        opt_params2.infid_target = e_step_vec[i];
        opt_params2.str_infid_target = e_step_vec_str[i];
        model2.optimizeStepsLoop(opt_params2);
    }

    // Params for reacing target infidelities
    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoInfidParams eps_params2;
        eps_params2.s_step = s_step_vec[i];
        eps_params2.s_string = s_step_vec_str[i];
        eps_params2.K_max = 100;
        eps_params2.infid_min = 1E-2;
        eps_params2.infid_N = 14;
        model2.degradingInfidLoop(eps_params2);
    }
}



void getSigma2(TwoSigmaSite& sigma2){
	// Params for Fidelity in time Algo
    vector<double> s_step_vec = {0.1, 0.6, 0.5, 0.4, 0.2, 0.8, 1.0};
    vector<string> s_step_vec_str = {"0.1SS", "0.6SS", "0.5SS", "0.4SS", "0.2SS", "0.8SS", "1.0SS"};

    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoLoopParams loop_params2;
        loop_params2.s_step = s_step_vec[i];
        loop_params2.s_string = s_step_vec_str[i];
        loop_params2.K = 100;
        loop_params2.infid_target = 1E-14;
        loop_params2.str_infid_target = "1E-14";
        sigma2.basicModelLoop(loop_params2);
    }

    // Params for optimization for given epsilon
    vector<double> e_step_vec = {1E-3, 1E-5, 1E-7, 1E-9, 1E-11, 1E-13};
    vector<string> e_step_vec_str = {"1E-3SS", "1E-5SS", "1E-7SS", "1E-9SS", "1E-11SS", "1E-13SS"};

    for (int i = 0; i < e_step_vec.size(); i++) {
        AlgoOptParams opt_params2;
        opt_params2.s_min = 0.05;
        opt_params2.s_bin = 0.05;
        opt_params2.s_N = 25;
        opt_params2.K_max = 100;
        opt_params2.infid_target = e_step_vec[i];
        opt_params2.str_infid_target = e_step_vec_str[i];
        sigma2.optimizeStepsLoop(opt_params2);
    }

    // Params for reacing target infidelities
    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoInfidParams eps_params2;
        eps_params2.s_step = s_step_vec[i];
        eps_params2.s_string = s_step_vec_str[i];
        eps_params2.K_max = 100;
        eps_params2.infid_min = 1E-2;
        eps_params2.infid_N = 8;
        sigma2.degradingInfidLoop(eps_params2);
    }
	//plot_with_python({"data0.2SS.csv", "data0.4SS.csv", "data0.5SS.csv", "data0.6SS.csv", "data0.8SS.csv", "data1.0SS.csv"}, "plotSig2.png", "python3", "plot.py");
    //plot_with_python({"data1E-3SSopt.csv", "data1E-5SSopt.csv", "data1E-7SSopt.csv", "data1E-9SSopt.csv", "data1E-11SSopt.csv", "data1E-13SSopt.csv"}, "plotOptiSig2.png", "python3", "plot.py", "opti");
    //plot_with_python({"data0.8SSeps.csv", "data0.6SSeps.csv", "data0.5SSeps.csv", "data0.2SSeps.csv"}, "plot_eps_sig2.png", "python3", "plot.py", "infid_fit_trace");

}

void getChainN(ChainModel& chain){
	// Params for Fidelity in time Algo
    vector<double> s_step_vec = { 0.2, 0.21, 0.22, 0.23, 0.242};//, 0.25, 0.255, 0.26, 0.265, 0.275};
    vector<string> s_step_vec_str = {"0.2Ch2", "0.21Ch2", "0.22Ch2", "0.23Ch2", "0.242Ch2"};//, "0.25Ch2", "0.255Ch2", "0.26Ch2", "0.265Ch2", "0.275Ch2"};

    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoLoopParams loop_params2;
        loop_params2.s_step = s_step_vec[i];
        loop_params2.s_string = s_step_vec_str[i];
        loop_params2.K = 300;
        loop_params2.infid_target = 1E-14;
        loop_params2.str_infid_target = "1E-14";
        //chain.basicModelLoop(loop_params2);
    }
    //plot_with_python_E({"data0.22Ch2.csv", "data0.23Ch2.csv", "data0.242Ch2.csv", "data0.25Ch2.csv", "data0.255Ch2.csv", "data0.26Ch2.csv", "data0.265Ch2.csv", "data0.275Ch2.csv"}, "Elevels4.csv", "plotCh42.png",  "python3", "plot.py", "multiE");

   // Params for optimization for given epsilon
    vector<double> e_step_vec = {1E-5, 1E-7, 1E-9, 1E-11};
    vector<string> e_step_vec_str = {"1E-5Ch2", "1E-7Ch2", "1E-9Ch2", "1E-11Ch2", "1E-13Ch2"};

    for (int i = 0; i < e_step_vec.size(); i++) {
        AlgoOptParams opt_params2;
        opt_params2.s_min = 0.06;
        opt_params2.s_bin = 0.005;
        opt_params2.s_N = 14;
        opt_params2.K_max = 900;
        opt_params2.infid_target = e_step_vec[i];
        opt_params2.str_infid_target = e_step_vec_str[i];
        chain.optimizeStepsLoop(opt_params2);
    }

    // Params for reacing target infidelities
    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoInfidParams eps_params2;
        eps_params2.s_step = s_step_vec[i];
        eps_params2.s_string = s_step_vec_str[i];
        eps_params2.K_max = 200;
        eps_params2.infid_min = 1E-2;
        eps_params2.infid_N = 9;
        //chain.degradingInfidLoop(eps_params2);
    }
    //plot_with_python({"data0.2Ch2.csv", "data0.25Ch2.csv", "data0.3Ch2.csv", "data0.125Ch2.csv", "data0.175Ch2.csv", "data0.275Ch2.csv"}, "plotCh42.png",  "python3", "plot.py", "multi");
    plot_with_python({ "data1E-5Ch2opt.csv", "data1E-7Ch2opt.csv", "data1E-9Ch2opt.csv", "data1E-11Ch2opt.csv"}, "plotOptiCh4255.png", "python3", "plot.py", "opti");
    //plot_with_python({  "data0.2Ch2eps.csv", "data0.21Ch2eps.csv", "data0.22Ch2eps.csv", "data0.23Ch2eps.csv"}, "plot_eps_Ch425.png", "python3", "plot.py", "infid_fit_trace");

}

void scanSpectrum(){
	string filename = "ExactEnergies3.txt";

	vector<vector<double>> FullEnergies;
	vector<double> g_vals;
	double g_step = 0.05;
	double g0 = 0.4;
	for (int i = 0; i < 18; i++){
		double g = g0 + i*g_step;
		g_vals.push_back(g);
		ChainModel chainN = ChainModel(g, 3, false);
		vector<double> ei = chainN.ExactEnergies();
		FullEnergies.push_back(ei);

	}
	write_csv_Espectrum(filename, g_vals, FullEnergies);
	plot_Espectrum("ExactEnergies3.txt", "OptimumSteps3.txt", "E3spectrum_gscan.png", "python3", "plot.py" );
}

///////////////////////// Main ////////////////////////////
int main(){
    cout << fixed << setprecision(14);
	//scanSpectrum();
	///plot_Espectrum("ExactEnergies3.txt", "OptimumSteps3.txt", "E3spectrum_gscan.png", "python3", "plot.py" );
	if(false){
		ChainModel chainN = ChainModel(1.0, 3, false);
		AlgoOptParams opt_params2 = {0.2, 0.05, 40, 300, 1E-11, "1E-11CGap3"};
		chainN.optimizeStepsLoop(opt_params2);
    	plot_with_python({ "data1E-11CGap3opt.csv"}, "plotGap3.png", "python3", "plot.py", "opti");
	}



    // ------------ Generate new Data or Plot existing Data -------------
    bool Data1 = false;
    bool Data2 = false;
	bool DataSig2 = false;
	bool DataChain = false;
    bool Plot = false;

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
	if (DataSig2 == true){
		TwoSigmaSite sigma2(1.0);
		sigma2.printSummary();
		getSigma2(sigma2);
	}
	if (DataChain == true){
		ChainModel chainN = ChainModel(0.45, 3, false);
		getChainN(chainN);
	}

    if (Plot == true) {
        //plot_with_python_S("data0.5D.csv", "plot.png", "python3", "plot.py");

        //plot_with_python({"data0.2.csv", "data0.4.csv", "data0.5.csv", "data0.6.csv", "data0.8.csv", "data1.0.csv"}, "plot.png", "python3", "plot.py");
        //plot_with_python({ "data0.1D2.csv", "data0.6D2.csv", "data0.5D2.csv", "data0.2D2.csv"}, "plotN2.png", "python3", "plot.py"); // creates plot_energy.png + plot_fidelity.png
        //plot_with_python({"data0.2S.csv", "data0.4S.csv", "data0.5S.csv", "data0.6S.csv", "data0.8S.csv", "data1.0S.csv"}, "plotSig.png", "python3", "plot.py");


        //plot_with_python({"data1E-3opt.csv", "data1E-5opt.csv", "data1E-7opt.csv", "data1E-9opt.csv", "data1E-11opt.csv", "data1E-13opt.csv"}, "plotOpti.png", "python3", "plot.py", "opti");
        //plot_with_python({"data1E-3optD2.csv", "data1E-5optD2.csv", "data1E-7optD2.csv", "data1E-9optD2.csv", "data1E-11optD2.csv", "data1E-13optD2.csv"}, "plotOptiD2.png", "python3", "plot.py", "opti");
        //plot_with_python({"data1E-3Sopt.csv", "data1E-5Sopt.csv", "data1E-7Sopt.csv", "data1E-9Sopt.csv", "data1E-11Sopt.csv", "data1E-13Sopt.csv"}, "plotOptiSig.png", "python3", "plot.py", "opti");


        //plot_with_python({"data0.8eps.csv", "data0.6eps.csv","data0.5eps.csv", "data0.4eps.csv", "data0.2eps.csv"}, "plot_eps_fit.png", "python3", "plot.py", "infid_fit_trace");
        //plot_with_python({"data0.8D2eps.csv", "data0.6D2eps.csv","data0.5D2eps.csv", "data0.4D2eps.csv"}, "plot_epsD2_fit.png", "python3", "plot.py", "infid_fit_trace");
        //plot_with_python({"data0.8Seps.csv", "data0.6Seps.csv","data0.5Seps.csv", "data0.4Seps.csv", "data0.2Seps.csv"}, "plot_eps_sig.png", "python3", "plot.py", "infid_trace");

    }


    return 0;
}



