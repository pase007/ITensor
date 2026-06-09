//
// Created by Pascal Knoll on 07.05.26.
//

#include "GetData.h"
#include <fstream> // ADDED correlation-length data scan
#include <limits>  // ADDED correlation-length data scan

void getData(bool Data1, bool Data2, bool DataSig2, bool DataChain, bool MPS1, bool MPSAdaptive, bool Plot) {
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
		ChainModel chainN = ChainModel(1.0, 5, false);
		getChainN(chainN);
	}

	// MPS algos
	if (MPS1 == true) {
		getChainMPS(0.8, 7, false);
	}
	if (MPSAdaptive == true) {
		getChainMPSAdaptive(0.8, 7, false);
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
}

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
    vector<double> s_step_vec = { 0.19, 0.17, 0.15, 0.13, 0.11 };//, 0.25, 0.255, 0.26, 0.265, 0.275};
    vector<string> s_step_vec_str = {"0.19Ch2","0.17Ch2", "0.15Ch2", "0.13Ch2", "0.11Ch2"};//, "0.25Ch2", "0.255Ch2", "0.26Ch2", "0.265Ch2", "0.275Ch2"};

    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoLoopParams loop_params2;
        loop_params2.s_step = s_step_vec[i];
        loop_params2.s_string = s_step_vec_str[i];
        loop_params2.K = 100;
        loop_params2.infid_target = 1E-13;
        loop_params2.str_infid_target = "1E-13";
        //chain.basicModelLoop(loop_params2);
    }
    //plot_with_python_E(s_step_vec_str, "Elevels4.csv", "plotCh4205.png",  "python3", "plot.py", "multiE");

   // Params for optimization for given epsilon
    vector<double> e_step_vec = {1E-3, 1E-5, 1E-7, 1E-9};
    vector<string> e_step_vec_str = { "1E-2Ch2", "1E-5Ch2", "1E-7Ch2", "1E-9Ch2"};

	AlgoOptParamsVec opt_params2;
	opt_params2.s_min = 0.135;
	opt_params2.s_bin = 0.005;
	opt_params2.s_N = 14;
	opt_params2.K_max = 200;
	opt_params2.infid_target = e_step_vec;
	opt_params2.str_infid_target = e_step_vec_str;
	//chain.optimizeStepsLoopVec(opt_params2);
    //plot_with_python({ "data1E-2Ch2opt.csv", "data1E-5Ch2opt.csv", "data1E-7Ch2opt.csv", "data1E-9Ch2opt.csv"}, "plotOptiCh42post.png", "python3", "plot.py", "opti");

    // Params for reacing target infidelities
    for (int i = 0; i < s_step_vec.size(); i++) {
        AlgoInfidParams eps_params2;
        eps_params2.s_step = s_step_vec[i];
        eps_params2.s_string = s_step_vec_str[i];
        eps_params2.K_max = 200;
        eps_params2.infid_min = 1E-1;
        eps_params2.infid_N = 10;
        //chain.degradingInfidLoop(eps_params2);
    }
    //plot_with_python({"data0.2Ch2.csv", "data0.25Ch2.csv", "data0.3Ch2.csv", "data0.125Ch2.csv", "data0.175Ch2.csv", "data0.275Ch2.csv"}, "plotCh42.png",  "python3", "plot.py", "multi");
    //plot_with_python({  "data0.009Ch2eps.csv", "data0.01Ch2eps.csv", "data0.015Ch2eps.csv", "data0.11Ch2eps.csv", "data0.0095Ch2eps.csv"}, "plot_eps_Ch42post.png", "python3", "plot.py", "infid_fit_trace");

}



// ------------------------------------------ Spectrum and Energy and Gap scans ---------------------------------------
void scanSpectrum(int NrSites){
    string filename = "ExactEnergies4.txt";
	string filenameGap = "gapEnergiesDiag4.txt";

    vector<vector<double>> FullEnergies;
	vector<EnergyAnalysis> gap_energies;
	int N = 16;
	gap_energies.reserve(N);
    vector<double> g_vals;
    double g_step = 0.05;
    double g0 = 0.65;

    for (int i = 0; i < N; i++){
        double g = g0 + i*g_step;
        g_vals.push_back(g);
        ChainModel chainN = ChainModel(g, NrSites, false);

        vector<double> ei = chainN.ExactEnergies();
        //FullEnergies.push_back(ei);

    	double E0 = chainN.E0();
    	double E1 = chainN.E1();
    	double gap = chainN.gap();

    	gap_energies.push_back(EnergyAnalysis{E0, E1, gap});
    }
    //write_csv_Espectrum(filename, g_vals, FullEnergies);
	write_csv_EspectrumGap(filenameGap, g_vals, gap_energies);

	plot_Espectrum_gap(filenameGap, "SpectrumGapDiag4.png", "python3", "plot.py");
    //plot_Espectrum("ExactEnergies3.txt", "OptimumSteps3.txt", "E3spectrum_gscan.png", "python3", "plot.py" );
    //plot_Espectrum_Pert("ExactEnergies2.txt", "E2SpectrumPerturbation.png", "python3", "plot.py"  );
}

void scanSpectrumDMRG(int NrSites) {
    string filename = "gapEnergies" + to_string(NrSites) + ".txt";
    vector<double> g_vals;
    vector<EnergyAnalysis> gap_energies;
	int N = 16;
	gap_energies.reserve(N);

    double g_step = 0.05;
    double g0 = 0.65;
    for (int i = 0; i < N; i++) {
    	cout << "Step: " << i << endl;
        double g = g0 + i*g_step;
        g_vals.push_back(g);
        SigmaChainAnalysis chain(g, NrSites, false);
        double E0 = chain.groundEnergy();
        double E1 = chain.excEnergy();
        double gap = chain.gap();
        double var0 = chain.groundVariance();
        double var1 = chain.excVariance();
        gap_energies.push_back(EnergyAnalysis{E0, E1, gap, var0, var1});
    }
    write_csv_EspectrumGap(filename, g_vals, gap_energies);
    plot_Espectrum_gap(filename, "SpectrumGap" + to_string(NrSites) + ".png", "python3", "plot.py");

}

//Main code to optain mass gap for various G as N-->inf
void getChainGapConvergence() {
	string filenameGap = "gapConv.txt";
	string filenameGapPlot = "gapConv3.txt";
	bool data = !true;
	int Nval = 63;
	int N0 = 3;
	double g0 = 0.7;
	double g_step = 0.1;
	double g_max = 1.3;
	int g_nr = 2;//int((g_max - g0)/g_step)+1;

	vector<int> sysN;
	vector<double> g_vals;
	vector<vector<EnergyAnalysis>> gap_energies(g_nr);

	for (int i = 0; i < g_nr; i++) {
		if (!data) continue;
		double g = g0 + i*g_step;
		g_vals.push_back(g);
		cout << "--------- g = " << g << " ----------" << "\n" << "----------------------------------------" << endl;
		for (int N = N0; N <= Nval; N++) {
			cout << "system size: N = " << N << endl;
			if (g==g0) sysN.push_back(N);
			SigmaChainAnalysis chain(g, N, false);
			double E0 = chain.groundEnergy();
			double E1 = chain.excEnergy();
			double gap = chain.gap();
			gap_energies[i].push_back(EnergyAnalysis{E0, E1, gap});
		}
	}
	if (data) write_csv_gapConvergence(filenameGap, g_vals, sysN, gap_energies);
	if (!data) {
		plot_Espectrum_gap(filenameGapPlot, "gapConvergence2_fit.png", "python3", "plot.py", "gapConvergenceFit");
		plot_Espectrum_gap(filenameGapPlot, "energyConvergence.png", "python3", "plot.py", "energyConvergence");
	}
}

// --- ADDED correlation-length data scan ---
void scanCorrelationLengthDMRG(int NrSites, int comp_i, int comp_j, bool computeExcited) {
	string filename = "correlationLengths" + to_string(NrSites) + ".txt";
	const string data_dir = "./Data/";
	ofstream out(data_dir + filename);
	if(!out) {
		throw runtime_error("Could not open file for writing: " + filename);
	}

	out << setprecision(14);
	out << "g,N,E0,E1,gap,xi_gap,xi_corr,m_corr,i0,comp_i,comp_j,r,corr,abs_corr,log_abs_corr,rMin,rMax\n";

	int gN = 8;
	double g_step = 0.05;
	double g0 = 0.85;

	int i0 = NrSites / 4;
	int rMin = 1;
	int rMax = NrSites/2; //std::min(NrSites - i0, NrSites / 3);
	if(rMax < rMin) {
		rMin = 1;
		rMax = NrSites - i0;
	}

	for(int step = 0; step < gN; ++step) {
		double g = g0 + step * g_step;
		cout << "CorrLength: current g = " << g << endl;
		cout << "  excited-state DMRG: "
		     << (computeExcited ? "on" : "off")
		     << endl;

		SigmaChainAnalysis chain(g, NrSites, false, true, computeExcited);
		double E0 = chain.groundEnergy();
		double E1 = chain.excEnergy();
		double gap = chain.gap();
		double xi_gap = (std::isfinite(gap) && gap != 0.0)
		              ? 1.0 / gap
		              : std::numeric_limits<double>::quiet_NaN();
		double xi_corr = std::numeric_limits<double>::quiet_NaN();

		vector<std::pair<int,double>> profile;
		for(int r = 1; r <= rMax; ++r) {
			if(r == 1 || r % 2 == 0 || r == rMax) {
				cout << "  measuring  r = " << r << "/" << rMax << " ...";
			}

			double corr = chain.connectedCorr(i0, i0 + r, comp_i, comp_j);
			profile.push_back({r, corr});

			double abs_corr = std::abs(corr);
			double log_abs_corr = (abs_corr > 0.0)
			                    ? std::log(abs_corr)
			                    : std::numeric_limits<double>::quiet_NaN();

			if(r == 1 || r % 2 == 0 || r == rMax) {
				cout << "  done C(r) = " << corr << endl;
			}
		}

		try {
			double sx = 0.0;
			double sy = 0.0;
			double sxx = 0.0;
			double sxy = 0.0;
			int n = 0;

			for(auto const& [r, corr] : profile) {
				if(r < rMin || r > rMax) continue;
				double abs_corr = std::abs(corr);
				if(abs_corr <= 1E-14 || !std::isfinite(abs_corr)) continue;

				double x = static_cast<double>(r);
				double y = std::log(abs_corr);
				sx += x;
				sy += y;
				sxx += x * x;
				sxy += x * y;
				++n;
			}

			if(n < 2) {
				throw runtime_error("not enough nonzero data points");
			}

			double denom = n * sxx - sx * sx;
			if(std::abs(denom) <= 1E-14) {
				throw runtime_error("singular linear fit");
			}

			double slope = (n * sxy - sx * sy) / denom;
			if(slope >= 0.0) {
				throw runtime_error("fitted slope is non-negative");
			}

			xi_corr = -1.0 / slope;
			cout << "  fitted xi_corr = " << xi_corr
			     << ", m_corr = " << 1.0 / xi_corr
			     << endl;
		} catch(std::exception const& e) {
			cout << "fitCorrelationLength failed for g = " << g << ": " << e.what() << endl;
		}
		double m_corr = (std::isfinite(xi_corr) && xi_corr != 0.0)
		              ? 1.0 / xi_corr
		              : std::numeric_limits<double>::quiet_NaN();

		for(auto const& [r, corr] : profile) {
			double abs_corr = std::abs(corr);
			double log_abs_corr = (abs_corr > 0.0)
			                    ? std::log(abs_corr)
			                    : std::numeric_limits<double>::quiet_NaN();

			out << g << ","
			    << NrSites << ","
			    << E0 << ","
			    << E1 << ","
			    << gap << ","
			    << xi_gap << ","
			    << xi_corr << ","
			    << m_corr << ","
			    << i0 << ","
			    << comp_i << ","
			    << comp_j << ","
			    << r << ","
			    << corr << ","
			    << abs_corr << ","
			    << log_abs_corr << ","
			    << rMin << ","
			    << rMax << "\n";
		}
		out.flush();
	}
}

void scanLaplaceBeltramiCorrelationLengthDMRG(int NrSites) {
	scanCorrelationLengthDMRG(NrSites, 0, 0);
}

void plotCorrelationLengthDMRG(int NrSites) {
	string filename = "correlationLengths" + to_string(NrSites) + ".txt";
	string plotname = "correlationLengths" + to_string(NrSites) + ".png";
	plot_Espectrum_gap(filename, plotname, "python", "plot.py", "corrProfile");
}
// --- END ADDED correlation-length data scan ---




// ------------------------------------------ MPS Code: Shedule and Optimization Loop ---------------------------------------
void getChainMPS(double g, int N, bool PBC) {
	// Params for Fidelity in time Algo
	vector<double> s_step_vec = {0.15, 0.18, 0.2, 0.24};//, 0.25, 0.255, 0.26, 0.265, 0.275};
	vector<string> s_step_vec_str = {"0.15MPS", "0.18MPS", "0.2MPS", "0.24MPS" };//, "0.25Ch2", "0.255Ch2", "0.26Ch2", "0.265Ch2", "0.275Ch2"};
	double schedule_factor = 0.9;

	if (s_step_vec.size() != s_step_vec_str.size()) {
		throw runtime_error("getChainMPS: s_step_vec and s_step_vec_str must have the same length");
	}

	for (size_t i = 0; i < s_step_vec.size(); i++) {
		const double s_step = s_step_vec[i];
		const string s_string = s_step_vec_str[i];

		ChainModelMPS local_mps(g, N, PBC);

		AlgoLoopParams loop_params2;
		loop_params2.s_step = s_step;
		loop_params2.s_string = s_string;
		loop_params2.K = 6;
		loop_params2.infid_target = 1E-5;
		loop_params2.str_infid_target = "1E-5";
		loop_params2.schedule_factor = schedule_factor;

		local_mps.basicModelLoop(loop_params2);
	}
    plot_with_python_E(s_step_vec_str, "Elevels6.csv", "plotMPS4.png",  "python3", "plot.py", "multiE");
}

void getChainMPSAdaptive(double g, int N, bool PBC) {
	vector<double> s_step_vec = {0.18, 0.26, 0.35};
	vector<string> s_step_vec_str = {"0.18MPSopt", "0.26MPSopt", "0.35MPSopt"};
	vector<double> s_candidates = {};// {0.05, 0.08, 0.1, 0.12, 0.15, 0.18, 0.2, 0.24};
	bool refine_s = false;

	if (s_step_vec.size() != s_step_vec_str.size()) {
		throw runtime_error("getChainMPSAdaptive: s_step_vec and s_step_vec_str must have the same length");
	}

	for (size_t i = 0; i < s_step_vec.size(); ++i) {
		ChainModelMPS chainMPS(g, N, PBC);

		AlgoLoopParams loop_params;
		loop_params.s_step = s_step_vec[i];
		loop_params.s_string = s_step_vec_str[i];
		loop_params.K = 5;
		loop_params.infid_target = 1E-5;
		loop_params.str_infid_target = "1E-5";
		loop_params.s_candidates = s_candidates;
		loop_params.refine_s = refine_s;

		chainMPS.adaptiveSModelLoop(loop_params);
	}

	plot_with_python_E(s_step_vec_str, "Elevels6.csv", "plotAdaptiveMPS.png",  "python3", "plot.py", "multiE");
}
