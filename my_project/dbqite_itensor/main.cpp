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
#include "ChainModelMPS.h"
#include "GetData.h"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <future>
#include <stdexcept>
#include <sstream>
#include <string>

using namespace std;
///////////////////////// Main ////////////////////////////
namespace {

string safe_float_token(double x) {
	ostringstream os;
	os << fixed << setprecision(6) << x;
	string s = os.str();
	for(char& c : s) {
		if(c == '-') c = 'm';
		if(c == '.') c = 'p';
	}
	return s;
}

void print_usage(const char* exe) {
	cerr << "Usage:\n"
	     << "  " << exe << "\n"
	     << "  " << exe << " single-gap <g> <N> [output_dir]\n\n"
	     << "  " << exe << " gap-range <g> <N0> <N> [output_dir]\n\n"
	     << "  " << exe << " corr-scan <N> [comp_i comp_j] [gap]\n"
	     << "  " << exe << " corr-plot <N>\n\n"
	     << "Examples:\n"
	     << "  " << exe << " single-gap 0.7 48\n"
	     << "  " << exe << " single-gap 0.7 48 results/raw\n"
	     << "  " << exe << " gap-range 0.7 3 20 results/raw\n"
	     << "  " << exe << " corr-scan 50 2 2\n"
	     << "  " << exe << " corr-scan 50 2 2 gap\n"
	     << "  " << exe << " corr-plot 50\n";
}

EnergyAnalysis compute_gap(double g, int N) {
	if(N < 2) {
		throw invalid_argument("N must be at least 2 for the chain DMRG calculation");
	}

	SigmaChainAnalysis chain(g, N, false);
	return EnergyAnalysis{
		chain.groundEnergy(),
		chain.excEnergy(),
		chain.gap(),
		chain.groundVariance(),
		chain.excVariance()
	};
}

int run_single_gap(double g, int N, const string& output_dir) {
	cout << "Running single-gap job: g=" << g << ", N=" << N << "\n";
	const auto result = compute_gap(g, N);

	const string filename = "sigma_gap_g" + safe_float_token(g) + "_N" + to_string(N) + ".csv";
	const string path = output_dir + "/" + filename;
	const string tmp_path = path + ".tmp";

	ofstream out(tmp_path);
	if(!out) {
		throw runtime_error("Could not open output file for writing: " + tmp_path);
	}

	out << "g,N,E0,E1,gap,var0,var1\n";
	out << setprecision(14)
	    << g << ","
	    << N << ","
	    << result.E0 << ","
	    << result.E1 << ","
	    << result.gap << ","
	    << result.var0 << ","
	    << result.var1 << "\n";
	out.close();

	if(!out) {
		throw runtime_error("Failed while writing output file: " + tmp_path);
	}
	if(std::rename(tmp_path.c_str(), path.c_str()) != 0) {
		throw runtime_error("Could not move temporary output file to final path: " + path);
	}

	cout << "Wrote " << path << "\n";
	return 0;
}

int run_gap_range(double g, int N0, int N, const string& output_dir) {
	if(N0 < 2) {
		throw invalid_argument("N0 must be at least 2 for the chain DMRG calculation");
	}
	if(N < N0) {
		throw invalid_argument("N must be greater than or equal to N0");
	}

	cout << "Running gap range job: g=" << g << ", N=" << N0 << ".." << N << "\n";

	const string filename = "sigma_gap_g" + safe_float_token(g)
	                      + "_N" + to_string(N0) + "-" + to_string(N) + ".csv";
	const string path = output_dir + "/" + filename;
	const string tmp_path = path + ".tmp";

	ofstream out(tmp_path);
	if(!out) {
		throw runtime_error("Could not open output file for writing: " + tmp_path);
	}

	out << "g,N,E0,E1,gap,var0,var1\n";
	out << setprecision(14);
	for(int n = N0; n <= N; ++n) {
		cout << "  N=" << n << "\n";
		const auto result = compute_gap(g, n);
		out << g << ","
		    << n << ","
		    << result.E0 << ","
		    << result.E1 << ","
		    << result.gap << ","
		    << result.var0 << ","
		    << result.var1 << "\n";
		out.flush();
	}
	out.close();

	if(!out) {
		throw runtime_error("Failed while writing output file: " + tmp_path);
	}
	if(std::rename(tmp_path.c_str(), path.c_str()) != 0) {
		throw runtime_error("Could not move temporary output file to final path: " + path);
	}

	cout << "Wrote " << path << "\n";
	return 0;
}

} // namespace

int main(int argc, char* argv[]){
    cout << fixed << setprecision(14);
	if(argc > 1) {
		const string mode = argv[1];
		try {
			if(mode == "single-gap") {
				if(argc < 4 || argc > 5) {
					print_usage(argv[0]);
					return 2;
				}
				const double g = stod(argv[2]);
				const int N = stoi(argv[3]);
				const string output_dir = (argc == 5) ? argv[4] : "./Data";
				return run_single_gap(g, N, output_dir);
			}
			if(mode == "gap-range") {
				if(argc < 5 || argc > 6) {
					print_usage(argv[0]);
					return 2;
				}
				const double g = stod(argv[2]);
				const int N0 = stoi(argv[3]);
				const int N = stoi(argv[4]);
				const string output_dir = (argc == 6) ? argv[5] : "./Data";
				return run_gap_range(g, N0, N, output_dir);
			}
			if(mode == "corr-scan") {
				if(argc != 3 && argc != 5 && argc != 6) {
					print_usage(argv[0]);
					return 2;
				}
				const int N = stoi(argv[2]);
				const int comp_i = (argc >= 5) ? stoi(argv[3]) : 2;
				const int comp_j = (argc >= 5) ? stoi(argv[4]) : 2;
				const bool computeExcited = (argc == 6 && string(argv[5]) == "gap");
				scanCorrelationLengthDMRG(N, comp_i, comp_j, computeExcited);
				return 0;
			}
			if(mode == "corr-plot") {
				if(argc != 3) {
					print_usage(argv[0]);
					return 2;
				}
				const int N = stoi(argv[2]);
				plotCorrelationLengthDMRG(N);
				return 0;
			}

			print_usage(argv[0]);
			return 2;
		} catch(const exception& e) {
			cerr << "Error: " << e.what() << "\n";
			return 1;
		}
	}

	// --------------- Spectrum Analysis ---------------
	//scanSpectrum(4);
	//scanSpectrumDMRG(20);
	//scanSpectrumDMRG(4);
	//getChainGapConvergence();
	//scanCorrelationLengthDMRG(20, 2, 2);
	//plot_Espectrum("ExactEnergies4.txt", "OptimumSteps.txt", "Espectrum_gscan.png", "python3", "plot.py" );
	//scanCorrelationLengthDMRG(50, 2, 2);
	//plotCorrelationLengthDMRG(50);
	// Use the CLI instead, for example:
	//   ./dbqite_itensor corr-scan 50 2 2
	//   ./dbqite_itensor corr-plot 50



	// --------------- Forgot what that is ---------------
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
	bool MPS1 = false;
	bool MPSAdaptive = false;
    bool Plot = false;
	getData(Data1, Data2, DataSig2, DataChain, MPS1, MPSAdaptive, Plot);


    return 0;
}
