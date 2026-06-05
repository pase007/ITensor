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
	     << "Examples:\n"
	     << "  " << exe << " single-gap 0.7 48\n"
	     << "  " << exe << " single-gap 0.7 48 results/raw\n";
}

int run_single_gap(double g, int N, const string& output_dir) {
	if(N < 2) {
		throw invalid_argument("N must be at least 2 for the chain DMRG calculation");
	}

	cout << "Running single-gap job: g=" << g << ", N=" << N << "\n";
	SigmaChainAnalysis chain(g, N, false);

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
	    << chain.groundEnergy() << ","
	    << chain.excEnergy() << ","
	    << chain.gap() << ","
	    << chain.groundVariance() << ","
	    << chain.excVariance() << "\n";
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
	getChainGapConvergence();
	//scanCorrelationLengthDMRG(20, 2, 2);
	//plot_Espectrum("ExactEnergies4.txt", "OptimumSteps.txt", "Espectrum_gscan.png", "python3", "plot.py" );




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
