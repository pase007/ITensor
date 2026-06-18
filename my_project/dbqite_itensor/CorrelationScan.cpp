#include "CorrelationScan.h"
#include "SigmaChainAnalysis.h"
#include "data.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

void scanCorrelationLengthDMRGRange(int NrSites,
                                    int comp_i,
                                    int comp_j,
                                    double g0,
                                    double g_step,
                                    int gN,
                                    const string& out_file,
                                    bool computeExcited) {
	if(gN <= 0) {
		throw invalid_argument("gN must be positive");
	}
	string filename = out_file.empty()
	                ? "correlationLengths" + to_string(NrSites) + ".txt"
	                : out_file;
	const string data_dir = "./Data/";
	ofstream out(data_dir + filename);
	if(!out) {
		throw runtime_error("Could not open file for writing: " + filename);
	}

	out << setprecision(14);
	out << "g,N,E0,E1,gap,xi_gap,xi_corr,m_corr,i0,comp_i,comp_j,r,corr,abs_corr,log_abs_corr,rMin,rMax\n";

	int i0 = 25;//NrSites / 4;
	int rMin = 1;
	int rMax = NrSites-2*i0;//NrSites/2; //std::min(NrSites - i0, NrSites / 3);
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

void scanCorrelationLengthDMRG(int NrSites, int comp_i, int comp_j, bool computeExcited) {
	scanCorrelationLengthDMRGRange(
		NrSites,
		comp_i,
		comp_j,
		0.45,
		0.05,
		18,
		"correlationLengths" + to_string(NrSites) + ".txt",
		computeExcited
	);
}

void scanLaplaceBeltramiCorrelationLengthDMRG(int NrSites) {
	scanCorrelationLengthDMRG(NrSites, 0, 0);
}

void plotCorrelationLengthDMRG(int NrSites) {
	string filename = "correlationLengths" + to_string(NrSites) + ".txt";
	string plotname = "correlationLengths" + to_string(NrSites) + ".png";
	plot_Espectrum_gap(filename, plotname, "python", "plot.py", "corrProfile");
}
