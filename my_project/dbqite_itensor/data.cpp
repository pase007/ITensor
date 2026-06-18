//
// Created by Pascal Knoll on 22.02.26.
//
#include "data.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>   // std::system
#include <cerrno>
#include <sys/stat.h>

namespace {
void ensure_directory(const string& path) {
    if(::mkdir(path.c_str(), 0775) != 0 && errno != EEXIST) {
        throw runtime_error("Could not create directory: " + path);
    }
}
}

// ---------------------------------------------------------------------------
// --------------------------- Save Data  Handling ---------------------------
// Put data into csv-file
void write_csv(const string& filename, const vector<Row>& rows){
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out){
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "k,Energy,Fidelity\n";
    out << setprecision(14); // keep good numeric precision

    for(const auto& r : rows){
        out << r.k << "," << r.energy << "," << r.fidelity << "\n";
    }
}

void write_csv_E(const string& filename, const vector<double>& Elevels){
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out){
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "Energy,\n";
    out << setprecision(14); // keep good numeric precision

    for(const auto& E : Elevels){
        out << E << "," << 0.0 << "\n";
    }
}

void write_csv_Espectrum(const string& filename, const vector<double> g_vals, const vector<vector<double>>& Espectrum){
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out){
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "g,Energy\n";
    out << setprecision(14); // keep good numeric precision
	for (int i = 0; i < Espectrum.size(); ++i) {
		out << g_vals[i];
		for(const auto& E : Espectrum[i]){
        	out << "," << E;
    	}
		out << "\n";
	}
}

void write_csv_EspectrumGap(const string& filename, const vector<double> g_vals, const vector<EnergyAnalysis>& Espectrum){
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out){
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "g,E0,E1,gap,var0,var1\n";
    out << setprecision(14); // keep good numeric precision
    for (int i = 0; i < g_vals.size(); ++i) {
        const auto& E = Espectrum.at(i);
        out << g_vals.at(i) << ","
            << E.E0 << ","
            << E.E1 << ","
            << E.gap << ","
            << E.var0 << ","
            << E.var1 << "\n";
    }
}

void write_csv_gapConvergence(const string& filename, const vector<double>& g_vals, const vector<int> N_vals, const vector<vector<EnergyAnalysis>>& Espectrum){
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out){
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "g,N,E0,E1,gap\n";
    out << setprecision(14); // keep good numeric precision
	for (int i = 0; i < g_vals.size()-1; i++) {
    	out << g_vals[i] << ",";
	}
	out << g_vals.back() << "\n";

	// --- ADDED variable-length gap convergence format ---
	// Row 2 is now: N_start,Nmax_for_g0,Nmax_for_g1,...
	// This allows each g curve to have a different number of N values.
	if(N_vals.empty()) {
		throw runtime_error("write_csv_gapConvergence: N_vals must contain at least the start value");
	}
	int N_start = N_vals.front();
	out << N_start;
	for(size_t i = 0; i < g_vals.size(); ++i) {
		int N_max = N_start + static_cast<int>(Espectrum.at(i).size()) - 1;
		out << "," << N_max;
	}
	out << "\n";
	// --- END ADDED variable-length gap convergence format ---

	for (int i = 0; i < Espectrum.size(); i++) {
		for (int j = 0; j < Espectrum[i].size(); j++) {
        	const auto& E = Espectrum[i].at(j);
        	out << E.E0 << "," << E.E1 << "," << E.gap << "\n";
    	}
	}
}

// Save function for optimization loop
void write_csv_opt(const string& filename, const vector<RowOpt>& rows) {
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out) {
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "s,k,Infidelity\n";
    out << setprecision(14);
    for(const auto& r : rows) {
        out << r.s_step << "," << r.k << "," << r.infidelity << "\n";
    }
}

// Save function for epsilon loop
void write_csv_eps(const string& filename, const vector<RowInfid>& rows) {
    const string data_dir = "./Data/";
    ofstream out(data_dir + filename);
    if(!out) {
        throw runtime_error("Could not open file for writing: " + filename);
    }
    out << "eps,k,Infidelity\n";
    out << setprecision(14);
    for(const auto& r : rows) {
        out << r.infid_step << "," << r.k << "," << r.infidelity << "\n";
    }
}

void write_csv_unitary_trace(const string& filename, const vector<DBQITEUnitaryTraceRow>& rows) {
    const string data_dir = "./DataCounts/";
    ensure_directory(data_dir);
    ofstream out(data_dir + filename);
    if(!out) {
        throw runtime_error("Could not open file for writing: " + filename);
    }

    out << "k,Infidelity,U0,U0_dag,A,A_dag,R0,R0_dag,total_unitaries\n";
    out << setprecision(14);
    for(const auto& r : rows) {
        out << r.k << ","
            << r.infidelity << ","
            << r.counts.U0 << ","
            << r.counts.U0_dag << ","
            << r.counts.A << ","
            << r.counts.A_dag << ","
            << r.counts.R0 << ","
            << r.counts.R0_dag << ","
            << r.counts.total() << "\n";
    }
}

void write_csv_unitary_opt_trace(const string& filename, const vector<DBQITEUnitaryOptTraceRow>& rows) {
    const string data_dir = "./DataCounts/";
    ensure_directory(data_dir);
    ofstream out(data_dir + filename);
    if(!out) {
        throw runtime_error("Could not open file for writing: " + filename);
    }

    out << "s,k,Infidelity,U0,U0_dag,A,A_dag,R0,R0_dag,total_unitaries\n";
    out << setprecision(14);
    for(const auto& r : rows) {
        out << r.s_step << ","
            << r.k << ","
            << r.infidelity << ","
            << r.counts.U0 << ","
            << r.counts.U0_dag << ","
            << r.counts.A << ","
            << r.counts.A_dag << ","
            << r.counts.R0 << ","
            << r.counts.R0_dag << ","
            << r.counts.total() << "\n";
    }
}


// ---------------------------------------------------------------------------
// --------------------------- Python Plot Handling --------------------------
// Plot energy and fidelity for single data file
int plot_with_python_S(const string& csv_file,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode){
    // Quote arguments to survive spaces in paths.
    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir ="./Plots/";

    string cmd =
        q(py_exec) + " " + q(script) + " " + mode + " " + q(plot_dir + out_png) + " " + q(data_dir + csv_file);

    cout << "\n[plot] Running: " << cmd << "\n";
    int rc = system(cmd.c_str());
    if(rc != 0){
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    }
    return rc;
}

// Overload for multiple plotted lines
int plot_with_python(const vector<string>& csv_files,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode){

    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir ="./Plots/";

    ostringstream cmd;
    cmd << q(py_exec) << " " << q(script) << " " << mode << " " << q(plot_dir + out_png);

    for(const auto& f : csv_files)
        cmd << " " << q(data_dir + f);

    cout << "\n[plot] Running: " << cmd.str() << "\n";
    int rc = system(cmd.str().c_str());
    if(rc != 0)
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    return rc;
}

// Overload for multiple plotted lines
int plot_with_python_E(const vector<string>& csv_files,
                     const string& csv_Efile,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode){

    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir ="./Plots/";
    const string plot_dirE ="./Plots/E";

    ostringstream cmd;
    cmd << q(py_exec) << " " << q(script) << " " << mode << " " << q(plot_dir + out_png) << " " << q(plot_dirE + out_png) << " " << q(data_dir + csv_Efile);

    for(const auto& f : csv_files)
        cmd << " " << q(data_dir + "data" + f + ".csv");

    cout << "\n[plot] Running: " << cmd.str() << "\n";
    int rc = system(cmd.str().c_str());
    if(rc != 0)
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    return rc;
}

int plot_Espectrum(const string& csv_file,
                     const string& s_file,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode){
    // Quote arguments to survive spaces in paths.
    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir ="./Plots/";

    string cmd =
        q(py_exec) + " " + q(script) + " " + mode + " " + q(plot_dir + out_png) + " " + q(data_dir + csv_file) + " " + q(data_dir + s_file);

    cout << "\n[plot] Running: " << cmd << "\n";
    int rc = system(cmd.c_str());
    if(rc != 0){
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    }
    return rc;
}

int plot_Espectrum_Pert(const string& csv_file,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode){
    // Quote arguments to survive spaces in paths.
    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir ="./Plots/";

    string cmd =
        q(py_exec) + " " + q(script) + " " + mode + " " + q(plot_dir + out_png) + " " + q(data_dir  + csv_file);

    cout << "\n[plot] Running: " << cmd << "\n";
    int rc = system(cmd.c_str());
    if(rc != 0){
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    }
    return rc;
}

int plot_Espectrum_gap(const string& csv_file,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode){
    // Quote arguments to survive spaces in paths.
    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir ="./Plots/";

    string cmd =
        q(py_exec) + " " + q(script) + " " + mode + " " + q(plot_dir + out_png) + " " + q(data_dir  + csv_file);

    cout << "\n[plot] Running: " << cmd << "\n";
    int rc = system(cmd.c_str());
    if(rc != 0){
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    }
    return rc;
}

int plot_mass_difference(const string& left_csv_file,
                     const string& right_csv_file,
                     const string& out_png,
                     const string& left_column,
                     const string& right_column,
                     const string& py_exec,
                     const string& script,
                     const string& mode){
    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";
    const string plot_dir = "./Plots/";

    string cmd =
        q(py_exec) + " " + q(script) + " " + mode + " " + q(plot_dir + out_png)
        + " " + q(data_dir + left_csv_file)
        + " " + q(data_dir + right_csv_file)
        + " " + q(left_column)
        + " " + q(right_column);

    cout << "\n[plot] Running: " << cmd << "\n";
    int rc = system(cmd.c_str());
    if(rc != 0){
        cerr << "[plot] Python plotting failed (exit code " << rc << ")\n";
    }
    return rc;
}

int merge_correlation_files(const string& out_csv_file,
                     const vector<string>& input_csv_files,
                     const string& py_exec,
                     const string& script,
                     const string& mode){
    if(input_csv_files.empty()) {
        cerr << "[mergeCorr] No input files provided\n";
        return 2;
    }

    auto q = [](const string& s){
        ostringstream os;
        os << "\"";
        for(char c : s) { if(c == '"') os << '\\'; os << c; }
        os << "\"";
        return os.str();
    };
    const string data_dir = "./Data/";

    ostringstream cmd;
    cmd << q(py_exec) << " " << q(script) << " " << mode << " " << q(data_dir + out_csv_file);
    for(const auto& input : input_csv_files) {
        cmd << " " << q(data_dir + input);
    }

    cout << "\n[mergeCorr] Running: " << cmd.str() << "\n";
    int rc = system(cmd.str().c_str());
    if(rc != 0){
        cerr << "[mergeCorr] Python merge failed (exit code " << rc << ")\n";
    }
    return rc;
}
