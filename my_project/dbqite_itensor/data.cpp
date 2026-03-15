//
// Created by Pascal Knoll on 22.02.26.
//
#include "data.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>   // std::system

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


// ---------------------------------------------------------------------------
// --------------------------- Python Plot Handling --------------------------
// Plot energy and fidelity for single data file
int plot_with_python(const string& csv_file,
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
