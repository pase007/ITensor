//
// Created by Pascal Knoll on 22.02.26.
//

#ifndef MY_PROJECT_DATA_H
#define MY_PROJECT_DATA_H
#pragma once
#include "StructFile.h"
#include <string>
#include <vector>
using namespace std;

// Write CSV with header: k,Energy,Fidelity
void write_csv(const string& filename, const vector<Row>& rows);
void write_csv_E(const string& filename, const vector<double>& Elevels);
void write_csv_Espectrum(const string& filename, const vector<double> g_vals, const vector<vector<double>>& Espectrum);
void write_csv_EspectrumGap(const string& filename, const vector<double> g_vals, const vector<EnergyAnalysis>& Espectrum);
void write_csv_gapConvergence(const string& filename, const vector<double>& g_vals, const vector<int> N_vals, const vector<vector<EnergyAnalysis>>& Espectrum);
void write_csv_opt(const string& filename, const vector<RowOpt>& rows);
void write_csv_eps(const string& filename, const vector<RowInfid>& rows);

// Run python plotter (default: python3 plot.py data.csv plot.png)
int plot_with_python_S(const string& csv_file,
                     const string& out_png = "plot.png",
                     const string& py_exec = "python3",
                     const string& script  = "plot.py",
                     const string& mode = "single");

int plot_with_python(const vector<string>& csv_files,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode = "multi");

int plot_with_python_E(const vector<string>& csv_files,
                     const string& csv_Efile,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode);

int plot_Espectrum(const string& csv_file,
                     const string& s_file,
                     const string& out_png = "plot.png",
                     const string& py_exec = "python3",
                     const string& script  = "plot.py",
                     const string& mode = "spectrum");

int plot_Espectrum_Pert(const string& csv_file,
                     const string& out_png = "plot.png",
                     const string& py_exec = "python3",
                     const string& script  = "plot.py",
                     const string& mode = "spectrumPert");

int plot_Espectrum_gap(const string& csv_file,
                     const string& out_png = "plot.png",
                     const string& py_exec = "python3",
                     const string& script  = "plot.py",
                     const string& mode = "spectrumEgap");
#endif //MY_PROJECT_DATA_H