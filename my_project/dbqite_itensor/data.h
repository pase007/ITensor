//
// Created by Pascal Knoll on 22.02.26.
//

#ifndef MY_PROJECT_DATA_H
#define MY_PROJECT_DATA_H
#pragma once
#include <string>
#include <vector>
using namespace std;

struct Row{
    int k;
    double energy;
    double fidelity;
};

struct RowOpt {
    double s_step;
    int k;
    double infidelity;
};


// Write CSV with header: k,Energy,Fidelity
void write_csv(const string& filename, const vector<Row>& rows);
void write_csv_opt(const string& filename, const vector<RowOpt>& rows);

// Run python plotter (default: python3 plot.py data.csv plot.png)
int plot_with_python(const string& csv_file,
                     const string& out_png = "plot.png",
                     const string& py_exec = "python3",
                     const string& script  = "plot.py",
                     const string& mode = "single");

int plot_with_python(const vector<string>& csv_files,
                     const string& out_png,
                     const string& py_exec,
                     const string& script,
                     const string& mode = "multi");


#endif //MY_PROJECT_DATA_H