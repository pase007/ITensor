//
// Created by Pascal Knoll on 22.02.26.
//

#include "data.h"
#include <fstream>
#include <iostream>
#include <cstdlib>   // std::system

void write_csv_kEF(std::string const& filename,
                   std::vector<int> const& k,
                   std::vector<double> const& E,
                   std::vector<double> const& F){
    if(k.size() != E.size() || k.size() != F.size())
    {
        throw std::runtime_error("write_csv_kEF: vector sizes do not match");
    }

    std::ofstream out(filename);
    if(!out) throw std::runtime_error("Could not open " + filename);

    out << "k,energy,fidelity\n";
    for(size_t i = 0; i < k.size(); ++i)
    {
        out << k[i] << "," << E[i] << "," << F[i] << "\n";
    }
    out.flush();
}

static std::string shell_quote(std::string s)
{
    // minimal POSIX single-quote escaping
    std::string r = "'";
    for(char c : s)
    {
        if(c == '\'') r += "'\\''";
        else r += c;
    }
    r += "'";
    return r;
}

int run_plot(std::string const& python_cmd,
             std::string const& script,
             std::string const& csv_file)
{
    // Call: python plot_matrix.py data.csv
    std::string cmd = python_cmd + " " + shell_quote(script) + " " + shell_quote(csv_file);
    std::cout << "Plot command: " << cmd << "\n";
    return std::system(cmd.c_str());
}