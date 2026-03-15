//
// Created by Pascal Knoll on 22.02.26.
//

#ifndef MY_PROJECT_DATA_H
#define MY_PROJECT_DATA_H
#pragma once
#include <string>
#include <vector>

void write_csv_kEF(std::string const& filename,
                   std::vector<int> const& k,
                   std::vector<double> const& E,
                   std::vector<double> const& F);

int run_plot(std::string const& python_cmd,
             std::string const& script,
             std::string const& csv_file);
#endif //MY_PROJECT_DATA_H