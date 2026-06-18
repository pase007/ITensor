#ifndef MY_PROJECT_CORRELATIONSCAN_H
#define MY_PROJECT_CORRELATIONSCAN_H

#include <string>

void scanCorrelationLengthDMRG(int NrSites,
                               int comp_i = 2,
                               int comp_j = 1,
                               bool computeExcited = false);
void scanCorrelationLengthDMRGRange(int NrSites,
                               int comp_i,
                               int comp_j,
                               double g0,
                               double g_step,
                               int gN,
                               const std::string& out_file,
                               bool computeExcited = false);
void scanLaplaceBeltramiCorrelationLengthDMRG(int NrSites);
void plotCorrelationLengthDMRG(int NrSites);

#endif // MY_PROJECT_CORRELATIONSCAN_H
