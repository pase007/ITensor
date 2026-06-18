//
// Created by Pascal Knoll on 07.05.26.
//

#ifndef MY_PROJECT_GETDATA_H
#define MY_PROJECT_GETDATA_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <future>
#include <stdexcept>
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
#include "SigmaChainAnalysis.h"
using namespace std;

void getData1(TestModel& model);
void getData2(TwoQubitTestModel& model2);
void getSigma2(TwoSigmaSite& sigma2);
void getChainN(ChainModel& chain);
void getData(bool Data1, bool Data2, bool DataSig2, bool DataChain, bool MPS1, bool MPSAdaptive, bool Plot);

// Energy and GAp Analysis
void scanSpectrumDMRG(int NrSites);
void scanSpectrum(int NrSites);
void getChainGapConvergence();
// MPS based code with different theta shedules
void getChainMPS(double g, int N, bool PBC);
void getChainMPSAdaptive(double g, int N, bool PBC);


#endif //MY_PROJECT_GETDATA_H
