//
// Created by Pascal Knoll on 15.03.26.
//

#ifndef MY_PROJECT_BUILDINGHAMILTONIANS_H
#define MY_PROJECT_BUILDINGHAMILTONIANS_H
#pragma once
#include "itensor/all.h"
#include "functions.h"
#include <complex>
using namespace std;
using itensor::Index;
using itensor::ITensor;
using cplx = complex<double>;

// Building Hamiltonians
enum class OneQubitOp{
    I,
    X,
    Y,
    Z
};

struct BuildingHamiltonians{
    static ITensor make(const Index& s, OneQubitOp op);
};


#endif //MY_PROJECT_BUILDINGHAMILTONIANS_H