//
// Created by Pascal Knoll on 15.03.26.
//

#ifndef MY_PROJECT_BUILDINGHAMILTONIANS_H
#define MY_PROJECT_BUILDINGHAMILTONIANS_H
#pragma once
#include "itensor/all.h"
#include "functions.h"
#include "StructFile.h"
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
    Z,
    Hd
};

ITensor opFromArray2(Index const& s, cplx a00, cplx a01, cplx a10, cplx a11);
ITensor opFromMatrix4(Index const& s, array<cplx,16> const& M);
ITensor hadamard4(Index const& s);

struct BuildingHamiltonians{
    static ITensor make(const Index& s, OneQubitOp op);
};
ITensor pauliTensor4(Index const& s, char l, char r);
ITensor diagOp4(Index const& s, double d1, double d2, double d3, double d4);
ITensor linearCombination4(Index const& s, vector<PauliTerm> const& terms);

array<cplx,4> pauli2(char which);
array<cplx,16> kron2x2(array<cplx,4> const& A, array<cplx,4> const& B);

#endif //MY_PROJECT_BUILDINGHAMILTONIANS_H