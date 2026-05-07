//
// Created by Pascal Knoll on 15.03.26.
//

#include "BuildingHamiltonians.h"

ITensor BuildingHamiltonians::make(const Index& s, OneQubitOp op){
    switch(op){
        case OneQubitOp::I:
            return idGate(s);

        case OneQubitOp::X:
            return opFromArray(s, 0.0, 1.0, 1.0, 0.0);

        case OneQubitOp::Y:
            return opFromArray(s, 0.0, cplx(0.0,-1.0), cplx(0.0,1.0), 0.0);

        case OneQubitOp::Z:
            return opFromArray(s, 1.0, 0.0, 0.0, -1.0);

        case OneQubitOp::Hd:
            return opFromArray(s, 1.0/sqrt(2.0), 1.0/sqrt(2.0), 1.0/sqrt(2.0), -1.0/sqrt(2.0));

        default:
            itensor::error("Unknown one-qubit operator");
            return ITensor();
    }
}

ITensor opFromMatrix4(Index const& s, array<cplx,16> const& M){
    if(dim(s) != 4)
        throw runtime_error("opFromMatrix4 requires Index of dimension 4");

    auto sp = prime(s);
    auto sd = dag(s);    // Store this once!

    ITensor Op(sp, sd);

    for(int r = 1; r <= 4; ++r)
        for(int c = 1; c <= 4; ++c)
        {
            Op.set(sp(r), sd(c), M[4*(r-1) + (c-1)]);  // Use stored sd
        }

    return Op;
}

ITensor hadamard4(Index const& s){
    double a = 0.5; // = 1/2
    return opFromMatrix4(s, {
         a,  a,  a,  a,
         a, -a,  a, -a,
         a,  a, -a, -a,
         a, -a, -a,  a
    });
}

array<cplx,4> pauli2(char which){
    switch(which){
        case 'I': return {1.0, 0.0, 0.0, 1.0};
        case 'X': return {0.0, 1.0, 1.0, 0.0};
        case 'Y': return {0.0, cplx(0.0,-1.0), cplx(0.0,1.0), 0.0};
        case 'Z': return {1.0, 0.0, 0.0, -1.0};
        case 'H':return {1.0/sqrt(2.0), 1.0/sqrt(2.0), 1.0/sqrt(2.0), -1.0/sqrt(2.0)};
        default:
            throw runtime_error("Unknown Pauli label");
    }
}

ITensor pauliTensor4(Index const& s, char l, char r){
    auto A = pauli2(l);
    auto B = pauli2(r);
    auto M = kron2x2(A, B);
    return opFromMatrix4(s, M);
}

ITensor linearCombination4(Index const& s, vector<PauliTerm> const& terms){
    ITensor H(dag(s), prime(s));

    for(auto const& t : terms){
        H += t.coeff * pauliTensor4(s, t.left, t.right);
    }

    return H;
}

ITensor diagOp4(Index const& s, double d1, double d2, double d3, double d4){
    return opFromMatrix4(s, {
        d1, 0.0, 0.0, 0.0,
        0.0, d2, 0.0, 0.0,
        0.0, 0.0, d3, 0.0,
        0.0, 0.0, 0.0, d4
    });
}


ITensor opFromArray2(Index const& s, cplx a00, cplx a01, cplx a10, cplx a11){
    auto sp = prime(s);
    auto sd = dag(s);    // Store once
    ITensor U(sp, sd);
    U.set(sp(1), sd(1), a00);
    U.set(sp(1), sd(2), a01);
    U.set(sp(2), sd(1), a10);
    U.set(sp(2), sd(2), a11);
    return U;
}

array<cplx,16> kron2x2(array<cplx,4> const& A,array<cplx,4> const& B){
    array<cplx,16> K{};

    // A(r1,c1), B(r2,c2) -> K(2*r1+r2, 2*c1+c2)
    for(int r1 = 0; r1 < 2; ++r1)
        for(int c1 = 0; c1 < 2; ++c1)
            for(int r2 = 0; r2 < 2; ++r2)
                for(int c2 = 0; c2 < 2; ++c2){
                    int row = 2*r1 + r2;
                    int col = 2*c1 + c2;
                    K[4*row + col] = A[2*r1 + c1] * B[2*r2 + c2];
                }

    return K;
}

