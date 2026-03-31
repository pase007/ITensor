//
// Created by Pascal Knoll on 25.03.26.
//
#include "OneSigmaSite.h"

OneSigmaSite::OneSigmaSite(double g)
    : BaseModel(4, "Qutit"){
    g_ = g;
    Id_s_ = diagOp4(s_, 1.0, 1.0, 1.0, 1.0);

    // Hamiltonian
    H0_  = diagOp4(s_, 0, pow(g_,2), pow(g_,2), pow(g_,2));
    j1_ = 1/sqrt(3) * pauliTensor4(s_, 'I', 'Y');
    j2_ = 1/sqrt(3) * pauliTensor4(s_, 'Y', 'Z');
    j3_ = 1/sqrt(3) * pauliTensor4(s_, 'Y', 'X');

    H_ = H0_;

    // Start state
    phi_ = ITensor(s_);
    phi_.set(s_(1),  1.0 / 2.0);
    phi_.set(s_(2), -1.0 / 2.0);
    phi_.set(s_(3),  1.0 / 2.0);
    phi_.set(s_(4), -1.0 / 2.0);

    // taget: ground state
    ket0_ = ITensor(s_);
    ket0_.set(s_(1), 1.0);
    ket0_.set(s_(2), 0.0);
    ket0_.set(s_(3), 0.0);
    ket0_.set(s_(4), 0.0);

    checkHamiltonian();
}