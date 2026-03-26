//
// Created by Pascal Knoll on 25.03.26.
//
#include "OneSigmaSite.h"

OneSigmaSite::OneSigmaSite(double g)
    : BaseModel(4, "Qutit"){
    g_ = g;
    H0_  = diagOp4(s_, 0, pow(g_,2), pow(g_,2), pow(g_,2));

    H_1 = linearCombination4(s_, {
        {3/(4*pow(g_,2)*sqrt(3)), 'I', 'Y'},
        {3/(4*pow(g_,2)*sqrt(3)), 'Y', 'Z'},
        {3/(4*pow(g_,2)*sqrt(3)), 'Y', 'X'}});

    H_ = H_0;

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