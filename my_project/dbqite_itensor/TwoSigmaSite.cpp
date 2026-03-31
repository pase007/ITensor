//
// Created by Pascal Knoll on 26.03.26.
//

#include "TwoSigmaSite.h"
TwoSigmaSite::TwoSigmaSite(double g)
    : BaseModel2(Index(4,"Site1"), Index(4,"Site2")){
    g_ = g;

    // Local free Hamiltonian
    H0_s1_  = diagOp4(s1_, 0, pow(g_,2), pow(g_,2), pow(g_,2));
    H0_s2_  = diagOp4(s2_, 0, pow(g_,2), pow(g_,2), pow(g_,2));

    // Interaction parts for each site
    j1_s1_ = 1/sqrt(3) * pauliTensor4(s1_, 'I', 'Y');
    j2_s1_ = 1/sqrt(3) * pauliTensor4(s1_, 'Y', 'Z');
    j3_s1_ = 1/sqrt(3) * pauliTensor4(s1_, 'Y', 'X');
    j1_s2_ = 1/sqrt(3) * pauliTensor4(s2_, 'I', 'Y');
    j2_s2_ = 1/sqrt(3) * pauliTensor4(s2_, 'Y', 'Z');
    j3_s2_ = 1/sqrt(3) * pauliTensor4(s2_, 'Y', 'X');

    // Build the full Hamiltonian H0 + Hi
    buildHamiltonian_();
    U0_ = hadamard4(s1_) * hadamard4(s2_);

    // Start state psi0 = |++>
    psi1_ = ITensor(s1_);
    psi2_ = ITensor(s2_);
    psi1_.set(s1_(1), 1.0);
    psi1_.set(s1_(2), 0.0);
    psi1_.set(s1_(3), 0.0);
    psi1_.set(s1_(4), 0.0);
    psi2_.set(s2_(1), 1.0);
    psi2_.set(s2_(2), 0.0);
    psi2_.set(s2_(3), 0.0);
    psi2_.set(s2_(4), 0.0);

    psi0_ = psi1_ * psi2_;

    // taget: ground state
    pair<double, vector<ITensor>> gs = groundSpace2(1E-14);
    E0_   = gs.first;
    groundspace_ = gs.second;
    checkGroundStateDegeneracy(1E-14);

    // Reflection state |00>
    ITensor p1 = ITensor(s1_);
    ITensor p2 = ITensor(s2_);
    p1.set(s1_(1), 1.0);
    p1.set(s1_(2), 0.0);
    p1.set(s1_(3), 0.0);
    p1.set(s1_(4), 0.0);
    p2.set(s2_(1), 1.0);
    p2.set(s2_(2), 0.0);
    p2.set(s2_(3), 0.0);
    p2.set(s2_(4), 0.0);
    p0_ = p1 * p2;

}


void TwoSigmaSite::buildHamiltonian_() {
    ITensor I1 = diagOp4(s1_, 1.0, 1.0, 1.0, 1.0);
    ITensor I2 = diagOp4(s2_, 1.0, 1.0, 1.0, 1.0);
    ITensor H0 = H0_s1_ * I2 + I1 * H0_s2_;
    H_Ik_ = -3.0/(4.0*pow(g_,2))*(j1_s1_*j1_s2_ + j2_s1_*j2_s2_ + j3_s1_*j3_s2_);
    H_ =  H0 + H_Ik_;
}

