//
// Created by Pascal Knoll on 15.03.26.
//
#include "TwoQubitTestModel.h"

TwoQubitTestModel::TwoQubitTestModel(OneQubitOp op1, OneQubitOp op2)
    : s1_(2,"Qubit1"), s2_(2,"Qubit2"){
    ITensor H1 = BuildingHamiltonians::make(s1_, op1);
    ITensor H2 = BuildingHamiltonians::make(s2_, op2);
    ITensor I1 = BuildingHamiltonians::make(s1_, op2);
    ITensor I2 = BuildingHamiltonians::make(s2_, op1);

    H_ = ((H1 * H2) + 2*(I1 * I2))/2;
    //H_ = transformToMatrix(H_);
    //cout << "transformed H" << endl;


    // taget: ground state |-->
    ket01_ = ITensor(s1_);
    ket01_.set(s1_(1), 1.0/sqrt(2.0));
    ket01_.set(s1_(2), -1.0/sqrt(2.0));
    ket02_ = ITensor(s2_);
    ket02_.set(s2_(1), 1.0/sqrt(2.0));
    ket02_.set(s2_(2), -1.0/sqrt(2.0));

    ket0_ = ket01_ * ket02_;
    //phi_ = transformToVector(phi_);
    //cout << "transformed phi" << endl;

    // start state |psi0> = |00>
    psi1_ = ITensor(s1_);
    psi1_.set(s1_(1), 1.0);
    psi1_.set(s1_(2), 0.0);
    psi2_ = ITensor(s2_);
    psi2_.set(s2_(1), 1.0);
    psi2_.set(s2_(2), 0.0);

    psi_ = psi1_ * psi2_;
}


void TwoQubitTestModel::printSummary() const{
    cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    cout << "H = " << H_ << "\n";
    //cout << "phi = " << phi_ << "\n";
    //cout << "ket0 = " << ket0_ << "\n";
}

void TwoQubitTestModel::print_state() const {
    //cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    //cout << "state Phi = " << phi_ << "\n";
}

ITensor TwoQubitTestModel::transformToMatrix(ITensor const& A) const {
    vector<Index> ket_inds;
    vector<Index> bra_inds;

    for(auto const& i : inds(A)){
        if(primeLevel(i) == 0)
            ket_inds.push_back(i);
    }

    reverse(ket_inds.begin(), ket_inds.end());

    for(auto const& k : ket_inds){
        bra_inds.push_back(prime(k));
    }

    auto [cket, ket_combined] = combiner(ket_inds);
    auto [cbra, bra_combined] = combiner(bra_inds);

    ITensor M = cbra * (cket * A);

    // only prime the first combined index
    auto bra_p = prime(bra_combined);
    M.replaceInds({bra_combined}, {bra_p});

    //cout << "bra prime level = " << primeLevel(bra_p) << "\n";
    //cout << "ket prime level = " << primeLevel(ket_combined) << "\n";

    //cout << inds(M) << endl;
    //cout << "M = " << M << "\n";

    return M;
}

ITensor TwoQubitTestModel::transformToVector(ITensor const& psi) const {
    vector<Index> ket_inds;

    for(auto const& i : inds(psi)) {
        if(primeLevel(i) == 0)
            ket_inds.push_back(i);
    }

    reverse(ket_inds.begin(), ket_inds.end());

    auto [cket, ket_combined] = combiner(ket_inds);
    auto V = cket * psi;

    //cout << "V = " << V << "\n";
    //cout << inds(V) << endl;
    //cout << "V prime level = " << primeLevel(inds(V)[0]) << "\n";
    return V;
}