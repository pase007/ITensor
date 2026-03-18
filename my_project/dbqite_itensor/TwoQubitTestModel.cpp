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

    H_ = -(H1 * H2) - (I1 * I2);


    // target state |-> = (|0> - |1>) / sqrt(2)
    phi1_ = ITensor(s1_);
    phi1_.set(s1_(1),  1.0);
    phi1_.set(s1_(2), 0.0);
    phi2_ = ITensor(s2_);
    phi2_.set(s2_(1),  1.0);
    phi2_.set(s2_(2), 0.0);

    phi_ = phi1_ * phi2_;

    // |psi0> = |0>
    ket01_ = ITensor(s1_);
    ket01.set(s(1), 1.0/sqrt(2.0));
    ket01.set(s(2), 1.0/sqrt(2.0));
    ket02_ = ITensor(s2_);
    ket02.set(s(1), 1.0/sqrt(2.0));
    ket02.set(s(2), 1.0/sqrt(2.0));

    ket0_ = ket01_ * ket02_;
}


void TwoQubitTestModel::printSummary() const{
    cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    cout << "H = " << H_ << "\n";
    //cout << "phi = " << phi_ << "\n";
    //cout << "ket0 = " << ket0_ << "\n";
}

void TwoQubitTestModel::print_state() const {
    cout << "Site index: " << s1_ << " and " << s2_ << "\n";
    cout << "state Phi = " << phi_ << "\n";
}

void TwoQubitTestModel::printAsMatrix(ITensor const& A) const {
    // Collect ket and bra indices
    vector<Index> ket_inds;
    vector<Index> bra_inds;

    for(auto& i : inds(A)){
        if(primeLevel(i) == 0)
            ket_inds.push_back(i);
    }

    // flip order
    reverse(ket_inds.begin(), ket_inds.end());

    // build bra in same order
    for(auto& k : ket_inds){
        bra_inds.push_back(prime(k));
    }


    // Create combiners
    auto [cket, ket_combined] = combiner(ket_inds);
    auto [cbra, bra_combined] = combiner(bra_inds);

    // Combine indices
    auto M = cbra * (cket * A);
    cout << "M = " << M << "\n";
}

void TwoQubitTestModel::printAsVector(ITensor const& psi) const {
    vector<Index> ket_inds;

    for(auto const& i : inds(psi)) {
        if(primeLevel(i) == 0)
            ket_inds.push_back(i);
    }

    reverse(ket_inds.begin(), ket_inds.end());

    auto [cket, ket_combined] = combiner(ket_inds);
    auto V = cket * psi;

    cout << "V = " << V << "\n";
    PrintData(V);
}