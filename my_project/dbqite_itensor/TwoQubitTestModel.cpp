//
// Created by Pascal Knoll on 15.03.26.
//
#include "TwoQubitTestModel.h"

TwoQubitTestModel::TwoQubitTestModel(OneQubitOp op1, OneQubitOp op2)
    : BaseModel2(Index(2,"Qubit1"), Index(2,"Qubit2")){
    ITensor H1 = BuildingHamiltonians::make(s1_, op1);
    ITensor H2 = BuildingHamiltonians::make(s2_, op2);
    ITensor I1 = BuildingHamiltonians::make(s1_, op2);
    ITensor I2 = BuildingHamiltonians::make(s2_, op1);

    H_ = ((H1 * H2) + (I1 * I2))/2;
    //cout << "H = " << H_ << endl;


    // taget: ground state |-->
    ket01_ = ITensor(s1_);
    ket01_.set(s1_(1), 1.0/sqrt(2.0));
    ket01_.set(s1_(2), -1.0/sqrt(2.0));
    ket02_ = ITensor(s2_);
    ket02_.set(s2_(1), 1.0/sqrt(2.0));
    ket02_.set(s2_(2), -1.0/sqrt(2.0));

    ket0_ = ket01_ * ket02_;

    // start state |psi0> = |00>
    psi1_ = ITensor(s1_);
    psi1_.set(s1_(1), 1.0);
    psi1_.set(s1_(2), 0.0);
    psi2_ = ITensor(s2_);
    psi2_.set(s2_(1), 1.0);
    psi2_.set(s2_(2), 0.0);

    psi0_ = psi1_ * psi2_;
}























