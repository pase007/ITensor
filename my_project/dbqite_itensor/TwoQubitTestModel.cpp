//
// Created by Pascal Knoll on 15.03.26.
//
#include "TwoQubitTestModel.h"

TwoQubitTestModel::TwoQubitTestModel(OneQubitOp op1, OneQubitOp op2)
    : s1_(2,"Qubit1"), s2_(2,"Qubit2"){
    ITensor H1 = BuildingHamiltonians::make(s1_, op1);
    ITensor H2 = BuildingHamiltonians::make(s2_, op2);

    H_ = H1 * H2;
}