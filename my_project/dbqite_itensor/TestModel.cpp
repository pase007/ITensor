//
// Created by Pascal Knoll on 14.03.26.
//

#include "TestModel.h"
#include <iostream>
#include <cmath>

#include "TestModel.h"

TestModel::TestModel()
    : BaseModel(2, "Qubit"){
    H_  = opFromArray(s_, 0.0, cplx(1.0, 0.0),  cplx(1.0, 0.0), 0.0);
    H2_ = opFromArray(s_, 0.0, cplx(0.0, -1.0), cplx(0.0, 1.0), 0.0);
    H3_ = opFromArray(s_, cplx(1.0, 0.0), 0.0, 0.0, cplx(-1.0, 0.0));

	// Start state
    phi_ = ITensor(s_);
    phi_.set(s_(1),  1.0 / std::sqrt(2.0));
    phi_.set(s_(2), -1.0 / std::sqrt(2.0));

	// taget: ground state
    ket0_ = ITensor(s_);
    ket0_.set(s_(1), 1.0);
    ket0_.set(s_(2), 0.0);

    checkHamiltonian();
}




