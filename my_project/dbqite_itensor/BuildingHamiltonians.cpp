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

        default:
            itensor::error("Unknown one-qubit operator");
            return ITensor();
    }
}