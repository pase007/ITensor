//
// Created by Pascal Knoll on 26.03.26.
//

#include "TwoSigmaSite.h"
TwoSigmaSiteModel::TwoSigmaSiteModel(double g)
    : BaseModel2(Index(4,"Site1"), Index(4,"Site2")){
    g_ = g;
    H0_  = diagOp4(s_, 0, pow(g_,2), pow(g_,2), pow(g_,2));


}