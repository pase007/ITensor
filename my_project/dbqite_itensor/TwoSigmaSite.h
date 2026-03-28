//
// Created by Pascal Knoll on 26.03.26.
//

#ifndef MY_PROJECT_TWOSIGMASITE_H
#define MY_PROJECT_TWOSIGMASITE_H
#include "BaseModel2.h"
class TwoSigmaSite : public BaseModel2 {
protected:
    double g_;

    ITensor j1_s1_;
    ITensor j1_s2_;
    ITensor j2_s1_;
    ITensor j2_s2_;
    ITensor j3_s1_;
    ITensor j3_s2_;

    ITensor H_Ik_;

    ITensor H0_s1_;
    ITensor H0_s2_;

public:
    TwoSigmaSite(double g);
    const double g() const { return g_; }
    void buildHamiltonian_();
};
#endif //MY_PROJECT_TWOSIGMASITE_H