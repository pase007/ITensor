//
// Created by Pascal Knoll on 25.03.26.
//

#ifndef MY_PROJECT_ONESIGMASITE_H
#define MY_PROJECT_ONESIGMASITE_H
#include "BaseModel.h"

class OneSigmaSite : public BaseModel {
protected:
    double g_;
    ITensor j1_s_;
    ITensor j2_s_;
    ITensor j3_s_;
public:
    OneSigmaSite(double g);
    const double g() const { return g_; }
};

#endif //MY_PROJECT_ONESIGMASITE_H
