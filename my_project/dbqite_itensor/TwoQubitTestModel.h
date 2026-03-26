//
// Created by Pascal Knoll on 15.03.26.
//

#ifndef MY_PROJECT_TWOQUBITTESTMODEL_H
#define MY_PROJECT_TWOQUBITTESTMODEL_H
#include "BaseModel2.h"
class TwoQubitTestModel : public BaseModel2 {
public:
    TwoQubitModel(OneQubitOp op1, OneQubitOp op2);
};
#endif //MY_PROJECT_TWOQUBITTESTMODEL_H