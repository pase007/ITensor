//
// Created by Pascal Knoll on 15.03.26.
//

#ifndef MY_PROJECT_STRUCTFILE_H
#define MY_PROJECT_STRUCTFILE_H
#include "itensor/all.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
using namespace std;
using cplx = complex<double>;

struct BigCount {
    string value = "0";

    BigCount() = default;
    BigCount(unsigned long long n) : value(to_string(n)) {}
    explicit BigCount(string v) : value(std::move(v)) {
        normalize();
    }

    void normalize() {
        auto first_nonzero = value.find_first_not_of('0');
        if(first_nonzero == string::npos) {
            value = "0";
        } else if(first_nonzero > 0) {
            value.erase(0, first_nonzero);
        }
    }
};

inline BigCount operator+(const BigCount& left, const BigCount& right) {
    string a = left.value;
    string b = right.value;
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());

    string result;
    int carry = 0;
    const size_t n = max(a.size(), b.size());
    result.reserve(n + 1);
    for(size_t i = 0; i < n; ++i) {
        int digit = carry;
        if(i < a.size()) digit += a[i] - '0';
        if(i < b.size()) digit += b[i] - '0';
        result.push_back(char('0' + (digit % 10)));
        carry = digit / 10;
    }
    if(carry > 0) result.push_back(char('0' + carry));
    reverse(result.begin(), result.end());
    return BigCount(result);
}

inline BigCount operator*(unsigned int factor, const BigCount& count) {
    if(factor == 0 || count.value == "0") return BigCount(0);

    string result;
    int carry = 0;
    result.reserve(count.value.size() + 2);
    for(auto it = count.value.rbegin(); it != count.value.rend(); ++it) {
        int digit = (*it - '0') * static_cast<int>(factor) + carry;
        result.push_back(char('0' + (digit % 10)));
        carry = digit / 10;
    }
    while(carry > 0) {
        result.push_back(char('0' + (carry % 10)));
        carry /= 10;
    }
    reverse(result.begin(), result.end());
    return BigCount(result);
}

inline BigCount operator*(const BigCount& count, unsigned int factor) {
    return factor * count;
}

inline ostream& operator<<(ostream& os, const BigCount& count) {
    os << count.value;
    return os;
}

// Algorithm Parameters
struct AlgoLoopParams{
    double s_step;
    string s_string;
    int K;
    double infid_target;
    string str_infid_target;
    double schedule_factor = 1.0;
    vector<double> s_candidates = {};
    bool refine_s = false;
};

struct AlgoOptParams{
    double s_min;
    double s_bin;
    double s_N;
    int K_max;
    double infid_target;
    string str_infid_target;
};

struct AlgoOptParamsVec{
    double s_min;
    double s_bin;
    double s_N;
    int K_max;
    vector <double> infid_target;
    vector <string> str_infid_target;
};

struct AlgoInfidParams{
    double s_step;
    string s_string;
    int K_max;
    double infid_min;
    int infid_N;
};


//Helper for Hamiltonians
struct PauliTerm{
    cplx coeff;
    char left;
    char right;
};


// Data Containers
struct Row{
    int k;
    double energy;
    double fidelity;
};

struct RowOpt {
    double s_step;
    int k;
    double infidelity;
};

struct RowInfid {
    double infid_step;
    int k;
    double infidelity;
};

struct DBQITEUnitaryCounts {
    BigCount U0 = 0;
    BigCount U0_dag = 0;
    BigCount A = 0;
    BigCount A_dag = 0;
    BigCount R0 = 0;
    BigCount R0_dag = 0;

    BigCount total() const {
        return U0 + U0_dag + A + A_dag + R0 + R0_dag;
    }

    DBQITEUnitaryCounts daggered() const {
        DBQITEUnitaryCounts result;
        result.U0 = U0_dag;
        result.U0_dag = U0;
        result.A = A_dag;
        result.A_dag = A;
        result.R0 = R0_dag;
        result.R0_dag = R0;
        return result;
    }

    void advanceOneDBQITEStep() {
        const DBQITEUnitaryCounts old = *this;
        const DBQITEUnitaryCounts old_dag = old.daggered();

        U0 = 2 * old.U0 + old_dag.U0;
        U0_dag = 2 * old.U0_dag + old_dag.U0_dag;
        A = 2 * old.A + old_dag.A + 1;
        A_dag = 2 * old.A_dag + old_dag.A_dag + 1;
        R0 = 2 * old.R0 + old_dag.R0 + 1;
        R0_dag = 2 * old.R0_dag + old_dag.R0_dag;
    }
};

struct DBQITEUnitaryTraceRow {
    int k;
    double infidelity;
    DBQITEUnitaryCounts counts;
};

struct DBQITEUnitaryOptTraceRow {
    double s_step;
    int k;
    double infidelity;
    DBQITEUnitaryCounts counts;
};

struct EnergyAnalysis {
    double E0;
    double E1;
    double gap;
    double var0 = std::numeric_limits<double>::quiet_NaN();
    double var1 = std::numeric_limits<double>::quiet_NaN();
};



#endif //MY_PROJECT_STRUCTFILE_H
