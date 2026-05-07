//
// Created by Pascal Knoll on 04.05.26.
//

#ifndef MY_PROJECT_TEST_DMRG_H
#define MY_PROJECT_TEST_DMRG_H
#include "itensor/all.h"
using namespace itensor;

void test_dmrg() {
    int N = 2;
    vector<Index> sites;
    for(int i = 0; i < N; ++i) {
        sites.push_back(Index(2, "S=" + to_string(i)));
    }

    // Build trivial MPS: |00>
    MPS psi(sites);
    for(int i = 1; i <= N; ++i) {
        auto s = sites[i-1];
        ITensor A(s);
        A.set(s(1), 1.0);
        psi.setA(i, A);
    }

    // Build trivial MPO: identity
    MPO H(N);
    for(int i = 1; i <= N; ++i) {
        auto s = sites[i-1];
        auto sp = prime(s);
        ITensor W(sp, s);
        for(int j = 1; j <= dim(s); ++j) {
            W.set(sp(j), s(j), 1.0);
        }
        H.setA(i, W);
    }

    // Try DMRG
    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 10;
    auto [E, psi_gs] = dmrg(H, psi, sweeps, {"Quiet", true});
    cout << "Energy: " << E << "\n";

}
#endif //MY_PROJECT_TEST_DMRG_H