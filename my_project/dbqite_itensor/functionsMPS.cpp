//
// Created by Pascal Knoll on 05.05.26.
//
#include "functionsMPS.h"

// ---------------------------------------- Calculations with MPOs ------------------------------------
MPO composeMPO(MPO const& A, MPO const& B){
    auto args = Args("Cutoff",1E-12,"MaxDim",20, "Verbose", true);
    MPO C;
    // prime only physical site indices of B
    auto Bp = prime(B, "Site");
    nmultMPO(A, Bp, C, args);

    // bring result back to standard MPO form: s', s
    C.mapPrime(2,1,"Site");

    return C;
}

MPO dagMPO(MPO const& U){
    MPO Ud(length(U));

    for(int i = 1; i <= length(U); ++i){
        auto T = dag(U.A(i));

        // Restore standard MPO convention: Site' , Site
        T.swapPrime(0,1,"Site");
        Ud.setA(i,T);
    }
    return Ud;
}

MPO buildIdentityMPO(SiteSet const& sitesVec, int NrSites){
    MPO I(NrSites);

    std::vector<Index> links;
    links.reserve(NrSites-1);

    for(int b = 1; b < NrSites; ++b)
        links.emplace_back(1, "Link,Id,l=" + std::to_string(b));

    if(NrSites == 1){
        I.setA(1, idGate(sitesVec(1)));
        return I;
    }

    {
        auto W = idGate(sitesVec(1));
        W *= setElt(links.at(0)(1));
        I.setA(1,W);
    }

    for(int i = 2; i <= NrSites-1; ++i){
        auto W = idGate(sitesVec(i));
        W *= setElt(links.at(i-2)(1));
        W *= setElt(links.at(i-1)(1));
        I.setA(i,W);
    }

    {
        auto W = idGate(sitesVec(NrSites));
        W *= setElt(links.back()(1));
        I.setA(NrSites,W);
    }
    return I;
}

MPO addMPO(MPO const& A, MPO const& B, Cplx a, Cplx b){
    int N = length(A);
    MPO C(N);

    if(length(B) != N) Error("addMPO: MPO lengths do not match");

    // New combined links
    std::vector<Index> cLinks;
    cLinks.reserve(N-1);

    for(int i = 1; i < N; ++i){
        auto la = linkIndex(A,i);
        auto lb = linkIndex(B,i);

        cLinks.emplace_back(dim(la) + dim(lb), "Link,Add,l=" + std::to_string(i));
    }

    for(int i = 1; i <= N; ++i){
        auto Ai = A.A(i);
        auto Bi = B.A(i);

        ITensor Ci;

        if(N == 1){
            Ci = a*Ai + b*Bi;
            C.setA(1,Ci);
            return C;
        }

        if(i == 1){
            auto la = linkIndex(A,1);
            auto lb = linkIndex(B,1);
            auto lc = cLinks.at(0);

            Ci = ITensor(inds(Ai)); // placeholder not ideal
            Ci = ITensor();

            for(int n = 1; n <= dim(la); ++n){
                auto part = Ai * setElt(la(n));
                part *= setElt(lc(n));
                Ci += a * part;
            }

            for(int n = 1; n <= dim(lb); ++n){
                auto part = Bi * setElt(lb(n));
                part *= setElt(lc(dim(la)+n));
                Ci += b * part;
            }
        }else if(i == N){
            auto la = linkIndex(A,N-1);
            auto lb = linkIndex(B,N-1);
            auto lc = cLinks.at(N-2);

            Ci = ITensor();

            for(int n = 1; n <= dim(la); ++n){
                auto part = Ai * setElt(la(n));
                part *= setElt(lc(n));
                Ci += part;
            }

            for(int n = 1; n <= dim(lb); ++n){
                auto part = Bi * setElt(lb(n));
                part *= setElt(lc(dim(la)+n));
                Ci += part;
            }
        }else{
            auto laL = linkIndex(A,i-1);
            auto laR = linkIndex(A,i);
            auto lbL = linkIndex(B,i-1);
            auto lbR = linkIndex(B,i);

            auto lcL = cLinks.at(i-2);
            auto lcR = cLinks.at(i-1);

            Ci = ITensor();

            // A block
            for(int l = 1; l <= dim(laL); ++l)
            for(int r = 1; r <= dim(laR); ++r){
                auto part = Ai * setElt(laL(l)) * setElt(laR(r));
                part *= setElt(lcL(l)) * setElt(lcR(r));
                Ci += part;
            }

            // B block
            int offL = dim(laL);
            int offR = dim(laR);

            for(int l = 1; l <= dim(lbL); ++l)
            for(int r = 1; r <= dim(lbR); ++r){
                auto part = Bi * setElt(lbL(l)) * setElt(lbR(r));
                part *= setElt(lcL(offL+l)) * setElt(lcR(offR+r));
                Ci += part;
            }
        }

        C.setA(i,Ci);
    }

    return C;
}



// ------------------------------------- Unitary defect calculation -----------------------------------
double unitarityDefectOnState(MPO const& U, MPS const& psi){
    auto args = Args("Method","DensityMatrix", "Cutoff",1E-12, "MaxDim",200);
    MPO Udag = dagMPO(U);

    MPS phi = applyMPO(U, psi, args);
    MPS back = applyMPO(Udag, phi, args);

    double nPsi  = real(innerC(psi, psi));
    double nBack = real(innerC(back, back));
    double ov    = real(innerC(psi, back));

    double d2 = nBack + nPsi - 2.0 * ov;

    if(d2 < 0.0 && std::abs(d2) < 1E-12) d2 = 0.0;

    return std::sqrt(d2);
}

MPO reunitarizeMPO(MPO const& U, SiteSet const& sitesVec, int NrSites){
    MPO Ud  = dagMPO(U);
    MPO UdU = composeMPO(Ud,U);

    MPO I = buildIdentityMPO(sitesVec,NrSites);

    MPO X = addMPO(I, UdU, Cplx(3.0,0.0), Cplx(-1.0,0.0));

    MPO Unew = composeMPO(U, X);

    // multiply first tensor by 1/2
    auto T = Unew.A(1);
    T *= 0.5;
    Unew.setA(1,T);

    return Unew;
}

void testTwoSiteGate(SiteSet const& sitesVec, vector<OneSigmaSite> const& models, int i, double theta, double g){
    auto h = bondHamiltonian(models, i, g);
    auto G = expHermitian(h, Cplx(0.0, theta));

    auto s1 = sitesVec(i);
    auto s2 = sitesVec(i+1);

    // -------------------------------
    // Hermiticity of h
    // -------------------------------
    auto hdag = dag(h);
    hdag.swapPrime(0,1,"Site");

    double herm_def = norm(h - hdag);

    // -------------------------------
    // Unitarity of G
    // Need Gdag * G as operator product.
    // -------------------------------
    auto Gdag = dag(G);
    Gdag.swapPrime(0,1,"Site");

    auto Gp = prime(G, "Site");      // Site' -> Site'', Site -> Site'
    auto GG = Gdag * Gp;             // contracts middle Site' legs
    GG.mapPrime(2,1,"Site");         // restore Site'' -> Site'

    auto I1 = idGate(s1);
    auto I2 = idGate(s2);
    auto Id = I1 * I2;

    double unit_def = norm(GG - Id);

    cout << "\n--- Two-site gate test ---\n";
    cout << "bond i = " << i << "\n";
    cout << "theta  = " << theta << "\n";
    cout << "Hermiticity defect ||h-h†|| = " << herm_def << "\n";
    cout << "Unitarity defect ||G†G-I|| = " << unit_def << "\n";
}



// ---------------------------- Helper functions for Prime handling and BondDim printing -------------------
void normalizeMPOPrimeLevels(MPO& K){
    for(int i = 1; i <= length(K); ++i){
        auto T = K.A(i);
        T.mapPrime(3,1,"Site");
        K.setA(i,T);
    }
}

int maxBondDim(MPO const& M){
    int m = 1;
    for(int b = 1; b < length(M); ++b){
        auto l = linkIndex(M,b);
        m = std::max(m, static_cast<int>(dim(l)));
    }
    return m;
}


// ---------------------------------- Helper fucntion to make MPS --------------------------------------
MPS buildProductStateHelper(SiteSet const& sitesVec, int NrSites, IndexSet const& links) {
    MPS p0 = MPS(NrSites);
    // Site 1: (right_bond, s)
    {
        auto s = sitesVec(1);
        ITensor A(links[0], s);
        A.set(links[0](1), s(1), 1.0);
        A.set(links[0](1), s(2), 0.0);
        A.set(links[0](1), s(3), 0.0);
        A.set(links[0](1), s(4), 0.0);
        p0.setA(1, A);
    }

    // Middle sites: (left_bond, right_bond, s)
    for(int i = 2; i < NrSites; ++i) {
        auto s = sitesVec(i);
        ITensor A(links[i-2], links[i-1], s);
        A.set(links[i-2](1), links[i-1](1), s(1), 1.0);
        A.set(links[i-2](1), links[i-1](1), s(2), 0.0);
        A.set(links[i-2](1), links[i-1](1), s(3), 0.0);
        A.set(links[i-2](1), links[i-1](1), s(4), 0.0);
        p0.setA(i, A);
    }

    // Last site: (left_bond, s)
    {
        auto s = sitesVec(NrSites);
        ITensor A(links[NrSites-2], s);
        A.set(links[NrSites-2](1), s(1), 1.0);
        A.set(links[NrSites-2](1), s(2), 0.0);
        A.set(links[NrSites-2](1), s(3), 0.0);
        A.set(links[NrSites-2](1), s(4), 0.0);
        p0.setA(NrSites, A);
    }
    return p0;
}



// ------------------------------- Helper function for Refelction GAte R0 ---------------------------------
// ---


// ---------------------------------- Helper function to make warmUp U0 --------------------------------------
// ---


// -------------------------------- Helper functions for Trotterization of A ------------------------------------
pair<ITensor,ITensor> gateToMPOBond(SiteSet const& sitesVec, ITensor G, int i, string tag){
    auto s1 = sitesVec(i);
    auto s2 = sitesVec(i+1);
    auto left_inds = IndexSet(prime(s1), s1);

    auto [U,S,V] = svd(G, left_inds, {"Cutoff",1E-14});

    auto W1 = U * S;
    auto W2 = V;
    return {W1, W2};
}

ITensor bondHamiltonian(vector<OneSigmaSite> const& models, int i, double g){
    int NrSites = models.size();

    auto const& mi = models.at(i-1);
    auto const& mj = models.at(i);

    double c = -3.0 / (4.0 * g * g);
    double wl = (i == 1) ? 1.0 : 0.5;
    double wr = (i == NrSites - 1) ? 1.0 : 0.5;

    ITensor h = wl * 0.5 * mi.H0() * mj.Id_s()
              + wr * 0.5 * mi.Id_s() * mj.H0()
              + c * ( mi.j1() * mj.j1()
                    + mi.j2() * mj.j2()
                    + mi.j3() * mj.j3());
    return h;
}

ITensor twoSiteGate(vector<OneSigmaSite> const& models, int i, double theta, double g) {
    auto h = bondHamiltonian(models, i, g);
    auto G = expHermitian(h, Cplx(0.0, theta));
    return G;
}

MPS applyTrotterLayerToMPS(MPS psi, vector<OneSigmaSite> const& models, double theta, bool oddLayer, int NrSites, double g, Args const& args){
    int start = oddLayer ? 1 : 2;

    for(int b = start; b < NrSites; b += 2){
        auto G = twoSiteGate(models, b, theta, g);

        psi.position(b);

        auto wf = psi.A(b) * psi.A(b+1);

        // Apply two-site gate
        wf = G * wf;

        // Remove primes on physical indices after gate application
        wf.noPrime("Site");

        // Split back into MPS tensors
        psi.svdBond(b, wf, Fromleft, args);
    }
    return psi;
}



// ---------------------------------- H as MPO helper functions --------------------------------------
ITensor makeFirstMPO_(Index const& wR, OneSigmaSite const& site, double coupling){
    auto Id = site.Id_s();
    auto H0 = site.H0();
    auto j1 = site.j1();
    auto j2 = site.j2();
    auto j3 = site.j3();

    ITensor W(wR);
    W += setElt(wR(5)) * Id;
    W += setElt(wR(1)) * H0;
    W += setElt(wR(2)) * (coupling * j1);
    W += setElt(wR(3)) * (coupling * j2);
    W += setElt(wR(4)) * (coupling * j3);
    return W;
}

ITensor makeLastMPO_(Index const& wL, OneSigmaSite const& site){
    auto Id = site.Id_s();
    auto H0 = site.H0();
    auto j1 = site.j1();
    auto j2 = site.j2();
    auto j3 = site.j3();

    ITensor W(wL);
    W += setElt(wL(5)) * H0;
    W += setElt(wL(1)) * Id;
    W += setElt(wL(2)) * j1;
    W += setElt(wL(3)) * j2;
    W += setElt(wL(4)) * j3;
    return W;
}

ITensor makeBulkMPO_(Index const& wL, Index const& wR, OneSigmaSite const& site, double coupling){
    auto Id = site.Id_s();
    auto H0 = site.H0();
    auto j1 = site.j1();
    auto j2 = site.j2();
    auto j3 = site.j3();

    // Build Matrix containing Energy terms and interaction
    ITensor W(wL, wR);
    W += setElt(wL(5)) * setElt(wR(5)) * Id;
    W += setElt(wL(1)) * setElt(wR(1)) * Id;
    W += setElt(wL(2)) * setElt(wR(1)) * j1;
    W += setElt(wL(3)) * setElt(wR(1)) * j2;
    W += setElt(wL(4)) * setElt(wR(1)) * j3;
    W += setElt(wL(5)) * setElt(wR(1)) * H0;
    W += setElt(wL(5)) * setElt(wR(2)) * (coupling * j1);
    W += setElt(wL(5)) * setElt(wR(3)) * (coupling * j2);
    W += setElt(wL(5)) * setElt(wR(4)) * (coupling * j3);
    return W;
}
