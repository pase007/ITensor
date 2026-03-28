//
// Created by Pascal Knoll on 21.02.26.
//
#include "functions.h"

// ---------------------------------
// --- Build Itensor from Matrix ---
// ---------------------------------

// Build 2x2 Itensor from some Array
ITensor opFromArray(Index const& s, cplx a00, cplx a01, cplx a10, cplx a11){
    auto sp = prime(s);
    ITensor U(sp, s);
    U.set(sp(1), s(1), a00);
    U.set(sp(1), s(2), a01);
    U.set(sp(2), s(1), a10);
    U.set(sp(2), s(2), a11);
    return U;
}


// ---------------------------------
// ---- Build ITensor for Gates ----
// ---------------------------------

// 2D Identity Gate
ITensor idGate(Index const& s){
    Index sp = prime(s);
    int d = dim(s);
    ITensor I(sp, s);
    for (int i=1; i<=d; i++) {
        I.set(sp(i), s(i), cplx(1.0, 0.0));
    }
    return I;
}

// N qubit Id generalization
ITensor idGateN(IndexSet const& sites) {
    ITensor I(1.0);
    for(auto const& s : sites) {
        I *= idGate(s);
    }
    return I;
}


// R0(theta) = diag(e^{i theta}, 1) Phase gate
ITensor phaseOnZero(Index const& s, double theta){
    Index sp = prime(s);
    cplx p = exp(cplx(0.0, theta));
    ITensor R0 = idGate(s);
    R0.set(sp(1), s(1), p);
    return R0;
}


// R_psi(theta) Phase gate
ITensor phaseOnState(ITensor const& psi_ref, Index const& s, double theta){
    ITensor I = idGate(s);
    ITensor P = projectorOnState(psi_ref, s);
    return I + (exp(cplx(0.0, theta)) - cplx(1.0,0.0)) * P;
}

// R_psi(theta) Phase gate
ITensor phaseOnStateN(ITensor const& psi_ref, IndexSet const& svec, double theta){
    ITensor I = idGateN(svec);
    ITensor P = projectorOnStateN(psi_ref);
    return I + (exp(cplx(0.0, theta)) - cplx(1.0,0.0)) * P;
}


// exp( i H theta)
ITensor exp_i_theta_H(ITensor const& H, double theta){
    //for(auto const& I : inds(H)){
        //cout << I << " H prime=" << primeLevel(I) << "\n";
    //}
    return expHermitian(H, cplx(0.0, theta));
}


// ---------------------------------
// ---  Build Itensor Operations ---
// ---------------------------------

//Apply gate and take care of index
ITensor applyGate(ITensor const& G, ITensor psi){
    ITensor out = G * psi;   // out has prime level 1 on site index
    //out = removeTrivial(out);
    out.noPrime();
    double nrm = norm(out);
    if(!std::isfinite(nrm) || nrm == 0.0) itensor::error("overflow/underflow in applyGate");
    out /= nrm;
    return out;
}


// Multiplicate two Matricies with right indexing
ITensor composeGate(ITensor L, ITensor R, Index const& s){
    Index sp  = prime(s);
    Index mid = prime(s,2);

    L = replaceInds(L, {s}, {mid});   // L(sp,mid)
    R = replaceInds(R, {sp}, {mid});  // R(mid,s)

    ITensor C = L * R;                   // (sp,s)
    return C;
}

ITensor composeGateN(ITensor L, ITensor R, IndexSet const& sites){
    vector<Index> ket_inds;
    vector<Index> bra_inds;
    vector<Index> mid_inds;

    for(Index const& s : sites){
        ket_inds.push_back(s);
        bra_inds.push_back(prime(s));
        mid_inds.push_back(prime(s,2));
    }

    // L(bra, ket) -> L(bra, mid)
    L.replaceInds(ket_inds, mid_inds);

    // R(bra, ket) -> R(mid, ket)
    R.replaceInds(bra_inds, mid_inds);

    // contracts over all mid indices
    ITensor C = canonGateN(L * R, sites);
    //cout << "C = " << C << endl;
    return C;
}

ITensor adjointGate(ITensor const& U){
    IndexSet Ind = inds(U);
    Index a = Ind[0], b = Ind[1];
    ITensor Udag = permute(U, b, a);
    Udag = dag(Udag);
    Udag = swapPrime(Udag, 1, 0);
    return Udag;
}

ITensor adjointGateN(ITensor const& U){
    vector<Index> perm_inds;
    IndexSet is = inds(U);

    // swap each neighboring pair: (s',s,s',s,...) -> (s,s',s,s',...)
    for(int j = 0; j < length(is); j += 2){
        perm_inds.push_back(is[j+1]);
        perm_inds.push_back(is[j]);
    }

    ITensor Udag = permute(U, perm_inds);
    Udag = dag(Udag);
    Udag = swapPrime(Udag, 0, 1);

    return Udag;
}


//Evaluate Unitary defect of Tensor
double unitary_Defect(ITensor const& U, Index const& s) {
    ITensor Udag = adjointGate(U);
    ITensor UUdag = composeGate(U, Udag, s);
    //cout << "This should be identity already" << UUdag << endl;
    UUdag = canonGate(UUdag,s);
    ITensor Id = idGate(s);
    double def = norm(UUdag - Id);
    return def;
}

//Evaluate Unitary defect of Tensor
double unitary_DefectN(ITensor const& U, IndexSet const& sites) {
    ITensor Udag = adjointGateN(U);
    ITensor UUdag = composeGateN(U, Udag, sites);
    //cout << "This should be identity already" << UUdag << endl;
    ITensor Id = idGateN(sites);
    double def = norm(UUdag - Id);
    return def;
}


// Build projector from a state
ITensor projectorOnState(ITensor const& psi, Index const& s) {
    Index sp = prime(s);
    ITensor psiDag = dag(psi);
    psiDag = replaceInds(psiDag, {s}, {sp});
    ITensor P = psiDag * psi;
    return P;
}

ITensor projectorOnStateN(ITensor const& psi) {
    // collect ket indices
    vector<Index> ket_inds;
    vector<Index> bra_inds;
    for(Index const& i : inds(psi)){
        if(primeLevel(i) == 0){
            ket_inds.push_back(i);
            bra_inds.push_back(prime(i));
        }
    }

    // build <psi|
    ITensor psiDag = dag(psi);
    psiDag.replaceInds(ket_inds, bra_inds);

    // build |psi><psi|
    ITensor P = psi * psiDag;

    return P;
}


ITensor canonGate(ITensor G, Index const& s){
    Index sp = prime(s);
    if(order(G) != 2) itensor::error("canonGate: gate must be rank-2");
    IndexSet IG = inds(G);
    Index a = IG[0];
    Index b = IG[1];
    ITensor Out = replaceInds(G, {a, b}, {sp, s});
    return Out;
}

ITensor canonGateN(ITensor const& G, IndexSet const& sites){
    vector<Index> ord;
    for(Index const& s : sites){
        ord.push_back(prime(s));
        ord.push_back(s);
    }
    return permute(G, ord);
}

// Reunitarize matrix U
ITensor reunitarize_polar_gate(ITensor const& U, Index const& s, double eps){
    // Build G = U*Udag
    Index sp = prime(s);
    ITensor G = composeGate(adjointGate(U), U, s);
    if(!(hasIndex(G,prime(s)) && hasIndex(G,s))) itensor::error("G must have indices (prime(s), s)");
    if(!hasIndex(G, sp) || !hasIndex(G, s))
        itensor::error("reunitarize_polar_gate: G does not have indices (prime(s), s)");

    // Diagonalize Hermitian G:  G = V D V
    ITensor D, V;
    diagHermitian(G, V, D);
    V = prime(V);
    ITensor Vdag = adjointGate(V);

    // Build D^{-1/2} on the indices as D
    IndexSet indsD = inds(D);      // IndexSet
    Index d  = indsD[0];        // Index
    Index dp = indsD[1];        // Index

    ITensor Dinv2(d, dp);
    for(int n = 1; n <= dim(d); ++n){
        cplx lamC = eltC(D, d(n), dp(n));
        double lam = real(lamC);
        //cout << "lam = " << lam << "\t";
        if(lam < eps) {
            lam = eps;
            //cout << "negative iegenvalue lam = " << lam << "\t";
        }
        Dinv2.set(d(n), dp(n), 1.0/sqrt(lam));
    }

    // Q = U ∘ Ginv2
    ITensor Ginv2 = V * Dinv2 * Vdag;
    ITensor Q = composeGate(U, Ginv2, s);
    return Q;
}

ITensor reunitarize_polar_gateN(ITensor const& U, IndexSet const& sins, double eps){
    // Adjoint and G
    ITensor Udag = adjointGateN(U);
    ITensor G = composeGateN(Udag, U, sins);

    // Optional sanity check:
    for(Index const& s : sins){
        if(!hasIndex(G,s) || !hasIndex(G,prime(s)))
            itensor::error("reunitarize_polar_gateN: G missing expected indices");
    }

    // Diagonalize Hermitian G
    ITensor D, V;
    diagHermitian(G, V, D);
    ITensor Vdag = dag(V);
    V = prime(V);

    // Build D^{-1/2}
    IndexSet is = inds(D);
    Index d  = is[0];
    Index dp = is[1];

    ITensor Dinv2(d,dp);
    for(int n = 1; n <= dim(d); ++n){
        auto lamC = eltC(D,d(n),dp(n));
        double lam = real(lamC);
        if(lam < eps) lam = eps;
        Dinv2.set(d(n),dp(n),1.0/sqrt(lam));
    }

    // Reconstruct G^{-1/2}
    ITensor Ginv2 = V * Dinv2 * Vdag;
    ITensor Q = composeGateN(U, Ginv2, sins);
    return Q;
}



// -------------------------------------------------
// --- Build functions for measuring Observables ---
// -------------------------------------------------

// Compute <psi|Op|psi>
double expectation(ITensor const& psi, ITensor const& Op){
    ITensor OpPsi = Op * psi;
    ITensor psiDag = prime(dag(psi));
    double expect = real( (psiDag * OpPsi).eltC() );
    return expect;
}


// |<phi|psi>|^2
double fidelity(ITensor const& psi, ITensor const& phi){
    cplx fid = (dag(phi) * psi).eltC();
    return abs(fid)*abs(fid);
}

double fidelityToSubspace(ITensor const& psi, vector<ITensor> const& basis){
    double F = 0.0;

    for(auto const& phi : basis){
        cplx ov = (dag(phi) * psi).eltC();
        F += norm(ov);
    }
    return F;
}


