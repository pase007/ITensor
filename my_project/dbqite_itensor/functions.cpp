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

// Build Itensor for general dimensions
//Itensor opFromDims(){}


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


// exp(i theta H)
ITensor exp_i_theta_H(ITensor const& H, double theta){
    ITensor A = expHermitian(H, cplx(0.0, theta));
    A.swapPrime(1,0);
    return A;
}


// ---------------------------------
// ---  Build Itensor Operations ---
// ---------------------------------

//Apply gate and take care of index
ITensor applyGate(ITensor const& G, ITensor psi){
    ITensor out = G * psi;   // out has prime level 1 on site index
    out = removeTrivial(out);
    out.noPrime();
    double nrm = norm(out);
    if(!std::isfinite(nrm) || nrm == 0.0) itensor::error("overflow/underflow in applyGate");
    out /= nrm;
    return out;
}


// Multiplicate two Matricies with right indexing
ITensor composeGate(ITensor L, ITensor R, Index const& s){
    auto sp  = prime(s);
    auto mid = prime(s,2);

    L = replaceInds(L, {s}, {mid});   // L(sp,mid)
    R = replaceInds(R, {sp}, {mid});  // R(mid,s)

    auto C = L * R;                   // (sp,s)
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


// Build projector from a state
ITensor projectorOnState(ITensor const& psi, Index const& s) {
    Index sp = prime(s);
    ITensor psiDag = dag(psi);
    psiDag = replaceInds(psiDag, {s}, {sp});
    ITensor P = psiDag * psi;
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


ITensor reunitarize_polar_svd(ITensor const& U, Index const& s){
    cout << "Step 1 - ";
    cout << "Step 2 - ";
    ITensor W, S, V;
    cout << "order(U) = " << order(U) << "\n";
    cout << inds(U);
    cout << "Step 3 - ";
    svd(U, W, S, V);
    cout << "Step 4 - ";
    auto Q = W * adjointGate(V);
    cout << "Step 5 - ";
    return canonGate(Q, s);
}
// -------------------------------------------------
// --- Build functions for measuring Observables ---
// -------------------------------------------------

// Compute <psi|Op|psi>
double expectation(ITensor const& psi, ITensor const& Op, Index const& s){
    Index sp = prime(s);
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


