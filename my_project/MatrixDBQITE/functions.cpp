//
// Created by Pascal Knoll on 21.02.26.
//
#include "functions.h"
#include <unsupported/Eigen/MatrixFunctions>

// exp(i theta |0><0|) = diagonal phase on |0>
Mat phase_on_zero(int dim, double theta){
    Mat R = Mat::Identity(dim, dim);
    R(0,0) = exp(cplx(0.0, theta));
    return R;
}

// exp(i theta H) = evolution
Mat evolution(double theta, const Mat& H) {
    Mat A  = (cplx(0.0,1.0) * theta * H).exp();
    return A;
}

// One DB-QITE recursion step
Mat dbqite_step(const Mat& H, const Mat& U, double s){
    double theta = sqrt(s);
    Mat A = evolution(theta, H);
    Mat Ad = A.adjoint();
    Mat R0 = phase_on_zero(H.rows(), theta);
    return A * U * R0 * U.adjoint() * Ad * U;
}

// Energy expectation value <psi|H|psi>
double energy(const Mat& H, const Vec& psi){
    cplx v = psi.adjoint() * H * psi;
    return real(v);
}

// Fidelity |<phi|psi>|^2
double fidelity(const Vec& target, const Vec& psi){
    cplx ov = target.adjoint() * psi;
    return norm(ov);
}


double unitarity_defect(const Mat& U){
    const int n = (int)U.rows();
    Mat D = (U.adjoint() * U).eval();
    D.diagonal().array() -= cplx(1.0, 0.0);
    return D.norm(); // Frobenius norm
}


// Returns a "reunitarized" matrix: U * (U†U)^(-1/2)
Mat reunitarize_polar(const Mat& U, double eps){
    const int n = (int)U.rows();
    if(U.cols() != n) {
        throw std::runtime_error("reunitarize_polar: U must be square");
    }

    // G = U† U  (Hermitian positive-definite if U is close to unitary)
    Mat G = (U.adjoint() * U).eval();

    // SelfAdjointEigenSolver works for Hermitian matrices
    Eigen::SelfAdjointEigenSolver<Mat> es(G);
    if(es.info() != Eigen::Success) {
        throw std::runtime_error("reunitarize_polar: eigensolver failed");
    }

    // Eigenvalues should be ~1. Clamp small/negative values due to numerical drift.
    Eigen::Matrix<double, Eigen::Dynamic, 1> evals = es.eigenvalues().real();

    Eigen::Matrix<double, Eigen::Dynamic, 1> inv_sqrt(n);
    for(int i = 0; i < n; ++i) {
        double lam = evals(i);
        if(lam < eps) lam = eps;                 // clamp
        inv_sqrt(i) = 1.0 / std::sqrt(lam);
    }

    // G^{-1/2} = V * diag(inv_sqrt) * V†
    Mat Vinv2 = es.eigenvectors()
              * inv_sqrt.asDiagonal()
              * es.eigenvectors().adjoint();

    // U <- U * G^{-1/2}
    return (U * Vinv2).eval();
}