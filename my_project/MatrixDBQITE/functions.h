//
// Created by Pascal Knoll on 21.02.26.
//
#ifndef MY_PROJECT_FUNCTIONS_H
#define MY_PROJECT_FUNCTIONS_H

#include <Eigen/Dense>
#include <complex>
#include <limits>
using namespace Eigen;
using cplx = std::complex<double>;
using Mat  = Matrix<cplx, Dynamic, Dynamic>;
using Vec  = Matrix<cplx, Dynamic, 1>;
using namespace std;

//functions
Mat phase_on_zero(int dim, double theta);

Mat evolution(double s, const Mat& H);

Mat dbqite_step(const Mat& H, const Mat& U, double s);

double energy(const Mat& H, const Vec& psi);

double fidelity(const Vec& phi, const Vec& psi);

double unitarity_defect(const Mat& U);

Mat reunitarize_polar(const Mat& U, double eps = 1e-8L);

#endif //MY_PROJECT_FUNCTIONS_H