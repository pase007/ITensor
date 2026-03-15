//
// Created by Pascal Knoll on 21.02.26.
// clang++ -std=c++17 main.cpp -I"$(brew --prefix eigen)/include/eigen3" -o dbqite
//
#include "functions.h"
#include "data.h"
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;
using cplx = complex<double>;
using Mat  = Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic>;
using Vec  = Eigen::Matrix<cplx, Eigen::Dynamic, 1>;

int main() {
    cout << fixed << setprecision(20);

    //Hamiltonian H
    Mat H = Mat(2,2);
    H << cplx(0.0,0.0), cplx(0.0,-1.0),
     cplx(0.0,1.0), cplx(0.0,0.0);

    //Target ground state/-energy of Hamiltonian
    Eigen::SelfAdjointEigenSolver<Mat> solver(H);
    Vec ground_state = solver.eigenvectors().col(0);
    double ground_energy = solver.eigenvalues()(0);

    cout << "Ground state" << ground_state << endl;
    cout << "Ground energy = " << ground_energy << "\n\n";

    //Set simulation params and start state
    double s = 0.1;
    int K = 50;

    Vec ket0(2);
    ket0 << 1.0,0.0;
    Mat U = Mat(2,2);
    U << cplx(1.0,0.0), cplx(0.0,0.0), cplx(0.0,0.0), cplx(1.0,0.0);
    cout << "U_def: " << unitarity_defect(U) << "\n";

    //Data container
    vector<int> ks;
    vector<double> Es, Fs;

    cout << "k\tEnergy\t\tFidelity\t\tUnitary defect\n";

    //recursion formula and evolution
    for (int k = 0; k < K; k++) {
        //state evolution
        Vec psi = U * ket0;
        psi /= psi.norm();
        double Ek = energy(H, psi);
        double Fk = fidelity(ground_state, psi);

        //show energy and fidelity of step
        ks.push_back(k);
        Es.push_back(Ek);
        Fs.push_back(Fk);
        cout << k << "\t"
                  << energy(H, psi) << "\t"
                  << fidelity(ground_state, psi)
                  << "\t";

        //do recursion step
        U = dbqite_step(H, U, s);
        double Udef = unitarity_defect(U);
        cout << "U_def: " << Udef << "\n";
        if (Udef > 1e-10){
            U = reunitarize_polar(U);
        }


    }
    const string csv = "matrix_data.csv";
    write_csv_kEF(csv, ks, Es, Fs);

    cout << "Done" << endl;

}