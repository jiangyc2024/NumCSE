//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Compute lmin from vector d, naive implementation
 * @param[in] d An n-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
void rankoneinvit(const VectorXd & d, const double & tol, double & lmin)
{
    VectorXd ev = d;
    lmin = 0;
    double lnew = d.cwiseAbs().minCoeff();

    while(abs(lnew-lmin)>tol*lmin) {
    }

    lmin = lnew;
}

/* @brief Compute lmin from vector d, optimized implementation
 * @param[in] d An n-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
void rankoneinvit_fast(const VectorXd & d, const double & tol, double & lmin)
{
    VectorXd ev=d;
    lmin=0;
    double lnew=d.cwiseAbs().minCoeff();

    VectorXd dinv=(1/d.array()).matrix();
    while (abs(lnew-lmin)>tol*lmin) {
    }

    lmin=lnew;
}

int main() {
    srand((unsigned int) time(0));
    double tol = 1e-3;
    double lmin;
    int n = 10;

    // Compute with both implementations
    VectorXd d = VectorXd::Random(n);
    std::cout << "Direct porting from MATLAB (naive implementation): " << std::endl;
    rankoneinvit(d,tol,lmin);
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "Fast implementation: " << std::endl;
    rankoneinvit_fast(d,tol,lmin);
    std::cout << "lmin = " << lmin << std::endl;

    // Compare runtimes of different implementations of rankoneinvit
    std::cout << "*** Runtime comparison of two implementations" << std::endl;
    unsigned int repeats = 3;
    Timer tm_slow, tm_fast;

    for(unsigned int p = 2; p <= 9; p++) {
        tm_slow.reset();
        tm_fast.reset();
        unsigned int n = pow(2,p);

        for(unsigned int r = 0; r < repeats; ++r) {
         // d = VectorXd::Random(n);
            d = VectorXd::LinSpaced(n,1,2);

	    tm_slow.start();
            rankoneinvit(d,tol,lmin);
	    tm_slow.stop();

	    tm_fast.start();
            rankoneinvit_fast(d,tol,lmin);
	    tm_fast.stop();
        }

        std::cout << "The slow method took: "
                  << tm_slow.min()
                  << " s for n = " << n << std::endl;
        std::cout << "The fast method took: "
                  << tm_fast.min()
                  << " s for n = " << n << std::endl;
    }
}
