//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Compute $l_{min}$ from vector $d$
 * Naive implementation
 * @param[in] d An $n$-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
void rankoneinvit(const VectorXd & d,
                  const double & tol, double & lmin) {
    // Initialization
    VectorXd ev = d;
    lmin = 0;
    double lnew = d.cwiseAbs().minCoeff();

    while(std::abs(lnew-lmin)>tol*lmin) {
        Timer tm_slow;
        tm_slow.start();

        lmin = lnew;
        MatrixXd M = d.asDiagonal();

        M += ev*ev.transpose();
        ev = M.lu().solve(ev);

        ev.normalize();
        lnew = ev.transpose()*M*ev;

        tm_slow.stop();
    }

    lmin = lnew;
}

/* @brief Compute $l_{min}$ from vector $d$
 * Optimized implementation
 * @param[in] d An $n$-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
void rankoneinvit_fast(const VectorXd & d,
                       const double & tol, double & lmin)
{
    // Initialization
    VectorXd ev=d;
    lmin=0;
    double lnew=d.cwiseAbs().minCoeff();

    VectorXd dinv=(1/d.array()).matrix();
    while (std::abs(lnew-lmin)>tol*lmin) {
        Timer tm_fast;
        tm_fast.start();

        lmin = lnew;
        VectorXd ev0 = ev;

		// Here we solve the linear system
		// with the Sherman-Morrison-Woodbury formula
		// in the case of rank-1 perturbations.
        // This holds from $M = diag(d) + ev*ev^t$
        VectorXd Aib = dinv.cwiseProduct(ev);
        double temp = ev.transpose()*Aib;
        ev = Aib/(1+temp);

        ev.normalize();
		// Better than the corresponding naive implementation.
        // This holds from $M = diag(d) + ev*ev^t$, too
        lnew = ev.transpose()*d.cwiseProduct(ev)
            + pow(ev.transpose()*ev0,2);

        tm_fast.stop();
    }

    lmin=lnew;
}

int main() {
    // Initialization
    srand((unsigned int) time(0));
    double tol = 1e-3;
    double lmin;
    int n = 10;

    // Compute with both implementations
    VectorXd d = VectorXd::Random(n);
    std::cout << "Direct porting from MATLAB "
              <<"(naive implementation): "
              << std::endl;
    rankoneinvit(d,tol,lmin);
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "Fast implementation: " << std::endl;
    rankoneinvit_fast(d,tol,lmin);
    std::cout << "lmin = " << lmin << std::endl;

    // Compare runtimes of different
    // implementations of rankoneinvit
    std::cout << "*** Runtime comparison of two implementations"
              << std::endl;
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
