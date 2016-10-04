//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <iomanip>
#include <cmath>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* \brief Multiply A*x
 * Build matrix with block structure
 * \param[in] d1 n-dim vector
 * \param[in] d2 n-dim vector
 * \param[in] c n-dim vector
 * \param[in] x A vector with dim 2n
 * \return A*x
 */
VectorXd multA(const VectorXd & d1, const VectorXd & d2,
               const VectorXd & c, const VectorXd & x) {
    int n = d1.size();
    assert(n == d2.size()
           && n == c.size()
           && 2*n == x.size()
           && "Size mismatch!");

    ArrayXd y(2*n);

    // TODO: efficiently compute A*x


    return y;
}

/* \brief Solve A*x = b for x
 * Use partial pivotisation and exploit matrix structure
 * \param[in] d1_ n-dim vector
 * \param[in] d2_ n-dim vector
 * \param[in] c_ n-dim vector
 * \param[in] b_ A r.h.s vector with dim 2n
 * \return x s.t. A*x = b
 */
VectorXd solveA(const VectorXd & d1_, const VectorXd & d2_,
                const VectorXd & c_, const VectorXd & b_) {
    int n = d1_.size();
    assert(n == d2_.size()
           && n == c_.size()
           && 2*n == b_.size()
           && "Size mismatch!");

    // TODO: efficiently solve system using partial pivotisation

    return VectorXd(2*n);

}

int main() {
    VectorXd d1, d2, c, b;
    int n = 20;

    // Helper function to setup vectors with prescribed size
    auto buildvectors = [&d1, &d2, &c, &b] (int n) {
        d1 = VectorXd::LinSpaced(n, 1, n);
        d2 = -d1;
        c = VectorXd::Ones(n);
        b.resize(2*n);
        b << d1, d1;
    };

    //// Testing correcntess of multA
    std::cout << "--> multA test" << std::endl;
    buildvectors(20);
    MatrixXd A(2*n,2*n);
    A << static_cast<MatrixXd>(d1.asDiagonal()), static_cast<MatrixXd>(c.asDiagonal()),
         static_cast<MatrixXd>(c.asDiagonal()),  static_cast<MatrixXd>(d2.asDiagonal());
    VectorXd ye = A*b;
    VectorXd yo = multA(d1, d2, c, b);
    double err = (ye - yo).norm();

    std::cout << "Error:    " << err
              << std::endl;

    //// Testing correctness of solveA
    std::cout << "--> solveA test" << std::endl;
    ye = A.partialPivLu().solve(b);
    yo = solveA(d1, d2, c, b);
    err = (ye - yo).norm();
    std::cout << "Error:    " << err
              << std::endl;
    double res = (A*yo - b).norm();
    std::cout << "Residual: " << res
              << std::endl;

    //// Runtime measurements
    int repeats = 10;
    // TODO: implement runtime measurements
}
