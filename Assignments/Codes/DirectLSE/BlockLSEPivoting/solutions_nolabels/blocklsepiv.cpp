//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
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

    auto & x1 = x.head(n).array();
    auto & x2 = x.tail(n).array();

    y << d1.array() * x1 + c.array() * x2,
         d2.array() * x2 + c.array() * x1;



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

    ArrayXd c1 = c_, c2 = c_, d1 = d1_, d2 = d2_, b = b_;

    double eps = std::numeric_limits<double>::epsilon();

    // For forward elimination + pivotisation we are
    // only required to loop the first half of the matrix
    // Loop over diagonal
    for(int k = 0; k < n; ++k) {
        // Check if need to pivot (i.e. swap two rows)
        // Luckily we only need to check two rows
        double maxk = std::max(std::abs(d1(k)), std::abs(c1(k)));
        double maxnk = std::max(std::abs(c2(k)), std::abs(d2(k)));
        if( std::abs( c1(k) ) / maxk // max relative pivot at row k
            >
            std::abs( d1(k) ) / maxnk // max relative pivot at rok k+n
            ) {
            // Matrix
            std::swap(d1(k), c2(k));
            std::swap(c1(k), d2(k));
            // R.h.s.
            std::swap(b(n), b(n+k));
        }

        // Check if matrix is almost singuloar
        double piv = d1(k);
        // Norm of the block from k,k to n-1,n-1
        double norm = std::abs(d1(k)) + std::abs(d2(k))
                    + std::abs(c1(k)) + std::abs(c2(k));
        if( piv < eps * norm ) {
            std::cout << "Warning: matrix nearly singular!" << std::endl;
        }

        // Multiplication facot:
        double fac = c2(k) / piv;

        // Actually perform substitution
        // Bottom Right poriton changes
        d2(k) -= c1(k) * fac;
        // R.h.s
        b(n+k) -= b(k) * fac;
    }

    // Now the system has the form:
    // | d1 | c  |   |   |   |   |
    // | 0  | d2 | * | x | = | b |
    // with d1, d2, c diagonal

    // Backward substitution

    // Lower potion
    b.tail(n) /= d2;
    // Upper portion
    b.head(n) = (b.head(n) - c1*b.tail(n)) / d1;

    return b;

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
    // Table header
    std::cout << std::setw(10) << "k"
              << std::scientific << std::setprecision(3)
              << std::setw(15) << "runtime [s]" << std::endl;
    // Loop from $2$ to $2^7$
    for(int k = 2; k <= (1<<7); k = k << 1) {
        buildvectors(k); // creates d1, d2, c, b

        // Repeat measurments
        Timer tm;
        for(int r = 0; r < repeats; ++r) {
            // Test run
            tm.start();
            yo = solveA(d1, d2, c, b);
            tm.stop();
        }

        // Print measurments
        std::cout << std::setw(10) << k
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << tm.duration() << std::endl;
    }
}
