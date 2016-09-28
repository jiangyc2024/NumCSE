//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <Eigen/Dense>
#include <iostream>

#include "timer.h"

using namespace Eigen;

/* \brief Use efficient implementation A*x = bb
 * \param[in] R MatrixXd is nxn and upper triangular
 * \param[in] v VectorXd is nx1
 * \param[in] u VectorXd is nx1
 * \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
 * \param[out] x solution A*bb = x
 */
void solvelse(const MatrixXd & R,
              const VectorXd & v, const VectorXd & u,
              const VectorXd & bb, VectorXd & x) {
    // TODO: Implement a block-gauss elimination
}

//! \brief Use Eigen's LU-solver to solve Ax = y
//! \param[in] R MatrixXd is nxn and upper triangular
//! \param[in] v VectorXd is nx1
//! \param[in] u VectorXd is nx1
//! \param[in] bb vector is (n+1)x1 and is stacked (b, \beta)^T =: b
//! \param[out] x solution A*bb = x
void solvelse_lu(const MatrixXd & R,
                 const VectorXd & v, const VectorXd & u,
                 const VectorXd & bb, VectorXd & x) {
    // TODO: solve the system using LU decomposition
}

int main() {

    ///// CORRECTNESS TESTS

    // Bunch of random vectors/matrices
    int n = 9;

    // Random test vectors
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(n,n).triangularView<Eigen::Upper>();
    Eigen::VectorXd v = Eigen::VectorXd::Random(n);
    Eigen::VectorXd u = Eigen::VectorXd::Random(n);
    Eigen::VectorXd bb = Eigen::VectorXd::Random(n+1);
    Eigen::VectorXd xe, xo;

    // Check that answer is correct
    std::cout << "--> Check correctness." << std::endl;
    solvelse(R, v, u, bb, xo);
    std::cout << "Block Gauss:" << std::endl
              << xo << std::endl;
    // Compare with eigen partial pivot LU-solve
    solvelse_lu(R, v, u, bb, xe);
    std::cout << "Eigen LU:" << std::endl
              << xe << std::endl;
    std::cout << "Error:" << std::endl
              << (xe-xo).norm() << std::endl;

}
