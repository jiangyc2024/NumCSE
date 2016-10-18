//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

using namespace Eigen;

/* \brief Use symmetric Gauss-Seidel iterations to solve the system Ax = b
 * \param[in] A system matrix to be decompsed (L + D + U)
 * \param[in] b r.h.s. vector
 * \param[in,out] x initial guess and last iterate (approximated solution)
 * \param[in] rtol relative tolerance for termination criteria
 */
void GSIt(const MatrixXd & A, const VectorXd & b,
          VectorXd & x, double rtol) {
    // TODO: Implement Gauss-Seidel iteration

    return;
}

int main(int, char**) {
    unsigned int n = 9;

    MatrixXd A = MatrixXd::Zero(n,n);
    for(unsigned int i = 0; i < n; ++i) {
        if(i > 0) A(i,i-1) = 2;
        A(i,i) = 3;
        if(i < n-1) A(i,i+1) = 1;
    }
    VectorXd b = VectorXd::Constant(n, 1);

    //// Subproblem b: test GSIt
    std::cout << "--> Test GSIt" << std::endl;

    // TODO: test the code GSIt using the given data
}
