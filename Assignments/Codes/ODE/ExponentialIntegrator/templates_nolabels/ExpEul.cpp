//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <vector>
#include <iomanip>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

//! \brief Function $\phi$ used in the Exponential Euler
//! single step method for an autonomous ODE.
MatrixXd phim(const MatrixXd & Z) {
    int n = Z.cols();
    assert( n == Z.rows() && "Matrix must be square.");
    MatrixXd C(2*n,2*n);
    C << Z, MatrixXd::Identity(n,n),
         MatrixXd::Zero(n,2*n);
    return C.exp().block(0, n, n, n);
}

//! \brief Calculate a single step of the exponential Euler method.
//! \tparam Function function object for r.h.s. function
//! \tparam Jacobian function object for Jacobian of r.h.s.
//! \param[in] y0 The initial state
//! \param[in] f The r.h.s function $f$
//! \param[in] df The Jacobian of $f$
//! \param[in] h The stepsize of the method.
//! \return A single step of the Exponential Euler method
template <class Function, class Jacobian>
VectorXd ExpEulStep(const VectorXd & y0,
                    const Function& f, const Jacobian & df,
                    double h) {
    // TODO: Implement one step of the Exponential Euler method
}

int main() {
    // TODO: Test the exponential Euler method with the logistic ODE
    // and determine the approximated order of algebraic convergence.
}
