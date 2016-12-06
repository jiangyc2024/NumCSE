//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>

#include "rkintegrator.hpp"

using namespace Eigen;

/*!
 * \brief errors Approximates the order of convergence of a scheme.
 * The scheme is a Runge Kutta scheme defined by A and b when applied
 * to the first order system y'=f(y), y(0)=y0.
 * The function ouputs error of the solutions at the time T.
 * \tparam Function Type for r.h.s function f.
 * \param f
 * \param T
 * \param y0
 * \param A
 * \param b
 */
template <class Function>
void errors(const Function &f, double T,
            const VectorXd &y0,
            const MatrixXd &A, const VectorXd &b) {

    RKIntegrator<VectorXd> rk(A, b);

    std::vector<double> error(15);
    // TODO: output error and order of the method
}

int main() {
    // TODO: output error and order of the method
}
