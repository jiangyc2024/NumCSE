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
 * \param f The r.h.s function for the ODE.
 * \param T Final time.
 * \param y0 Initial data.
 * \param A Butcher matrix $A$.
 * \param b Buthcer vector $b$.
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
    // Construct data for Butcher schemes
    MatrixXd A1 = MatrixXd::Zero(1,1);
    VectorXd b1(1);
    b1 << 1;

    MatrixXd A2 = MatrixXd::Zero(2,2);
    A2(1,0) = 1;
    VectorXd b2(2);
    b2 << .5, .5;

    MatrixXd A3 = MatrixXd::Zero(3,3);
    A3(1,0) = .5;
    A3(2,0) = -1;
    A3(2,1) = 2;
    VectorXd b3(3);
    b3 << 1./6, 2./3, 1./6;

    MatrixXd A4 = MatrixXd::Zero(4,4);
    A4(1,0) = .5;
    A4(2,1) = .5;
    A4(3,2) = 1;
    VectorXd b4(4);
    b4 << 1./6, 1./3, 1./3, 1./6;

    // TODO: call order for all combinations of ODE ans RK method
}
