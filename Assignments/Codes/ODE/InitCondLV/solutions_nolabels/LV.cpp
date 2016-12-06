//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "ode45.hpp"

using namespace Eigen;

/*!
 * \brief Compute the maps Phi and W at time T.
 * Use initial data given by u0 and v0.
 * \param[in] u0 First component.
 * \param[in] v0 Second component.
 * \param[in] T Final time.
 */
std::pair<Vector2d, Matrix2d> PhiAndW(double u0,
                                      double v0,
                                      double T) {
    std::pair<Vector2d, Matrix2d> PaW;
    auto f = [] (const VectorXd & w) {
        Eigen::VectorXd temp(6);
        temp(0) = (2. - w(1))*w(0);
        temp(1) = (w(0) - 1.)*w(1);
        temp(2) = (2. - w(1))*w(2) - w(0)*w(3);
        temp(3) = w(1)*w(2) + (w(0) - 1.)*w(3);
        temp(4) = (2. - w(1))*w(4) - w(0)*w(5);
        temp(5) = w(1)*w(4) + (w(0) - 1.)*w(5);
        return temp;
    };

    Eigen::VectorXd w0(6);
    w0 << u0, v0, 1., 0, 0, 1.;

    // Construct ode solver with r.h.s
    ode45<Eigen::VectorXd> O(f);
    // Set options
    O.options.rtol = 1e-14;
    O.options.atol = 1e-12;
    // Solve ODE
    auto sol = O.solve(w0, T);
    // Extract needed component
    VectorXd wT = sol.back().first;

    PaW.first  << wT(0), wT(1);
    PaW.second << wT(2), wT(4), wT(3), wT(5);
    return PaW;
}

int main(){
    Vector2d y;
    y << 3, 2;
    double T = 5;

    std::pair<Vector2d,Matrix2d> PaW = PhiAndW(y(0), y(1), T);
    Vector2d F = PaW.first - y;
    Matrix2d DF;

    while (F.norm() > 1e-5) {
        PaW = PhiAndW(y(0), y(1), T);
        F = PaW.first - y;
        DF = PaW.second - MatrixXd::Identity(2,2);
        y = y - DF.lu().solve(F);
    }

    std::cout << "The obtained initial condition is: "
         << std::endl << y << std::endl;
    PaW = PhiAndW(y(0), y(1), 100);

    std::cout << "y(100) = " << std::endl
              << PaW.first << std::endl;
}
