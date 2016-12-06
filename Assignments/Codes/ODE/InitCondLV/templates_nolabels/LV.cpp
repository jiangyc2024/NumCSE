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
    // TODO:
    return PaW;
}

int main(){
    Vector2d y;
    y << 3, 2;
    double T = 5;

    // TODO: Apply the Newton method to find initial data
    // giving solutions with period equal to 5.
}
