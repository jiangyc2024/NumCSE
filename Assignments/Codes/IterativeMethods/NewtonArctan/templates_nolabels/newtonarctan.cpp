//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>


using namespace Eigen;

/* @brief
 * @param[in]
 * @param[out]
 */
double newton_arctan(double x0_) {

    double x0 = x0_;

// TODO:

    return x0;
}

int main() {
	// Initialization
    double x0_ = 2; // Initial guess

    // Netwon's method
    double x0 = newton_arctan(x0_);
    double x1 = x0 - std::atan(x0)*(1+x0*x0);
    double x2 = x1 - std::atan(x1)*(1+x1*x1);

}
