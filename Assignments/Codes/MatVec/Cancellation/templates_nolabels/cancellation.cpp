//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

int main() {
    // TODO: Compute approximation of the derivative of sin(x)
    // Print the error of each computation

    // Plot
    mgl::Figure fig;
    fig.setlog(true, true);
    fig.legend();
    fig.title("Error of approximation of f'(x_0)");
    fig.xlabel("h");
    fig.ylabel("| f'(x_0) - g_i(x_0, h) |");
    // TODO: plot errors and save figure
}
