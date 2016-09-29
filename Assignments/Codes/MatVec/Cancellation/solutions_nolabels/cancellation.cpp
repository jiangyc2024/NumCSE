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


using namespace Eigen;

int main() {
    // Array of values of $h$
    ArrayXd h = ArrayXd::LinSpaced(21, -20, 0.)
        .unaryExpr([] (double i) {
            return std::pow(10., i);
        });
    // Dummy array where to evaluate the derivative (1.2)
    ArrayXd x = ArrayXd::Constant(h.size(), 1.2);

    // Derivative
    ArrayXd g1 = (sin(x +h) - sin(x)) / h; // naive
    ArrayXd g2 = 2 * cos(x+0.5*h) * sin(0.5 * h) / h; // better
    ArrayXd ex = cos(x); // exact

    //// Print error

    // Table header
    std::cout << std::setw(15) << "h"
              << std::setw(15) << "exact"
              << std::setw(15) << "cancellation"
              << std::setw(15) << "error"
              << std::setw(15) << "improved"
              << std::setw(15) << "error" << std::endl;
    for(unsigned int i = 0; i < h.size(); ++i) {
        // Table entries
        std::cout << std::setw(15) << h(i)
                  << std::setw(15) << ex(i)
                  << std::setw(15) << g1(i)
                  << std::setw(15) << std::abs(g1(i) - ex(i))
                  << std::setw(15) << g2(i)
                  << std::setw(15) << std::abs(g2(i) - ex(i)) << std::endl;
    }

}
