//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <iomanip>
#include <cmath>


int main() {

    // Lambda function, call with sinh(x)
    auto sinh_unstable = [] (double x) {
        double t = std::exp(x);
        return .5 * (t - 1./t);
    };

    // TODO: compute relative error for $10^{-k}, k  = 1,\dots,10$
}
