//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

#include "strassen.hpp"

int main() {
    // Store seed for rng
    unsigned seed = (unsigned) time(0);
     // seed random number generator
    srand(seed);

    // Minimum number of repetitions
    unsigned int repetitions = 10;

    // TODO: time Strassen and Eigen matrix multiplication
    // repeat 10 times and output the minimum runtime
    // use 3 digits and scientific notation
}
