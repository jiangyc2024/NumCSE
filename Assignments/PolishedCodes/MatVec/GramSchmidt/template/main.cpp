////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

#include "gramschmidt.hpp"

int main() {
  // Orthonormality test
  std::cout << "Gram-Schmidt test of subproblem (c) returns for n = 2: "
            << testGramSchmidt(2) << std::endl;
  std::cout << "Gram-Schmidt test of subproblem (c) return for n = 10: "
            << testGramSchmidt(10) << std::endl;
}
