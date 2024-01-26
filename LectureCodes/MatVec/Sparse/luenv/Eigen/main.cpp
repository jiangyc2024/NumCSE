///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "luenv.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>

using Eigen::MatrixXd;
using Eigen::Triplet;
using Eigen::SparseMatrix;
using envelope::luenv;

int main() {
  const int n = 7;
  typedef Triplet<double> triplet_t;
  std::vector<triplet_t> triplets;
  triplets.reserve(n);

  // example in lecture document ex:env
  for (int i = 0; i < n; ++i) {
    triplets.emplace_back(i, i, 1.0);  // diagonal
  }

  triplets.emplace_back(0, 2, 1.0);
  triplets.emplace_back(2, 0, 1.0);

  triplets.emplace_back(1, 4, 1.0);
  triplets.emplace_back(4, 1, 1.0);

  triplets.emplace_back(2, 6, 1.0);
  triplets.emplace_back(6, 2, 1.0);

  triplets.emplace_back(3, 4, 1.0);
  triplets.emplace_back(4, 3, 1.0);

  triplets.emplace_back(3, 6, 1.0);
  triplets.emplace_back(6, 3, 1.0);

  triplets.emplace_back(4, 5, 1.0);
  triplets.emplace_back(5, 4, 1.0);

  SparseMatrix<double> A(n, n);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();

  const MatrixXd B(A);
  std::cout << "A=\n" << B << std::endl;

  MatrixXd L = MatrixXd::Zero(n, n);
  MatrixXd U = MatrixXd::Zero(n, n);
  luenv(B, L, U);
  std::cout << "L=\n" << L << std::endl;
  std::cout << "U=\n" << U << std::endl;
  std::cout << "L*U=\n" << L*U << std::endl;
  return 0;
}
