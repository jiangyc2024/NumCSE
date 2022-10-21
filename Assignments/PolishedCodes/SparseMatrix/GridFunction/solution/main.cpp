//
// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
// Contributors: tille, jgacon, dcasati
// This file is part of the NumCSE repository.
//

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>

#include "gridfun.hpp"

// We use this as type of each index (will be int or similar)
using index_t = Eigen::Index;
// Use this to contain the "size" of our matrix
using shape_t = std::array<index_t, 2>;

int main() {
  // Size of the "grid" matrices $\mathbf{X}$ and $\mathbf{Y}$.
  constexpr index_t n = 16, m = 16;

  // The stencil matrix $\mathbf{S}$.
  Eigen::Matrix3d S;
  S << 0, 1, 0, 1, -4, 1, 0, 1, 0;
  // The grid matrices $\mathbf{X}$ and $\mathbf{Y}$.
  Eigen::MatrixXd X(n, m), Y;

  // A function of indices $(i,j)$ of the matrix $\mathbf{X}$.
  auto f = define_f(n, m);

  // Fill entries of  $\mathbf{X}$ given $f$.
  std::cout << "--> Evaluating f at the indices of X..." << std::endl;
  eval(X, f);
  std::cout << X << std::endl;

  // Construct $\mathbf{A}_small$.
  constexpr index_t n_small = 3, m_small = 3;
  std::cout << "--> Building matrix A..." << std::endl;

  Eigen::SparseMatrix<double> A_small =
      build_matrix(S, shape_t{n_small, m_small});
  std::cout << Eigen::MatrixXd(A_small) << std::endl;

  // Construct $\mathbf{A}$ for more tests.
  Eigen::SparseMatrix<double> A = build_matrix(S, shape_t{n, m});

  // Compute $\mathbf{A} vec(\mathbf{X})$.
  std::cout << "--> Computing A*vec(X)..." << std::endl;
  mult(A, X, Y);
  std::cout << Y << std::endl;

  // Compute $\mathbf{A}^{-1} vec(\mathbf{X})$.
  std::cout << "--> Solving A*vec(X) = vec(Y)..." << std::endl;
  solve(A, Y, X);
  std::cout << "--> The Frobenius norm of the matrix X is " << X.norm()
            << std::endl;
}