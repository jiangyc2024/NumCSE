#pragma once

#include "totriplets.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
SparseMatrix<double> initA(unsigned int n) {
  unsigned int rows = n * (n - 1) / 2;
  unsigned int cols = n - 1;
  SparseMatrix<double> A(rows, cols);
  // TO DO (0-1.b)
  // START
  // Tell Eigen the number of non-zeros per column
  A.reserve(VectorXi::Constant(n - 1, n - 1));

  // Loop through the $\cob{n-1}$ columns of the matrix
  for (unsigned int i = 0; i < cols; ++i) {
    unsigned int offset = cols * i - i * (i - 1) / 2;
    for (unsigned int k = 0; k < cols - i; ++k) {
      // Set the -1's in $\cob{\VC_k}$ from \prbeqref{eq:ppos}
      A.insert(offset + k, i) = -1;
      // Set the identity blocks of $\cob{\VC_k}$
      if (k < cols - i - 1) {
        A.insert(offset + k, i + k + 1) = 1;
      }
    }
  }
  // END
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
VectorXd solveExtendedNormalEquations(const MatrixXd &D) {
  VectorXd x;
  unsigned int n = D.cols();
  unsigned int m = n * (n - 1) / 2;
  // TO DO (0-1.c)
  // START
  // initialize right hand side vector $\cob{\Vb}$ as in \prbeqref{eq:ppos}
  VectorXd b_ext = VectorXd::Zero(m + n - 1);
  int curr = 0;
  for (unsigned int i = 0; i < n - 1; ++i) {
    for (unsigned int k = i + 1; k < n; ++k) {
      b_ext(curr) = D(i, k);
      ++curr;
    }
  }
  // initialize matrix A and convert to triplet format
  SparseMatrix<double> A = initA(n);
  std::vector<Eigen::Triplet<double>> tripletsA = convertToTriplets(A);
  // define triplets for extended normal equation
  std::vector<Eigen::Triplet<double>> tripletsAext;
  for (unsigned int k = 0; k < m; ++k) {
    // - I block
    tripletsAext.push_back(Eigen::Triplet<double>(k, k, -1.0));
  }
  for (auto t : tripletsA) {
    // A block
    tripletsAext.push_back(
        Eigen::Triplet<double>(t.row(), m + t.col(), t.value()));
    // A.transpose() block
    tripletsAext.push_back(
        Eigen::Triplet<double>(m + t.col(), t.row(), t.value()));
  }

  SparseMatrix<double> A_ext(n - 1 + m, n - 1 + m);
  A_ext.setFromTriplets(tripletsAext.begin(), tripletsAext.end());
  // solve extended system and extract the last n-1 entries
  // of the solution vector
  SparseLU<SparseMatrix<double>> solver(A_ext);
  x = solver.solve(b_ext).tail(n - 1);
  // END
  return x;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
VectorXd solveNormalEquations(const MatrixXd &D) {
  VectorXd x;
  // TO DO (0-1.e)
  // START
  unsigned int n = D.cols();
  unsigned int m = n * (n - 1) / 2;
  // initialize right hand side
  VectorXd b = VectorXd::Zero(m);
  int curr = 0;
  for (unsigned int i = 0; i < n - 1; ++i) {
    for (unsigned int k = i + 1; k < n; ++k) {
      b(curr) = D(i, k);
      ++curr;
    }
  }
  // Direct computation of A^T*b
  SparseMatrix<double> AT = initA(n).transpose();
  VectorXd ATb = AT * b;
  // Use explicit expression of the inverse of A^T*A
  // to compute the solution
  x = (MatrixXd::Identity(n - 1, n - 1) + MatrixXd::Ones(n - 1, n - 1)) * ATb /
      n;
  // END
  return x;
}

/* SAM_LISTING_END_2 */
