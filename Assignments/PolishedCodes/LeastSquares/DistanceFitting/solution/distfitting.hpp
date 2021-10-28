#ifndef DISTFITTING_HPP
#define DISTFITTING_HPP

#include <Eigen/Sparse>

#include "totriplets.hpp"

/**
 * @brief Initializes the system matrix A
 *
 * @param n number of points
 * @return Eigen::SparseMatrix<double> A
 */
/* SAM_LISTING_BEGIN_0 */
Eigen::SparseMatrix<double> initA(unsigned int n) {
  const unsigned int rows = n * (n - 1) / 2;
  const unsigned int cols = n - 1;
  Eigen::SparseMatrix<double> A(rows, cols);
  // TODO: (3-10.b) Initialize the sparse coefficient matrix for
  // the distance fitting problem
  // START
  // Tell Eigen the number of non-zeros per column
  A.reserve(Eigen::VectorXi::Constant(n - 1, n - 1));

  // Loop through the $\cob{n-1}$ columns of the matrix
  for (unsigned int i = 0; i < cols; ++i) {
    const unsigned int offset = cols * i - i * (i - 1) / 2;
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

/**
 * @brief Solves the extended normal equations for the given problem.
 *
 * @param D n times n matrix whose strict upper triangular part contains the
 * distances
 * @return Eigen::VectorXd solution vector
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solveExtendedNormalEquations(const Eigen::MatrixXd& D) {
  Eigen::VectorXd x;
  const unsigned int n = D.cols();
  const unsigned int m = n * (n - 1) / 2;
  // TODO: (3-10.c) Solve the extended normal equations with a sparse solver
  // that Eigen provides
  // START
  // initialize right hand side vector $\cob{\Vb}$ as in \prbeqref{eq:ppos}
  Eigen::VectorXd b_ext = Eigen::VectorXd::Zero(m + n - 1);
  unsigned int curr = 0;
  for (unsigned int i = 0; i < n - 1; ++i) {
    for (unsigned int k = i + 1; k < n; ++k) {
      b_ext(curr) = D(i, k);
      ++curr;
    }
  }
  // initialize matrix A and convert to triplet format
  Eigen::SparseMatrix<double> A = initA(n);
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

  Eigen::SparseMatrix<double> A_ext(n - 1 + m, n - 1 + m);
  A_ext.setFromTriplets(tripletsAext.begin(), tripletsAext.end());
  // solve extended system and extract the last n-1 entries
  // of the solution vector
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A_ext);
  x = solver.solve(b_ext).tail(n - 1);
  // END
  return x;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves the normal equations for the given problem.
 *
 * @param D n times n matrix whose strict upper triangular part contains the
 * distances
 * @return Eigen::VectorXd solution vector
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveNormalEquations(const Eigen::MatrixXd& D) {
  Eigen::VectorXd x;
  const unsigned int n = D.cols();
  const unsigned int m = n * (n - 1) / 2;
  // TODO: (3-10.e) Solve the normal equations by exploiting the
  // Sherman-Morrison-Woodbury formula, i.e. by using your results from the
  // previous subproblem
  // START
  // initialize right hand side
  Eigen::VectorXd b = Eigen::VectorXd::Zero(m);
  unsigned int curr = 0;
  for (unsigned int i = 0; i < n - 1; ++i) {
    for (unsigned int k = i + 1; k < n; ++k) {
      b(curr) = D(i, k);
      ++curr;
    }
  }
  // Direct computation of A^T*b
  Eigen::SparseMatrix<double> AT = initA(n).transpose();
  Eigen::VectorXd ATb = AT * b;

  // Alternative version

  // Initialize the extended D
  Eigen::MatrixXd D_full = Eigen::MatrixXd::Zero(n, n);
  D_full.triangularView<Eigen::Upper>() =
      D.triangularView<Eigen::Upper>();  // in case D is not just a upper
                                         // triangular matrix
  Eigen::MatrixXd D_fullT = D_full.transpose();
  D_full = D_full - D_fullT;
  // Calculate ATb in a more efficient way
  Eigen::VectorXd ATb_alt =
      D_full.transpose().block(0, 0, n - 1, n) * Eigen::VectorXd::Ones(n);

  // We use ATb on forward but ATb == ATb_alt

  // Use explicit expression of the inverse of A^T*A
  // to compute the solution
  x = (Eigen::MatrixXd::Identity(n - 1, n - 1) +
       Eigen::MatrixXd::Ones(n - 1, n - 1)) *
      ATb / n;
  // END
  return x;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Tests the normal equation methods.
 *
 * @param D n times n matrix whose strict upper triangular part contains the
 * distances
 * @return true if the results of both methods agree up to tol
 * @return false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
bool testNormalEquations(const Eigen::MatrixXd& D) {
  constexpr double tol = 1e-9;  // tolerance to check for
  // TODO: (3-10.f) Call your implementations of solveExtendedNormalEquations()
  // and solveNormalEquations() and return true, if their results agree.
  // START
  Eigen::VectorXd x_ext = solveExtendedNormalEquations(D);
  Eigen::VectorXd x_fast = solveNormalEquations(D);
  if ((x_ext - x_fast).norm() < x_ext.norm() * tol) return true;
  // END
  return false;
}
/* SAM_LISTING_END_3 */

#endif