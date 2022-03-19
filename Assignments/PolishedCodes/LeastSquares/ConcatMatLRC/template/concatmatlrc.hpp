/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#define _USE_MATH_DEFINES

// Function computing the factors of the economical SVD of a matrix created by
// concatenating two equal-size matrices given in low-rank factorization
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>
eco_svd_concat(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
               const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  // Matrix dimensions and consistency checks
  const int m = A1.rows();
  assert(A2.rows() == m);
  const int n = B1.rows();
  assert(B2.rows() == n);
  const int k = A1.cols();
  assert((k == A2.cols()) && (k == B1.cols()) && (k == B2.cols()));
  // Matrices/vector to be returned by the function
  Eigen::MatrixXd U(m, 2 * k);
  Eigen::VectorXd s(2 * k);
  Eigen::MatrixXd V(2 * n, 2 * k);
  // TODO : Implement the economical SVD of the concatenated matrix 
  // [A1*B1^T A2*B2^T]. The algorithm must have optimal asymptotic complexity.
  // START
  
  // END
  return {U, s, V};
}

// Test for correct implementation of eco_svd_concat
bool test_eco_svd_concat(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
                         const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2) {
  auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
  unsigned m = A1.rows();
  unsigned n = B1.rows();
  // A suitable tolerance for checks
  const double tau = 100 * m * n * std::numeric_limits<double>::epsilon();
  bool out;
  
  // TODO : Write a function that checks whether the output obtained from 
  // eco_svd_concat() complies with the requirements mentioned in exam PDF
  // START
  
  // END 
  return out;
}

// Best low-rank approximation of concatenated low-rank matrices
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
concat_low_rank_best(const Eigen::MatrixXd &A1, const Eigen::MatrixXd &B1,
                     const Eigen::MatrixXd &A2, const Eigen::MatrixXd &B2,
                     double tol) {
  assert((tol > 0) && (tol < 1.0));
  auto [U, s, V] = eco_svd_concat(A1, B1, A2, B2);
  Eigen::MatrixXd AM,BM;
  // TODO : Write a function that computes the best low rank approximation  
  // (minimal rank) of the matrix [A1*B1^T A2*B2^T] and returns it in a 
  // factorized form AM * BM^T
  // START 
 
  // END 
  return {AM, BM};
}
