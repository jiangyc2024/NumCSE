#ifndef MATMATCOO_HPP
#define MATMATCOO_HPP

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <algorithm>
#include <random>
#include <set>
#include <vector>

using Trip = Eigen::Triplet<double>;
using TripVec = std::vector<Trip>;

/**
 * @brief Build matrix $A$ as Eigen::MatrixXd from COO format
 *
 * @param A An $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @return Eigen::MatrixXd The $m \times n$ matrix from vector of triplets $A$
 */
Eigen::MatrixXd COO2Mat(const TripVec& A) {
  unsigned int m = 0, n = 0;
  for (auto const& a : A) {
    if (a.row() > m) {
      m = a.row();
    }
    if (a.col() > n) {
      n = a.col();
    }
  }
  ++m, ++n;  // first index is 0
  Eigen::MatrixXd A_mat = Eigen::MatrixXd::Zero(m, n);

  for (auto const& a : A) {
    A_mat(a.row(), a.col()) += a.value();
  }

  return A_mat;
}

/**
 * @brief Build random binary matrix $A$ as Eigen::MatrixXd
 *
 * @param m Number of desired rows for $A$
 * @param n Number of desired cols for $A$
 * @param d Maximum density ($nnz/(m*n)$) for $A$
 * @return Eigen::MatrixXd An $m \times n$ random binary matrix
 */
Eigen::MatrixXd randMat(unsigned int m, unsigned int n, double d) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, n);
  const unsigned int nnz = std::round(m * n * d);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<unsigned int> distr_row(0, m - 1);
  std::uniform_int_distribution<unsigned int> distr_col(0, n - 1);

  // we allow to draw a couple $(i, j)$ multiple times:
  // density 'd' is an upper bound for the final density of $A$
  for (unsigned int k = 0; k < nnz; ++k) {
    const unsigned int i = distr_row(gen);
    const unsigned int j = distr_col(gen);
    A(i, j) = 1;
  }

  return A;
}

/**
 * @brief Build matrix $A$ in COO format
 *
 * @param A An $m \times n$ matrix
 * @return TripVec The $nnz(A)$-dimensional vector of triplets forming $A$
 */
/* SAM_LISTING_BEGIN_0 */
TripVec Mat2COO(const Eigen::MatrixXd& A) {
  TripVec triplets;
  const unsigned int m = A.rows();
  const unsigned int n = A.cols();

  // TODO: (2-17.b) Convert A into a sparse matrix in COO format.
  // START

  // END
  return triplets;
}
/* SAM_LISTING_END_0 */

/**
 * @brief Compute matrix product $AB$ in COO format: Naive implementation.
 *
 * @param A The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param B The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @return TripVec The $nnz(C)$-dimensional vector of triplets forming matrix $C
 * = AB$
 */
/* SAM_LISTING_BEGIN_1 */
TripVec COOprod_naive(const TripVec& A, const TripVec& B) {
  TripVec C;

  // TODO: (2-17.c) Compute the matrix-matrix product, naive version.
  // START

  // END
  return C;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Compute matrix product $AB$ in COO format: Efficient implementation.
 *
 * @param A The $nnz(A)$-dimensional vector of triplets forming matrix $A$
 * @param B The $nnz(B)$-dimensional vector of triplets forming matrix $B$
 * @return TripVec The $nnz(C)$-dimensional vector of triplets forming matrix $C
 * = AB$x
 */
/* SAM_LISTING_BEGIN_2 */
TripVec COOprod_effic(TripVec& A, TripVec& B) {
  TripVec C;

  // TODO: (2-17.e) Compute the matrix-matrix product efficiently
  // START

  // END
  return C;
}
/* SAM_LISTING_END_2 */

#endif