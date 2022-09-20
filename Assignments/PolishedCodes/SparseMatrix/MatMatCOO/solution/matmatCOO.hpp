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

  // of course, if nnz(A) can be estimated a priori, space should be reserved in
  // triplets
  for (unsigned int i = 0; i < m; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      if (A(i, j) != 0) {
        Trip triplet(i, j, A(i, j));
        triplets.push_back(triplet);
      }
    }
  }
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
  for (auto const& a : A) {
    for (auto const& b : B) {
      if (a.col() == b.row()) {
        Trip triplet(a.row(), b.col(), a.value() * b.value());
        C.push_back(triplet);
      }
    }
  }
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

  // sort A along columns
  std::sort(A.begin(), A.end(),
            [](const Trip& a1, const Trip& a2) { return a1.col() < a2.col(); });
  // and B along rows
  std::sort(B.begin(), B.end(),
            [](const Trip& b1, const Trip& b2) { return b1.row() < b2.row(); });

  std::size_t i_A = 0, i_B = 0;
  std::set<unsigned int> intersect;
  while (i_A != A.size() && i_B != B.size()) {
    if (A[i_A].col() < B[i_B].row()) {
      ++i_A;
    } else if (B[i_B].row() < A[i_A].col()) {
      ++i_B;
    } else {
      intersect.insert(A[i_A].col());
      ++i_A;
      ++i_B;
    }
  }

  // iteration along intersection of nonempty columns in A and rows in B
  TripVec::iterator A_idx = A.begin();
  TripVec::iterator B_idx = B.begin();
  for (auto i : intersect) {
    A_idx = std::find_if(A_idx, A.end(),
                         [i](const Trip& a) { return a.col() == i; });
    B_idx = std::find_if(B_idx, B.end(),
                         [i](const Trip& b) { return b.row() == i; });

    TripVec::iterator A_it;
    TripVec::iterator B_it;
    for (A_it = A_idx; A_it != A.end(); ++A_it) {
      if (A_it->col() != i) {
        break;
      } else {
        for (B_it = B_idx; B_it != B.end(); ++B_it) {
          if (B_it->row() != i) {
            break;
          } else {
            // actual multiplication
            Trip triplet(A_it->row(), B_it->col(),
                         A_it->value() * B_it->value());
            C.push_back(triplet);
          }
        }
      }
    }
    A_idx = A_it;
    B_idx = B_it;
  }
  // END
  return C;
}
/* SAM_LISTING_END_2 */

#endif