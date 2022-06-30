#ifndef MATMATCOO_HPP
#define MATMATCOO_HPP

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <algorithm>
#include <set>
#include <vector>

using Trip = Eigen::Triplet<double>;
using TripVec = std::vector<Trip>;

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