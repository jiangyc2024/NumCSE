#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Structure holding a triplet in format (row, col, value)
 * For convenience, we provide a void constructor and a init constructor. All members are public
 * Used in the TripletMatrix class to provide a nice wrapper for a triplet
 * @tparam scalar Type of the triplet (e.g. double)
 */
/* SAM_LISTING_BEGIN_0 */
template <typename scalar>
struct Triplet {
  // TODO: write Triplet if you want to use this structure for TripletMatrix
};
/* SAM_LISTING_END_0 */

/* @brief Defines a matrix stored in triplet format (using the Triplet<scalar> class).
 * Triplets may be duplicated and in *any* order.
 * If there is a multiple triplet for (row,col) pair, we assume that the values are intended to be added together.
 * Dimensions are also stored to simplify the code.
 * @tparam scalar Type of the matrix and triplets (e.g. double)
 */
/* SAM_LISTING_BEGIN_1 */
template <typename scalar>
struct TripletMatrix {
  // TODO: put members here

  MatrixXd densify() const;
};
/* SAM_LISTING_END_1 */

/* @brief Defines a matrix stored in CRS format.
 * Dimensions are also stored to simplify the code.
 * @tparam scalar Type of the matrix and CRS vectors (e.g. double)
 */
/* SAM_LISTING_BEGIN_2 */
template <typename scalar>
struct CRSMatrix {
    // TODO: put members here

  MatrixXd densify() const;
};
/* SAM_LISTING_END_2 */

/* @brief Convert *this (TripletMatrix) to an Eigen (dense) matrix:
 * Loop over all triplets and add values to a zero matrix.
 * WARNING: May fill in many nonzeros if large n,m
 * @return Matrix of 'scalar' type
 */
/* SAM_LISTING_BEGIN_3 */
template <class scalar>
MatrixXd TripletMatrix<scalar>::densify() const {
  MatrixXd M;
// TODO: return the "dense" version of "*this" (TripletMatrix)

  return M;
}
/* SAM_LISTING_END_3 */

/* @brief Convert *this (CRSMatrix) to an Eigen (dense) matrix:
 * Loop over all rows and add values at col to zero matrix.
 * WARNING: May fill in many nonzeros if large n,m
 * @return Matrix of 'scalar' type
 */
/* SAM_LISTING_BEGIN_4 */
template <typename scalar>
MatrixXd CRSMatrix<scalar>::densify() const {
  MatrixXd M;
    // TODO: convert "*this" to a dense matrix (CRSMatrix)

  return M;
}
/* SAM_LISTING_END_4 */

/* @brief Structure holding a pair column index-value to be used in CRS format
 * Provides handy constructor and comparison operators.
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 */
/* SAM_LISTING_BEGIN_7 */
template <typename scalar>
struct ColValPair {
  ColValPair(size_t col_, scalar v_)
    : col(col_), v(v_) { }

  /* Comparison operator for std::sort and std::lower\_bound:
     Basic sorting operator < to use with std::functions (i.e. for ordering according to first component col).
     We keep the column values sorted, either by sorting after insertion of by sorted insertion.
     It returns true if this->col < other.col.
   */
  bool operator<(const ColValPair& other) const {
    return this->col < other.col;
  }

  size_t col; // Col index
  scalar v; // Scalar value at col
};
/* SAM_LISTING_END_7 */

/* @brief Convert 'preCRS' to proper CRS format (CRSMatrix)
 * @param[in] preCRS Vector of rows of pairs (col\_ind, val).
 * The original matrix does not have rows only made by zeros.
 * @param[out] C Proper CRSMatrix format
 */
/* SAM_LISTING_BEGIN_8 */
template <typename scalar>
void properCRS(const std::vector< std::vector< ColValPair<scalar> > >& preCRS, CRSMatrix<scalar>& C) {

  C.row_ptr.push_back(0);
  for(auto & row : preCRS) {
	// If one whole row is empty,
	// the same index is stored in 'row\_ptr' twice.
	C.row_ptr.push_back(C.row_ptr.back() + row.size());
	for(auto & col_val : row) {
	  C.col_ind.push_back(col_val.col);
      C.val.push_back(    col_val.v);
	}
  }
  C.row_ptr.pop_back();
}
/* SAM_LISTING_END_8 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or duplicated,
 * but it is imposed that the original matrix does not have rows only made by zeros.
 * In case of duplicated triplets, values are added toghether.
 * The output CRS matrix may be empty (in which case T = C) or may be already filled with values (in which case C += T).
 * This version inserts the pairs ColVal already sorted in the list
 * Complexity: Loop over all $k$ triplets and look up over all $n_i$ columns of row $i$,
 * performing a linear search and an insertion (complexity $O(n_i)$).
 * If $n_i$ is bounded by a small $n$, the complexity is $k*n$, otherwise $k^2$.
 */
/* SAM_LISTING_BEGIN_5 */
template <typename scalar>
void tripletToCRS_insertsort(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
    // TODO: conversion function (insert sort)
}
/* SAM_LISTING_END_5 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or duplicated,
 * but it is imposed that the original matrix does not have rows only made by zeros.
 * In case of duplicated triplets, values are added toghether.
 * The output CRS matrix may be empty (in which case T = C) or may be already filled with values (in which case C += T).
 * This version inserts all the triplets (push_back), then sorts each row and cumsum all the duplicated col values.
 * Complexity: Loop over all $k$ triplets and do a cheap (amortized) o(1) push_back (complexity $k*1$).
 * Then sorts an array of rows (i-th quicksort complexity $k * log(k)$).
 */
/* SAM_LISTING_BEGIN_6 */
template <typename scalar>
void tripletToCRS_sortafter(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
    // TODO: conversion function (sort after)
}
/* SAM_LISTING_END_6 */

/* @brief overload of operator << for output of Triplet Matrix (debug).
 * WARNING: Uses densify() so there may be a lot of fill-in
 * @param[in] o Standard output stream
 * @param[in] S Matrix in Triplet matrix format
 * @return o std::ostream s.t. you can write o << A << B;
 */
std::ostream & operator<<(std::ostream& o, const TripletMatrix<double>& S) {
  return o << S.densify();
}

/* @brief overload of operator << for output of CRS Matrix (debug).
 * WARNING: Uses densify() so there may be a lot of fill-in
 * @param[in] o Standard output stream
 * @param[in] S Matrix in CRS matrix format
 * @return o std::ostream s.t. you can write o << A << B;
 */
std::ostream & operator<<(std::ostream& o, const CRSMatrix<double>& S) {
  return o << S.densify();
}

int main() {
  // Initialization
  size_t nrows = 7, ncols = 5, ntriplets = 9;
  TripletMatrix<double> T;
  CRSMatrix<double> C, D;

    // TODO: construct T here

  for(size_t i = 0; i < ntriplets; ++i) {
    // TODO: Use this loop to push back triplets in your matrix
    // Insert triplet with arguments: (rand() % nrows, rand() % ncols, rand() % 1000))
  }

  std::cout << "***Test conversion with random matrices***" << std::endl;
  tripletToCRS_insertsort(T, C);
  tripletToCRS_sortafter( T, D);
    // TODO: if you implemented densify(), compute Frobenius norm of $T - C$

}
