#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

// TODO: write Triplet if you want to use this structure for TripletMatrix

/* @brief Defines a matrix stored in triplet format (using the Triplet<scalar> class.
 * Triplets may be duplicated and in *any* order. If there is a multiple triplet for (row,col) pair, we assume
 * that the values are intended to be added together
 * Also stores dimension
 * @tparam scalar represents the scalar type of the matrix (and of the triplet) (e.g. double)
 */
/* SAM_LISTING_BEGIN_1 */
template <class scalar>
struct TripletMatrix {
    // TODO: put members here

  /* @brief Converts *this to an eigen (dense) function
   * Loops over all triplets and add value to zero matrix
   * WARNING: May fill in a lot of nonzeros if n,m large
   * @return Matrix of Dynamic size and scalar type
   */
  Matrix<scalar, -1, -1> densify() const {
    // TODO: return the "dense" version of "*this"
  }
};
/* SAM_LISTING_END_1 */

/* @brief Structure holding a pair column index-value to be used in CRS format
 * Provides handy constructor and comparison operators.
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 */
/* SAM_LISTING_BEGIN_2 */
    // TODO: write ColValPair if you want to use this structure for CRSMatrix
/* SAM_LISTING_END_2 */

/* @brief Defines a matrix stored in CRS format (using the ColValPair<scalar> struct.
 * The row_pt contains the data, indexed by row and column position
 * Also stores dimension
 * @tparam scalar represents the scalar type of the matrix (and of the ColValPair) (e.g. double)
 */
/* SAM_LISTING_BEGIN_3 */
template <class scalar>
struct CRSMatrix {
    // TODO: imsert members here

  /* @brief Converts *this to an eigen (dense) function
   * Loops over all rows and add value at col to zero matrix
   * WARNING: May fill in a lot of nonzeros if n,m large
   * @return Matrix of Dynamic size and scalar type
   */
  Matrix<scalar, -1, -1> densify() const {
    // TODO: convert "*this" to a dense matrix
  }
};
/* SAM_LISTING_END_3 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format
 * No assumption is made on the triplets, may be unsorted and/or duplicated
 * In case of duplicated triplets, values are added toghether
 * The output CRS matrix may be empty (in that case T = C) or may be already filled with values (in which case C += T)
 * This version inserts the ColValParis already sorted in the list
 * Complexity: Loops over all triplets (lets say k) and look up over all columns of the row (say $n_i$), performing a
 * linear search and an insertion (complexity $O(n_i)$)
 * Assuming $n_i$ is bounded by $n$ small, complexity is $k*n$, otherwise $k^2$
 */
/* SAM_LISTING_BEGIN_4 */
template <class scalar>
void tripletToCRS(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
    // TODO: conversion function
}
/* SAM_LISTING_END_4 */


/* @brief overload of operator << for output of Triplet Matrix (debug).
 * WARNING: Uses densify() so there may be a lot of fill-in
 * @param o standard output stream
 * @param S matrix in Triplet matrix format
 * @return a ostream o, s.t. you can write o << A << B;
 */
std::ostream & operator<<(std::ostream& o, const TripletMatrix<double>& S) {
  return o << S.densify();
}

/* @brief overload of operator << for output of CRS Matrix (debug).
 * WARNING: Uses densify() so there may be a lot of fill-in
 * @param o standard output stream
 * @param S matrix in CRS matrix format
 * @return a ostream o, s.t. you can write o << A << B;
 */
std::ostream & operator<<(std::ostream& o, const CRSMatrix<double>& S) {
  return o << S.densify();
}

int main() {
  //// Correctness test
  std::size_t nrows = 7, ncols = 5, ntriplets = 9;

  TripletMatrix<double> T;
  CRSMatrix<double> C, D;

    // TODO: construct T here

  for(std::size_t i = 0; i < ntriplets; ++i) {
    // TODO: Use this loop to push back triplets in your matrix
    // Insert triplet with arguments: (rand() % nrows, rand() % ncols, rand() % 1000))
  }

  std::cout << "***Test conversion with random matrices***" << std::endl;
  tripletToCRS(T, C);
    // TODO: if you implemented densify(), compute Frobenius norm of T - C
  // std::cout << "T = " << std::endl
            // << T << std::endl;
  // std::cout << "C = " << std::endl
            // << C << std::endl;

/* SAM_LISTING_END_5 */
}
