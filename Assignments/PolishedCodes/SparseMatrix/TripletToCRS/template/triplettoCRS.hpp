#ifndef TRIPLETTOCRS_HPP
#define TRIPLETTOCRS_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

/**
 * @brief Defines a matrix stored in triplet format (using the Triplet<scalar>
 * class). Triplets may be duplicated and in *any* order. If there is a multiple
 * triplet for (row,col) pair, we assume that the values are intended to be
 * added together. Dimensions are also stored to simplify the code.
 *
 * @tparam SCALAR Type of the matrix values (e.g. double)
 */
/* SAM_LISTING_BEGIN_0 */
template <typename SCALAR>
struct TripletMatrix {
  std::size_t n_rows{0};  // Number of rows
  std::size_t n_cols{0};  // Number of columns
  std::vector<std::tuple<std::size_t, std::size_t, SCALAR>> triplets;
};
/* SAM_LISTING_END_0 */

/**
 * @brief Defines a matrix stored in CRS format.
 * Dimensions are also stored to simplify the code.
 *
 * @tparam SCALAR Type of the matrix values (e.g. double)
 */
/* SAM_LISTING_BEGIN_1 */
template <typename SCALAR>
struct CRSMatrix {
  std::size_t n_rows{0};             // Number of rows
  std::size_t n_cols{0};             // Number of columns
  std::vector<SCALAR> val;           // Value array
  std::vector<std::size_t> col_ind;  // Column indices
  std::vector<std::size_t> row_ptr;  // Row pointers
};
/* SAM_LISTING_END_1 */

/**
 * @brief Converts a TripletMatrix to a dense Eigen matrix.
 *
 * @tparam SCALAR Type of the matrix values (e.g. double)
 * @param M the TripletMatrix to convert
 * @return Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> the densified
 * matrix
 */
/* SAM_LISTING_BEGIN_2 */
template <typename SCALAR>
Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> densify(
    const TripletMatrix<SCALAR>& M) {
  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dense =
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>::Zero(M.n_rows,
                                                                  M.n_cols);

  // TODO: (2-13.a) Densify the TripletMatrix.
  // START

  // END
  return dense;
}
/* SAM_LISTING_END_2 */

/**
 * @brief Converts a CRSMatrix to a dense Eigen matrix.
 *
 * @tparam SCALAR Type of the matrix values (e.g. double)
 * @param M the CRSMatrix to convert
 * @return Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> the densified
 * matrix
 */
/* SAM_LISTING_BEGIN_3 */
template <typename SCALAR>
Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> densify(
    const CRSMatrix<SCALAR>& M) {
  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dense =
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>::Zero(M.n_rows,
                                                                  M.n_cols);

  // TODO: (2-13.a) Densify the CRSMatrix.
  // START

  // END
  return dense;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Converts a TripletMatrix to a CRSMatrix.
 *
 * @tparam SCALAR Type of the matrix values (e.g. double)
 * @param T the TripletMatrix to convert
 * @return CRSMatrix<SCALAR> the converted matrix
 */
/* SAM_LISTING_BEGIN_4 */
template <typename SCALAR>
CRSMatrix<SCALAR> tripletToCRS(const TripletMatrix<SCALAR>& T) {
  CRSMatrix<SCALAR> crs;
  crs.n_rows = T.n_rows;
  crs.n_cols = T.n_cols;

  // TODO: (2-13.b) Convert the TripletMatrix to a CRSMatrix. Make sure that you
  // handle repeated index pairs.
  // START

  // END
  return crs;
}
/* SAM_LISTING_END_4 */

#endif