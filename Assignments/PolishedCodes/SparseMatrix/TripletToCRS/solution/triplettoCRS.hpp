#ifndef TRIPLETTOCRS_HPP
#define TRIPLETTOCRS_HPP

#include <Eigen/Dense>
#include <vector>

/**
 * @brief Structure holding a triplet in format (row, col, value)
 * All members are public. Used in the TripletMatrix class to provide a nice
 * wrapper for a triplet.
 *
 * @tparam scalar Type of the triplet (e.g. double)
 */
/* SAM_LISTING_BEGIN_0 */
template <typename scalar>
struct Triplet {
  // default constructor
  Triplet() = delete;

  // constructor taking indices i and j and a value v
  Triplet(std::size_t i_, std::size_t j_, scalar v_) : i(i_), j(j_), v(v_) {}

  std::size_t i, j;  // row and col
  scalar v;          // value in (row, col)
};
/* SAM_LISTING_END_0 */

/**
 * @brief Defines a matrix stored in triplet format (using the Triplet<scalar>
 * class). Triplets may be duplicated and in *any* order. If there is a multiple
 * triplet for (row,col) pair, we assume that the values are intended to be
 * added together. Dimensions are also stored to simplify the code.
 *
 * @tparam scalar Type of the matrix and triplets (e.g. double)
 */
/* SAM_LISTING_BEGIN_1 */
template <typename scalar>
struct TripletMatrix {
  std::size_t rows, cols;  // sizes: nrows and ncols
  std::vector<Triplet<scalar>> triplets;

  Eigen::MatrixXd densify() const;
};
/* SAM_LISTING_END_1 */

/**
 * @brief Defines a matrix stored in CRS format.
 * Dimensions are also stored to simplify the code.
 *
 * @tparam scalar Type of the matrix and CRS vectors (e.g. double)
 */
/* SAM_LISTING_BEGIN_2 */
template <typename scalar>
struct CRSMatrix {
  std::size_t rows, cols;  // sizes: nrows and ncols
  std::vector<scalar> val;
  std::vector<std::size_t> col_ind;
  std::vector<std::size_t> row_ptr;

  Eigen::MatrixXd densify() const;
};
/* SAM_LISTING_END_2 */

/**
 * @brief
 *
 * @tparam scalar
 * @return Eigen::MatrixXd
 */
/* SAM_LISTING_BEGIN_3 */
template <typename scalar>
Eigen::MatrixXd TripletMatrix<scalar>::densify() const {
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(rows, cols);

  for (auto&& it : triplets) {
    M(it->i, it->j) += it->v;
  }

  return M;
}
/* SAM_LISTING_END_3 */

#endif