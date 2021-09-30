#ifndef TRIPLETTOCRS_HPP
#define TRIPLETTOCRS_HPP

#include <Eigen/Dense>
#include <algorithm>
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
  for (auto& triplet : M.triplets) {
    dense(std::get<0>(triplet), std::get<1>(triplet)) += std::get<2>(triplet);
  }
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
  std::vector<std::size_t> row_ptr_end = M.row_ptr;
  row_ptr_end.push_back(M.col_ind.size());

  for (std::size_t i = 0; i < M.row_ptr.size(); ++i) {
    for (std::size_t j = row_ptr_end[i]; j < row_ptr_end[i + 1]; ++j) {
      dense(i, M.col_ind[j]) = M.val[j];
    }
  }
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
  // sort the triplets according to their indices over a copy of the vector
  auto triplets_copy = T.triplets;
  auto sorting_predicate =
      [](const std::tuple<std::size_t, std::size_t, SCALAR>& a,
         const std::tuple<std::size_t, std::size_t, SCALAR>& b) {
        return (std::get<0>(a) < std::get<0>(b)) ||
               ((std::get<0>(a) == std::get<0>(b)) &&
                (std::get<1>(a) < std::get<1>(b)));
      };
  std::sort(triplets_copy.begin(), triplets_copy.end(), sorting_predicate);

  // count the number of non-zeroes per row
  crs.row_ptr = std::vector<std::size_t>(crs.n_rows + 1, 0);
  std::size_t k = 0;
  std::size_t last_j = -1;
  std::size_t last_i = -1;
  for (std::size_t i = 1; i <= crs.n_rows; ++i) {
    while (k < triplets_copy.size() && std::get<0>(triplets_copy[k]) == i - 1) {
      // do only when no duplicate index pair
      if (std::get<1>(triplets_copy[k]) != last_j ||
          std::get<0>(triplets_copy[k]) != last_i) {
        last_j = std::get<1>(triplets_copy[k]);
        last_i = std::get<0>(triplets_copy[k]);
        ++crs.row_ptr[i];
      }
      ++k;
    }
    crs.row_ptr[i] += crs.row_ptr[i - 1];
  }

  // now, we know nnz; reserve in advance
  crs.val.reserve(crs.row_ptr[crs.n_rows]);
  crs.col_ind.reserve(crs.row_ptr[crs.n_rows]);

  crs.val.push_back(std::get<2>(triplets_copy[0]));
  crs.col_ind.push_back(std::get<1>(triplets_copy[0]));
  std::size_t i = 0;
  for (k = 1; k < triplets_copy.size(); ++k) {
    // repeated index pairs
    if (std::get<0>(triplets_copy[k]) == std::get<0>(triplets_copy[k - 1]) &&
        std::get<1>(triplets_copy[k]) == std::get<1>(triplets_copy[k - 1])) {
      crs.val[i] += std::get<2>(triplets_copy[k]);
    } else {
      crs.val.push_back(std::get<2>(triplets_copy[k]));
      crs.col_ind.push_back(std::get<1>(triplets_copy[k]));
      ++i;
    }
  }
  // END
  return crs;
}
/* SAM_LISTING_END_4 */

#endif