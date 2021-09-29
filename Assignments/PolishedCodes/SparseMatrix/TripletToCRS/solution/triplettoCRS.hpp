#ifndef TRIPLETTOCRS_HPP
#define TRIPLETTOCRS_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <ostream>
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
  // TODO: (2-13.a – optional) Complete this struct if you want to use it for
  // TripletMatrix.
  // START
  // default constructor
  Triplet() = delete;

  // constructor taking indices i and j and a value v
  Triplet(std::size_t i_, std::size_t j_, scalar v_) : i(i_), j(j_), v(v_) {}

  std::size_t i, j;  // row and col
  scalar v;          // value in (row, col)
  // END
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
  // TODO: (2-13.a) Store sizes and indices using appropriate types.
  // START
  std::size_t rows, cols;  // sizes: nrows and ncols
  std::vector<Triplet<scalar>> triplets;
  // END

  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> densify() const;
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
  // TODO: (2-13.b) Define the member variables for a CRS matrix using
  // appropriate types.
  // START
  std::size_t rows, cols;  // sizes: nrows and ncols
  std::vector<scalar> val;
  std::vector<std::size_t> col_ind;
  std::vector<std::size_t> row_ptr;
  // END

  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> densify() const;
};
/* SAM_LISTING_END_2 */

/**
 * @brief Convert *this (TripletMatrix) to an Eigen (dense) matrix:
 * Loop over all triplets and add values to a zero matrix.
 * WARNING: May fill in many nonzeros if large n,m
 *
 * @tparam scalar Type of the matrix and triplets (e.g. double)
 * @return Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix of
 * double type
 */
/* SAM_LISTING_BEGIN_3 */
template <typename scalar>
Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>
TripletMatrix<scalar>::densify() const {
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> M;

  // TODO: (2-13.c – optional) Convert the sparse TripletMatrix to a dense
  // Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>.
  // START
  M = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows, cols);

  for (auto&& it : triplets) {
    M(it.i, it.j) += it.v;
  }
  // END
  return M;
}
/* SAM_LISTING_END_3 */

/**
 * @brief Convert *this (CRSMatrix) to an Eigen (dense) matrix:
 * Loop over all rows and add values at col to zero matrix.
 * WARNING: May fill in many nonzeros if large n,m
 *
 * @tparam scalar Type of the matrix and triplets (e.g. double)
 * @return Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix of
 * double type
 */
/* SAM_LISTING_BEGIN_4 */
template <typename scalar>
Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>
CRSMatrix<scalar>::densify() const {
  Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> M;

  // TODO: (2-13.c – optional) Convert the sparse CRSMatrix to a dense
  // Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>.
  // START
  M = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows, cols);
  std::vector<std::size_t> row_ptr_end = row_ptr;
  row_ptr_end.push_back(col_ind.size());

  for (std::size_t i = 0; i < row_ptr.size(); ++i) {
    for (std::size_t j = row_ptr_end[i]; j < row_ptr_end[i + 1]; ++j) {
      M(i, col_ind[j]) = val[j];
    }
  }
  // END
  return M;
}
/* SAM_LISTING_END_4 */

/**
 * @brief Structure holding a pair column index-value to be used in CRS format
 * Provides handy constructor and comparison operators.
 *
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 */
/* SAM_LISTING_BEGIN_7 */
template <typename scalar>
struct ColValPair {
  ColValPair(std::size_t col_, scalar v_) : col(col_), v(v_) {}

  /**
   * @brief Comparison operator for std::sort and std::lower\_bound: Basic
   * sorting operator < to use with std::functions (i.e. for ordering according
   * to first component col). We keep the column values sorted, either by
   * sorting after insertion of by sorted insertion.
   *
   * @param other the ColValPair to compare to
   * @return true this->col < other.col.
   * @return false this->col >= other.col.
   */
  bool operator<(const ColValPair& other) const {
    return this->col < other.col;
  }

  std::size_t col;  // col index
  scalar v;         // scalar value at col
};
/* SAM_LISTING_END_7 */

/**
 * @brief Convert 'preCRS' to proper CRS format (CRSMatrix)
 *
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 * @param preCRS Vector of rows of pairs (col\_ind, val).
 * The original matrix does not have rows only made by zeros.
 * @param C Proper CRSMatrix format
 */
/* SAM_LISTING_BEGIN_8 */
template <typename scalar>
void properCRS(const std::vector<std::vector<ColValPair<scalar>>>& preCRS,
               CRSMatrix<scalar>& C) {
  C.row_ptr.push_back(0);
  for (auto& row : preCRS) {
    // If one whole row is empty,
    // the same index is stored in 'row\_ptr' twice.
    C.row_ptr.push_back(C.row_ptr.back() + row.size());
    for (auto& col_val : row) {
      C.col_ind.push_back(col_val.col);
      C.val.push_back(col_val.v);
    }
  }
  C.row_ptr.pop_back();
}
/* SAM_LISTING_END_8 */

/**
 * @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or
 * duplicated, but it is imposed that the original matrix does not have rows
 * only made by zeros. In case of duplicated triplets, values are added
 * together. The output CRS matrix may be empty (in which case T = C) or may be
 * already filled with values (in which case C += T). This version inserts the
 * pairs ColVal already sorted in the list Complexity: Loop over all $k$
 * triplets and look up over all $n_i$ columns of row $i$, performing a linear
 * search and an insertion (complexity $O(n_i)$). If $n_i$ is bounded by a small
 * $n$, the complexity is $k*n$, otherwise $k^2$.
 *
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 * @param T input matrix
 * @param C output matrix in CRS format
 */
/* SAM_LISTING_BEGIN_5 */
template <typename scalar>
void tripletToCRS_insertsort(const TripletMatrix<scalar>& T,
                             CRSMatrix<scalar>& C) {
  // TODO: (2-13.d – alternative 1) Convert the TripletMatrix to a CRSMatrix by
  // inserting the pair (col, val) already sorted into the data structure. You
  // may use the provided struct ColValPair and the function properCRS().
  // START
  C.rows = T.rows;
  C.cols = T.cols;
  // Temporary "quasi-CRS" format for the conversion
  std::vector<std::vector<ColValPair<scalar>>> preCRS;
  preCRS.resize(C.rows);

  // Loop over all triplets
  for (auto triplet = T.triplets.begin(); triplet != T.triplets.end();
       ++triplet) {
    // Store row (containing all ColValPair for this row) inside row\_ptr
    std::vector<ColValPair<scalar>>& row =
        preCRS.at(triplet->i);  // Notice the reference!

    // Create a ColValPair that must be inserted somewhere
    ColValPair<scalar> col_val(triplet->j, triplet->v);

    // Find the place (as iterator) where to insert col\_val, i.e. the first
    // position $lb_it$ s.t col\_val < row
    // WARNING: Costly call, which returns an iterator to an element
    // (ColValPair) in row
    auto lb_it = std::lower_bound(row.begin(), row.end(), col_val);

    // If lower bound has already a col (and is not end()) with some value,
    // just add the value to the column, othervise insert it
    if (lb_it != row.end() && lb_it->col == triplet->j) {
      lb_it->v += triplet->v;
    } else {
      // WARNING: Costly call (loop over all rows in the worst-case scenario)
      row.insert(lb_it, col_val);
    }
  }

  // Convert 'preCRS' to proper CRS format (CRSMatrix)
  properCRS(preCRS, C);
  // END
}
/* SAM_LISTING_END_5 */

/**
 * @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or
 * duplicated, but it is imposed that the original matrix does not have rows
 * only made by zeros. In case of duplicated triplets, values are added
 * together. The output CRS matrix may be empty (in which case T = C) or may be
 * already filled with values (in which case C += T). This version inserts all
 * the triplets (push_back), then sorts each row and cumsum all the duplicated
 * col values. Complexity: Loop over all $k$ triplets and do a cheap (amortized)
 * o(1) push_back (complexity $k*1$). Then sorts an array of rows (i-th
 * quicksort complexity $k * log(k)$).
 *
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 * @param T input matrix
 * @param C output matrix in CRS format
 */
/* SAM_LISTING_BEGIN_6 */
template <typename scalar>
void tripletToCRS_sortafter(const TripletMatrix<scalar>& T,
                            CRSMatrix<scalar>& C) {
  // TODO: (2-13.d – alternative 2) Convert the TripletMatrix to a CRSMatrix by
  // inserting all triplets, sorting each row and cumulatively summing all the
  // duplicated column values. You may use the struct ColValPair and the
  // function properCRS().
  // START
  C.rows = T.rows;
  C.cols = T.cols;
  // Temporary "quasi-CRS" format for the conversion
  std::vector<std::vector<ColValPair<scalar>>> preCRS;
  preCRS.resize(C.rows);

  // Loop over all triplets and push them at the right place (cheap)
  for (auto triplet = T.triplets.begin(); triplet != T.triplets.end();
       ++triplet) {
    ColValPair<scalar> cp(triplet->j, triplet->v);
    preCRS.at(triplet->i).push_back(cp);
  }
  // Loop over all rows, sort them (according to col) and sum duplicated values
  for (auto row_it = preCRS.begin(); row_it != preCRS.end(); ++row_it) {
    // WARNING: costly call in this case
    std::sort(row_it->begin(), row_it->end());

    // Sum duplicated cols
    // NOTE: iterating backwards should be faster
    for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
      // Check that next col is not at the end and that has same col index
      if ((col_it + 1) != row_it->end() && (col_it + 1)->col == col_it->col) {
        // Last col with same value (starting from next col)
        auto last_col_it = col_it + 1;
        while (last_col_it != row_it->end() &&
               last_col_it->col == col_it->col) {
          col_it->v += last_col_it->v;
          ++last_col_it;
        }

        // Remove range from next col to last col with same value
        // WARNING: std::vector keeps the data as c-array, so mildly costly call
        // (from col\_it to end(), at most M)
        row_it->erase(col_it + 1, last_col_it);
      }
    }
  }

  // Convert 'preCRS' to proper CRS format (CRSMatrix)
  properCRS(preCRS, C);
  // END
}
/* SAM_LISTING_END_6 */

/**
 * @brief overload of operator << for output of Triplet Matrix (debug).
 * WARNING: Uses densify() so there may be a lot of fill-in
 *
 * @param o Standard output stream
 * @param S Matrix in Triplet matrix format
 * @return o std::ostream s.t. you can write o << A << B;
 */
std::ostream& operator<<(std::ostream& o, const TripletMatrix<double>& S) {
  return o << S.densify();
}

/**
 * @brief overload of operator << for output of CRS Matrix (debug).
 * WARNING: Uses densify() so there may be a lot of fill-in
 *
 * @param o Standard output stream
 * @param S Matrix in CRS matrix format
 * @return o std::ostream s.t. you can write o << A << B;
 */
std::ostream& operator<<(std::ostream& o, const CRSMatrix<double>& S) {
  return o << S.densify();
}

#endif