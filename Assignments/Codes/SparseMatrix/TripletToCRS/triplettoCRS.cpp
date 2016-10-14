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
#if SOLUTION
  // Default constructor
  Triplet()
  : i(0), j(0), v(0) { }

  // Constructor taking indexes i and j and a value v
  Triplet(size_t i_, size_t j_, scalar v_)
    : i(i_), j(j_), v(v_) { }

    size_t i, j; // Row and col
    scalar v; // Value in (row,col)
#else // TEMPLATE
  // TODO: write Triplet if you want to use this structure for TripletMatrix
#endif // TEMPLATE
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
#if SOLUTION
  size_t rows, cols; // Sizes: nrows and ncols
  std::vector<Triplet<scalar> > triplets;
#else // TEMPLATE
  // TODO: put members here
#endif // TEMPLATE

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
#if SOLUTION
  size_t rows, cols; // Sizes: nrows and ncols
  std::vector<scalar_t> val;
  std::vector<size_t> col_ind;
  std::vector<size_t> row_ptr;
#else // TEMPLATE
    // TODO: put members here
#endif // TEMPLATE

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
#if SOLUTION
  // Initialization
  M = MatrixXd::Zero(rows, cols);

  for(auto it = triplets.begin(); it != triplets.end(); ++it) {
    M(it->i, it->j) += it->v;
  }
#else // TEMPLATE
// TODO: return the "dense" version of "*this" (TripletMatrix)
#endif // TEMPLATE

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
#if SOLUTION
// Initialization
  M = MatrixXd::Zero(rows, cols);
  std::vector<size_t> row_ptr_end = row_ptr;
  row_ptr_end.push_back(row_ptr.size());
  
  for(size_t i=0; i<row_ptr.size(); ++i) {
	  for(size_t j=row_ptr_end[i]; j<row_ptr_end[i+1]; ++j) {
		  M(i, col_ind[j]) = val[j];
	  }
  }
#else // TEMPLATE
    // TODO: convert "*this" to a dense matrix (CRSMatrix)
#endif // TEMPLATE

  return M;
}
/* SAM_LISTING_END_4 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or duplicated.
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
#if SOLUTION
  // Temporary "quasi-CRS" format for the conversion
  std::vector< std::vector< ColValPair<scalar> > > preCRS;

  // Copy sizes and reserve memory for rows
  C.rows = T.rows;
  C.cols = T.cols;
  preCRS_col_ind.resize(C.rows);
  preCRS_val.resize(C.rows);

  // Loop over all triplets
  for(auto triplet_it = T.triplets.begin(); triplet_it != T.triplets.end(); ++triplet_it) {
    // Store row (containing ColValPairs for this row) inside row\_pt
  	std::vector<scalar >& row_col_ind = preCRS_col_ind.at(triplet_it->i); // Notice the reference!
  	std::vector<scalar >& row_val     = preCRS_val.at(triplet_it->i); // Notice the reference!
    // Find the place (as iterator) where to insert col\_val, i.e. the first place where col\_val < C.row\_pt.at(i)
    // WARNING: Costly call, which returns an iterator to a ColValPair in row
    auto lb_it = std::lower_bound(row_col_ind.begin(), row_col_ind.end(), triplet_it->j);
    // If lower bound has already a col (and is not end()) with some value, just add the value to the column, othervise insert it
	  size_t i = lb_it - row_col_ind.begin();
    if(lb_it != row_col_ind.end() && *lb_it == triplet_it->j) {
      row_val.at(i) += triplet_it->v;
    } else {
      // WARNING: Costly call (loops over all rows in the worst-case scenario)
      row_col_ind.insert(lb_it, triplet_it->j);
      row_val.insert(row_val.begin()+i, triplet_it->v);
    }
  }
  
  // Convert 'preCRS' to proper CRS format (CRSMatrix)
  C.row_ptr.push_back(0);
  for(size_t i=0; i<preCRS_col_ind.size(); ++i) {
	if(preCRS_col_ind.at(i).size() != 0) {
	  C.row_ptr.push_back(C.row_ptr.back() + preCRS_col_ind.at(i).size() + 1);
	  C.col_ind.insert(C.col_ind.end(), preCRS_col_ind.at(i).begin(), preCRS_col_ind.at(i).end());
      C.val.insert(C.val.end(), preCRS_val.at(i).begin(), preCRS_val.at(i).end());
	}
  }
  C.row_ptr.pop_back();
  
#else // TEMPLATE
    // TODO: conversion function (insert sort)
#endif // TEMPLATE
}
/* SAM_LISTING_END_5 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or duplicated.
 * In case of duplicated triplets, values are added toghether.
 * The output CRS matrix may be empty (in which case T = C) or may be already filled with values (in which case C += T).
 * This version inserts all the triplets (push_back), then sorts each row and cumsum all the duplicated col values.
 * Complexity: Loop over all $k$ triplets and do a cheap (amortized) o(1) push_back (complexity $k*1$).
 * Then sorts an array of rows (i-th quicksort complexity $k * log(k)$).
 */
/* SAM_LISTING_BEGIN_6 */
template <typename scalar>
void tripletToCRS_sortafter(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
#if SOLUTION
  // Copy dimensions and reserve known space
  C.rows = T.rows;
  C.cols = T.cols;
  C.row_ptr.resize(C.rows);

  // Loops over all triplets and push them at the ritgh place (cheap)
  for(auto triplets_it = T.triplets.begin(); triplets_it != T.triplets.end(); ++triplets_it) {
    ColValPair<scalar> cp(triplets_it->j, triplets_it->v);
    C.row_ptr.at(triplets_it->i).push_back(cp);
  }
  // Loops over all rows , sort them (according to col) and sum duplicated values
  for(auto row_it = C.row_ptr.begin(); row_it != C.row_ptr.end(); ++row_it) {
    // WARNING: costly call in this case
    std::sort(row_it->begin(), row_it->end());
    // Sum duplicated cols
    // NOTE: iterating backwards should be faster
    for(auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
      // Check that next col is not at the end and that has same col value
      if((col_it + 1) != row_it->end() && (col_it + 1)->col == col_it->col) {
        // Last col with same value (starting from next col)
        auto last_col_it = col_it + 1;
        while(last_col_it != row_it->end() && last_col_it->col == col_it->col) {
          col_it->v += last_col_it->v;
          ++last_col_it;
        }
        // Remove range from next col to last col with same value
        // WARNING: std::vector keeps the data as c-array, so mildly costly call (from col\_it to end(), at most M)
        row_it->erase(col_it+1, last_col_it);
      }
    }
  }
  
  TODO: Convert to proper CRS
  
#else // TEMPLATE
    // TODO: conversion function (sort after)
#endif // TEMPLATE
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
  CRSMatrix<double> C;

#if SOLUTION
  T.rows = nrows;
  T.cols = ncols;
  T.triplets.reserve(ntriplets); // Always reserve space if you can
#else // TEMPLATE
    // TODO: construct T here
#endif // TEMPLATE

  for(size_t i = 0; i < ntriplets; ++i) {
#if SOLUTION
    // Test unordered triplets, random and possibly repeated triplets
    T.triplets.push_back(Triplet<double>(rand() % nrows, rand() % ncols, rand() % 1000));
#else // TEMPLATE
    // TODO: Use this loop to push back triplets in your matrix
    // Insert triplet with arguments: (rand() % nrows, rand() % ncols, rand() % 1000))
#endif // TEMPLATE
  }

  std::cout << "***Test conversion with random matrices***" << std::endl;
  tripletToCRS(T, C);
#if SOLUTION
  std::cout << "--> Frobenius norm of T - C: " << (T.densify()-C.densify()).norm() << std::endl;
#else // TEMPLATE
    // TODO: if you implemented densify(), compute Frobenius norm of $T - C$
#endif // TEMPLATE

#if SOLUTION
  // Initialization
  bool noruntime = false; // Skip runtime measurments
  if(noruntime) {
      return 0;
  }
  
  std::cout << "***Runtime test***" << std::endl;

  // Play around with parameters (also introducing lambda functions)
  auto frows = [](size_t M) { return 2*M; };
  auto fcols = [](size_t M) { return M; };
  auto ftriplets = [](size_t M) { return M*5; };

  // Runtime test
  Timer timer;
  for(size_t M = 2; M < 1024; M *= 2) {
    std::cout << "Runtime for " << M << "x" << M/2 << " matrix (with nnz(A) <= " << M << "):" << std::endl;
    TripletMatrix<double> A;
    CRSMatrix<double> B, E;

    A.rows = frows(M); // nrows
    A.cols = fcols(M); //ncols
    A.triplets.reserve(ftriplets(M));
    for(size_t i = 0; i < ftriplets(M); ++i) {
      A.triplets.push_back(Triplet<double>(rand() % A.rows, rand() % A.cols, rand() % 1000));
    }

    insertsort_timer.start();
    tripletToCRS_insertsort(A, B);
    insertsort_timer.stop();

    sortafter_timer.start();
    tripletToCRS_sortafter(A, E);
    sortafter_timer.stop();
    std::cout << "InsertSort took: " << insertsort_timer.duration() << " s." << std::endl;
    std::cout << "SortAfter took:  " << sortafter_timer.duration()  << " s." << std::endl;
  }
#endif
}
