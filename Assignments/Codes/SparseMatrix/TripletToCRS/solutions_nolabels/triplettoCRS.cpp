//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
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
template <typename scalar>
struct Triplet {
  // Default constructor
  Triplet()
  : i(0), j(0), v(0) { }

  // Constructor taking indexes i and j and a value v
  Triplet(size_t i_, size_t j_, scalar v_)
    : i(i_), j(j_), v(v_) { }

    size_t i, j; // Row and col
    scalar v; // Value in (row,col)
};

/* @brief Defines a matrix stored in triplet format (using the Triplet<scalar> class).
 * Triplets may be duplicated and in *any* order.
 * If there is a multiple triplet for (row,col) pair, we assume that the values are intended to be added together.
 * Dimensions are also stored to simplify the code.
 * @tparam scalar Type of the matrix and triplets (e.g. double)
 */
template <typename scalar>
struct TripletMatrix {
  size_t rows, cols; // Sizes: nrows and ncols
  std::vector<Triplet<scalar> > triplets;

  MatrixXd densify() const;
};

/* @brief Defines a matrix stored in CRS format.
 * Dimensions are also stored to simplify the code.
 * @tparam scalar Type of the matrix and CRS vectors (e.g. double)
 */
template <typename scalar>
struct CRSMatrix {
  size_t rows, cols; // Sizes: nrows and ncols
  std::vector<scalar> val;
  std::vector<size_t> col_ind;
  std::vector<size_t> row_ptr;

  MatrixXd densify() const;
};

/* @brief Convert *this (TripletMatrix) to an Eigen (dense) matrix:
 * Loop over all triplets and add values to a zero matrix.
 * WARNING: May fill in many nonzeros if large n,m
 * @return Matrix of 'scalar' type
 */
template <class scalar>
MatrixXd TripletMatrix<scalar>::densify() const {
  MatrixXd M;
  // Initialization
  M = MatrixXd::Zero(rows, cols);

  for(auto it = triplets.begin(); it != triplets.end(); ++it) {
    M(it->i, it->j) += it->v;
  }

  return M;
}

/* @brief Convert *this (CRSMatrix) to an Eigen (dense) matrix:
 * Loop over all rows and add values at col to zero matrix.
 * WARNING: May fill in many nonzeros if large n,m
 * @return Matrix of 'scalar' type
 */
template <typename scalar>
MatrixXd CRSMatrix<scalar>::densify() const {
  MatrixXd M;
// Initialization
  M = MatrixXd::Zero(rows, cols);
  std::vector<size_t> row_ptr_end = row_ptr;
  row_ptr_end.push_back(col_ind.size());
  
  for(size_t i=0; i<row_ptr.size(); ++i) {
	  for(size_t j=row_ptr_end[i]; j<row_ptr_end[i+1]; ++j) {
		  M(i, col_ind[j]) = val[j];
	  }
  }

  return M;
}

/* @brief Structure holding a pair column index-value to be used in CRS format
 * Provides handy constructor and comparison operators.
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 */
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

/* @brief Convert 'preCRS' to proper CRS format (CRSMatrix)
 * @param[in] preCRS Vector of rows of pairs (col\_ind, val).
 * The original matrix does not have rows only made by zeros.
 * @param[out] C Proper CRSMatrix format
 */
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
template <typename scalar>
void tripletToCRS_insertsort(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
  // Initialization
  C.rows = T.rows;
  C.cols = T.cols;
  // Temporary "quasi-CRS" format for the conversion
  std::vector< std::vector< ColValPair<scalar> > > preCRS;
  preCRS.resize(C.rows);

  // Loop over all triplets
  for(auto triplet = T.triplets.begin(); triplet != T.triplets.end(); ++triplet) {

    // Store row (containing all ColValPair for this row) inside row\_ptr
    std::vector<ColValPair<scalar> >& row = preCRS.at(triplet->i); // Notice the reference!
  	
  	// Create a ColValPair that must be inserted somewhere
  	ColValPair<scalar> col_val(triplet->j, triplet->v);
  	
    // Find the place (as iterator) where to insert col\_val, i.e. the first position $lb_it$ s.t col\_val < row
    // WARNING: Costly call, which returns an iterator to an element (ColValPair) in row
    auto lb_it = std::lower_bound(row.begin(), row.end(), col_val);
    
    // If lower bound has already a col (and is not end()) with some value,
    // just add the value to the column, othervise insert it
    if(lb_it != row.end() && lb_it->col == triplet->j) {
	  lb_it->v += triplet->v;
	} else {
	  // WARNING: Costly call (loop over all rows in the worst-case scenario)
	  row.insert(lb_it, col_val);
	}
  }
  
  // Convert 'preCRS' to proper CRS format (CRSMatrix)
  properCRS(preCRS, C);
  
}

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format.
 * No assumption is made on the triplets, which may be unsorted and/or duplicated,
 * but it is imposed that the original matrix does not have rows only made by zeros.
 * In case of duplicated triplets, values are added toghether.
 * The output CRS matrix may be empty (in which case T = C) or may be already filled with values (in which case C += T).
 * This version inserts all the triplets (push_back), then sorts each row and cumsum all the duplicated col values.
 * Complexity: Loop over all $k$ triplets and do a cheap (amortized) o(1) push_back (complexity $k*1$).
 * Then sorts an array of rows (i-th quicksort complexity $k * log(k)$).
 */
template <typename scalar>
void tripletToCRS_sortafter(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
  // Initialization
  C.rows = T.rows;
  C.cols = T.cols;
  // Temporary "quasi-CRS" format for the conversion
  std::vector< std::vector< ColValPair<scalar> > > preCRS;
  preCRS.resize(C.rows);

  // Loop over all triplets and push them at the right place (cheap)
  for(auto triplet = T.triplets.begin(); triplet != T.triplets.end(); ++triplet) {	
	ColValPair<scalar> cp(triplet->j, triplet->v);
	preCRS.at(triplet->i).push_back(cp);
  }
  // Loop over all rows, sort them (according to col) and sum duplicated values
  for(auto row_it = preCRS.begin(); row_it != preCRS.end(); ++row_it) {

    // WARNING: costly call in this case
    std::sort(row_it->begin(), row_it->end());
    
    // Sum duplicated cols
    // NOTE: iterating backwards should be faster
    for(auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it) {
	
      // Check that next col is not at the end and that has same col index
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
  
  // Convert 'preCRS' to proper CRS format (CRSMatrix)
  properCRS(preCRS, C);
  
}

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

  T.rows = nrows;
  T.cols = ncols;
  T.triplets.reserve(ntriplets); // Always reserve space if you can

  for(size_t i = 0; i < ntriplets; ++i) {
    // Test unordered triplets, random and possibly repeated triplets
    T.triplets.push_back(Triplet<double>(rand() % nrows, rand() % ncols, rand() % 1000));
  }

  std::cout << "***Test conversion with random matrices***" << std::endl;
  tripletToCRS_insertsort(T, C);
  tripletToCRS_sortafter( T, D);
  std::cout << "--> Frobenius norm of T - C: " << (T.densify()-C.densify()).norm() << std::endl;
  std::cout << "--> Frobenius norm of T - D: " << (T.densify()-D.densify()).norm() << std::endl;

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
  Timer insertsort_timer, sortafter_timer;
  for(size_t n = 2; n < 1024; n *= 2) {
    std::cout << "Runtime for " << n << "x" << n/2 << " matrix (with nnz(A) <= " << n << "):" << std::endl;
    TripletMatrix<double> A;
    CRSMatrix<double> E, F;

    A.rows = frows(n); // nrows
    A.cols = fcols(n); //ncols
    A.triplets.reserve(ftriplets(n));
    for(size_t i = 0; i < ftriplets(n); ++i) {
      A.triplets.push_back(Triplet<double>(rand() % A.rows, rand() % A.cols, rand() % 1000));
    }

    insertsort_timer.start();
    tripletToCRS_insertsort(A, E);
    insertsort_timer.stop();

    sortafter_timer.start();
    tripletToCRS_sortafter(A, F);
    sortafter_timer.stop();
    std::cout << "InsertSort took: " << insertsort_timer.duration() << " s." << std::endl;
    std::cout << "SortAfter took:  " << sortafter_timer.duration()  << " s." << std::endl;
  }
}
