#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Structure holding a triplet in format (row, col, value)
 * For convenience, we provide a void constructor and a init constructor. All members are public
 * Used in the TripletMatrix class to provide a nice wrapper for a triplet
 * @tparam scalar represents the type of the triplet (e.g. double)
 */
/* SAM_LISTING_BEGIN_0 */
template <class scalar>
struct Triplet {
  // Default constructor
  Triplet()
  : i(0), j(0), v(0) { }

  // Constructor taking indexes i and j and a value v
  Triplet(std::size_t i_, std::size_t j_, scalar v_)
    : i(i_), j(j_), v(v_) { }

    std::size_t i, j; // Row and col
    scalar v; // value @ (row,col)
};
/* SAM_LISTING_END_0 */

/* @brief Defines a matrix stored in triplet format (using the Triplet<scalar> class.
 * Triplets may be duplicated and in *any* order. If there is a multiple triplet for (row,col) pair, we assume
 * that the values are intended to be added together
 * Also stores dimension
 * @tparam scalar represents the scalar type of the matrix (and of the triplet) (e.g. double)
 */
/* SAM_LISTING_BEGIN_1 */
template <class scalar>
struct TripletMatrix {
  std::size_t rows, cols; // Sizes: nrows and ncols
  std::vector<Triplet<scalar> > triplets;

  Matrix<scalar, -1, -1> densify();
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
/* @brief Structure holding a pair column index-value to be used in CRS format
 * Provides handy constructor and comparison operators.
 * @tparam scalar represents the scalar type of the value stored (e.g. double)
 */
template <class scalar>
struct ColValPair {
  ColValPair(std::size_t col_, scalar v_)
    : col(col_), v(v_) { }

  /* @brief Comparison operator for std::sort and std::lower\_bound
   * Basic sorting operator < for use with std::functions (i.e. for ordering according to first component (col)
   * We keep the column values sorted, either by sorting after insertion of by sorted insertion
   * @return true if this->col < other.col
   */
  bool operator<(const ColValPair& other) const {
    return this->col < other.col;
  }

  std::size_t col; // Col index
  scalar v; // Scalar value at col
};
/* SAM_LISTING_END_2 */

/* @brief Defines a matrix stored in CRS format (using the ColValPair<scalar> struct.
 * The row_pt contains the data, indexed by row and column position
 * Also stores dimension
 * @tparam scalar represents the scalar type of the matrix (and of the ColValPair) (e.g. double)
 */
/* SAM_LISTING_BEGIN_3 */
template <class scalar>
struct CRSMatrix {
  std::size_t rows, cols; // Size of the matrix rows, cols
  std::vector< std::vector< ColValPair<scalar> > > row_pt; // Vector containing, for each row, al vector of (col, value) pairs (CRS format)

  Matrix<scalar, -1, -1> densify();
};
/* SAM_LISTING_END_3 */

/* @brief Converts *this to an eigen (dense) function
 * Loops over all triplets and add value to zero matrix
 * WARNING: May fill in a lot of nonzeros if n,m large
 * @return Matrix of Dynamic size and scalar type
 */
/* SAM_LISTING_BEGIN_4 */
TripletMatrix::Matrix<scalar, -1, -1> densify() const {
  // Initialization
  Matrix<scalar, -1, -1> M = Matrix<scalar, -1, -1>::Zero(rows, cols);

  for(auto it = triplets.begin(); it != triplets.end(); ++it) {
    M(it->i, it->j) += it->v;
  }

  return M;
}
/* SAM_LISTING_END_4 */

/* @brief Converts *this to an eigen (dense) function
 * Loops over all rows and add value at col to zero matrix
 * WARNING: May fill in a lot of nonzeros if n,m large
 * @return Matrix of Dynamic size and scalar type
 */
/* SAM_LISTING_BEGIN_5 */
CRSMatrix::Matrix<scalar, -1, -1> densify() const {
// Initialization
  Matrix<scalar, -1, -1> M = Matrix<scalar, -1, -1>::Zero(rows, cols);

  std::size_t i = 0;
  for(auto it = row_pt.begin(); it != row_pt.end(); ++it) {
    for(auto it2 = it->begin(); it2 != it->end(); ++it2) {
      M(i, it2->col) = it2->v;
    }
    ++i;
  }

  return M;
}
/* SAM_LISTING_END_5 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format
 * No assumption is made on the triplets, may be unsorted and/or duplicated
 * In case of duplicated triplets, values are added toghether
 * The output CRS matrix may be empty (in that case T = C) or may be already filled with values (in which case C += T)
 * This version inserts the pairs ColVal already sorted in the list
 * Complexity: Loops over all triplets (lets say k) and look up over all columns of the row (say $n_i$), performing a
 * linear search and an insertion (complexity $O(n_i)$)
 * Assuming $n_i$ is bounded by $n$ small, complexity is $k*n$, otherwise $k^2$
 */
/* SAM_LISTING_BEGIN_6 */
template <class scalar>
void tripletToCRS(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
  // Copy sizes and reserve memory for rows
  C.rows = T.rows;
  C.cols = T.cols;
  C.row_pt.resize(C.rows);

  // Loop over all triplets
  for(auto triplet_it = T.triplets.begin(); triplet_it != T.triplets.end(); ++triplet_it) {
    // Store row (containing ColValPairs for this row) inside row\_pt
  	std::vector<ColValPair<scalar> >& row = C.row_pt.at(triplet_it->i); // Notice the reference!
    // Create a ColVal pair that must be inserted somewhere
    ColValPair<scalar> col_val(triplet_it->j, triplet_it->v);
    // Find the place (as iterator) where to insert col\_val, i.e. the first place where col\_val < C.row\_pt.at(i)
    // WARNING: Costly call
    // returns an iterator to a ColValPair in row
    auto lb_it = std::lower_bound(row.begin(), row.end(), col_val);
    // If lower bound has already a col (and is not end()) with some value, just add the value to the column, othervise insert it
    if(lb_it != row.end() && lb_it->col == triplet_it->j) {
      lb_it->v += triplet_it->v;
    } else {
      // WARNING: Costly call (loops over all rows in the worst-case scenario)
      row.insert(lb_it, col_val);
    }
  }
}
/* SAM_LISTING_END_6 */

/* @brief Converts a matrix given as triplet matrix to a matrix in CRS format
 * No assumption is made on the triplets, may be unsorted and/or duplicated
 * In case of duplicated triplets, values are added toghether
 * The output CRS matrix may be empty (in that case T = C) or may be already filled with values (in which case C += T)
 * This version inserts all the triplets (push_back), then sorts each row and cumsum all the duplicated col values
 * Complexity: Loops over all triplets (lets say k) and do a cheap (amortized) o(1) push back (complexity k*1). Then sorts an array
 * of rows (ith quicksort complexity k * log(k))
 */
/* SAM_LISTING_BEGIN_7 */
template <class scalar>
void tripletToCRS_sortafter(const TripletMatrix<scalar>& T, CRSMatrix<scalar>& C) {
  // Copy dimensions and reserve known space
  C.rows = T.rows;
  C.cols = T.cols;
  C.row_pt.resize(C.rows);

  // Loops over all triplets and push them at the ritgh place (cheap)
  for(auto triplets_it = T.triplets.begin(); triplets_it != T.triplets.end(); ++triplets_it) {
    ColValPair<scalar> cp(triplets_it->j, triplets_it->v);
    C.row_pt.at(triplets_it->i).push_back(cp);
  }
  // Loops over all rows , sort them (according to col) and sum duplicated values
  for(auto row_it = C.row_pt.begin(); row_it != C.row_pt.end(); ++row_it) {
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
}
/* SAM_LISTING_END_7 */

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

  T.rows = nrows;
  T.cols = ncols;
  T.triplets.reserve(ntriplets); // Always reserve space if you can

  for(std::size_t i = 0; i < ntriplets; ++i) {
    // Test unordered triplets, random and possibly repeated triplets
    T.triplets.push_back(Triplet<double>(rand() % nrows, rand() % ncols, rand() % 1000));
  }

  std::cout << "***Test conversion with random matrices***" << std::endl;
  tripletToCRS(T, C);
  std::cout << "--> Frobenius norm of T - C: " << (T.densify()-C.densify()).norm() << std::endl;

  // Alternative solution
  tripletToCRS_sortafter(T, D);
  std::cout << "--> Frobenius norm of T - D: " << (T.densify()-D.densify()).norm() << std::endl;

  // Big benefit of how we defined our functions: can add new triplets to the old ones:
  std::cout << "***Test addition with random matrices***" << std::endl;
  TripletMatrix<double> T2;
  T2.rows = nrows;
  T2.cols = ncols;
  T2.triplets.reserve(ntriplets); // Always reserve space if you can
  for(std::size_t i = 0; i < ntriplets; ++i) {
    // Test unordered triplets, random and possibly repeated triplets
    T2.triplets.push_back(Triplet<double>(rand() % nrows, rand() % ncols, rand() % 1000));
  }
  tripletToCRS(T2, C);
  tripletToCRS_sortafter(T2, D);
  std::cout << "T = " << std::endl
            << T << std::endl;
  std::cout << "T2 = " << std::endl
            << T2 << std::endl;
  std::cout << "C = " << std::endl
            << C << std::endl;
  std::cout << "D = " << std::endl
            << D << std::endl;

  // Runtime test
  bool noruntime = false; // Skip runtime measurments
  bool noaddition = true; // Do not test addition of new triplets

  if(noruntime) {
      return 0;
  }
  std::cout << "***Runtime test***" << std::endl;

  // Play around with this parameters (also introducing lambda functions)
  auto frows = [](std::size_t M) { return 2*M; };
  auto fcols = [](std::size_t M) { return M; };
  auto ftriplets = [](std::size_t M) { return M*5; };

  // Do some timings
  Timer sortafter_timer, insertsort_timer;
  for(std::size_t M = 2; M < 1024; M *= 2) {
    std::cout << "Runtime for " << M << "x" << M/2 << " matrix (with nnz(A) <= " << M << "):" << std::endl;
    TripletMatrix<double> A;
    CRSMatrix<double> B, E;

    A.rows = frows(M); // nrows
    A.cols = fcols(M); //ncols
    A.triplets.reserve(ftriplets(M));
    for(std::size_t i = 0; i < ftriplets(M); ++i) {
      A.triplets.push_back(Triplet<double>(rand() % A.rows, rand() % A.cols, rand() % 1000));
    }

    insertsort_timer.start();
    tripletToCRS(A, B);
    if(!noaddition){
      tripletToCRS(A, B);
    }
    insertsort_timer.stop();

    sortafter_timer.start();
    tripletToCRS_sortafter(A, E);
    if(!noaddition) {
      tripletToCRS_sortafter(A, E);
    }
    sortafter_timer.stop();
    std::cout << "InsertSort took: " << insertsort_timer.duration() << " s." << std::endl;
    std::cout << "SortAfter took:  " << sortafter_timer.duration()  << " s." << std::endl;
  }
}