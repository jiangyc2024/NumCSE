#ifndef ELLPACK_HPP
#define ELLPACK_HPP

////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

using Triplet_new = Triplet<double>;
using Triplets = std::vector<Triplet_new>;
using Vector = VectorXd;
using index_t = std::ptrdiff_t;

/* @brief Class representing a sparse matrix in Ellpack format
 * \param[in] triplets Vector of Eigen triplets
 * \param[in] m Number of rows of matrix $A$
 * \param[in] n Number of columns of matrix $A$
 */

/* SAM_LISTING_BEGIN_0 */
class EllpackMat {
public:
  EllpackMat(const Triplets &triplets, index_t m, index_t n);
  double operator()(index_t i, index_t j) const;
  void mvmult(const VectorXd &x, VectorXd &y) const;
  void mtvmult(const VectorXd &x, VectorXd &y) const;
private:
  std::vector<double> val; //< Vector containing values
  // corresponding to entries in 'col'
  std::vector<index_t> col; //< Vector containing column
  // indices of the entries in 'val'.
  // The position of a cell in 'val' and 'col'
  // is determined by its row number and original position in 'triplets'
  index_t maxcols; //< Number of non-empty columns
  index_t m, n;    //< Number of rows, number of columns
};
/* SAM_LISTING_END_0 */

/* @brief Retrieve value of cell $(i,j)$
 * \param[in] i Row index
 * \param[in] j Column index
 * \param[out] val Value of cell $(i,j)$
 */

/* SAM_LISTING_BEGIN_1 */
double EllpackMat::operator()(index_t i, index_t j) const {
  assert(0 <= i && i < m && 0 <= j && j < n && "Index out of bounds!");

  for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
	  if (col.at(l) == j)
		  return val.at(l);
  }
  return 0;
}
/* SAM_LISTING_END_1 */

/* @brief Constructor of class EllpackMat from vector of triplets
 * \param[in] triplets Vector of Eigen triplets
 * \param[in] m Number of rows of matrix $A$
 * \param[in] n Number of columns of matrix $A$
 */

/* SAM_LISTING_BEGIN_2 */
EllpackMat::EllpackMat(const Triplets &triplets, index_t m, index_t n)
    : m(m), n(n) {
  // TODO (3-15.a): implement the constructor for the class EllpackMat
  // START

  // END
}
/* SAM_LISTING_END_2 */

/* @brief Compute $y = Ax$ exploiting the Ellpack format
 * \param[in] x Vector of matrix-vector multiplication
 * \param[out] y Vector from $y = Ax$
 */
/* SAM_LISTING_BEGIN_3 */
void EllpackMat::mvmult(const VectorXd &x, VectorXd &y) const {
  assert(x.size() == n && "Incompatible vector x size!");
  assert(y.size() == m && "Incompatible vector y size!");
  // TODO (3-15.b) : implement the multiplication $\cob{A^{\top}*x}$ using the
  // class EllpackMat, with optimal complexity.

  // START
  
  // END
}
/* SAM_LISTING_END_3 */
/* @brief Compute $y = A_t x$ exploiting the Ellpack format, where $A$ is
 * transposed \param[in] x Vector of matrix-vector multiplication \param[out] y
 * Vector from $y = Ax$
 */
/* SAM_LISTING_BEGIN_4 */
void EllpackMat::mtvmult(const VectorXd &x, VectorXd &y) const {
  assert(x.size() == m && "Incompatible vector x size!");
  assert(y.size() == n && "Incompatible vector y size!");
  // TODO (3-15.c) : implement the multiplication $\cob{A^{\top}*x}$ using the class
  // EllpackMat, with optimal complexity.

  // START
  
  // END
}
/* SAM_LISTING_END_4 */

#endif
