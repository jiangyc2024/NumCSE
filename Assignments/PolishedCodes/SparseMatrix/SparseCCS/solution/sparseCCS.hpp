#ifndef SPARSECCS_HPP
#define SPARSECCS_HPP

//
// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
// Contributors: tille, jgacon, dcasati
// This file is part of the NumCSE repository.
//

#include <Eigen/Core>

/**
 * @brief Converts a sparse matrix into compressed column storage format.
 *
 * @param A the matrix to convert
 * @param val value vector
 * @param row_ind row index vector
 * @param col_ptr column pointer vector
 */
/* SAM_LISTING_BEGIN_0 */
void CCS(const Eigen::MatrixXd& A, Eigen::VectorXd& val,
         Eigen::VectorXd& row_ind, Eigen::VectorXd& col_ptr) {
  // number of rows and columns
  const unsigned int m = A.rows();
  const unsigned int n = A.cols();

  // TODO: (2-14.a) Implement a function that stores the matrix A
  // in the vectors (val,row_ind,col_ptr) in CCS format
  // Do not use Eigen methods for Sparse matrices.
  // START

  // number of nonzeroes in the matrix
  unsigned int nnz = 0;
  for (unsigned int i = 0; i < m; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      if (A(i, j) != 0) ++nnz;
    }
  }

  // initialization
  val.resize(nnz);
  row_ind.resize(nnz);
  col_ptr.resize(n);

  // store $A$ in CCS format
  unsigned int index = 0;
  for (unsigned int j = 0; j < n; ++j) {  // col iterator
    // update 'col\_ptr'
    col_ptr(j) = index;

    for (unsigned int i = 0; i < m; ++i) {  // row iterator
      if (A(i, j) != 0) {
        // record the value to 'val'
        val(index) = A(i, j);
        // record row index to 'row\_ind'
        row_ind(index) = i;
        ++index;
      }
    }
  }
  // END
}
/* SAM_LISTING_END_0 */

#endif