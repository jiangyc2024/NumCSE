#ifndef SPARSECCS_HPP
#define SPARSECCS_HPP

//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;



/* \brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Unstable GS algorithm: output is prone to cancellation issues.
 * \param[in] $\mathbf{A}$ Matrix of linearly independent columns
 * \return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
/* SAM_LISTING_BEGIN_0 */
void CCS(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)
{
  // Number of rows and columns
  int m = A.rows();
  int n = A.cols();
  //TODO (3-14.a): Implement a function that stores the matrix A 
  // in the vectors (val,row_ind,col_ptr) in CCS format
  // Do not use Eigen methods for Sparse matrices.
  
  //START
  // Number of nonzero entries
  int nnz = 0;
  for(int i=0; i<m; ++i) { // Row iterator
    for(int j=0; j<n; ++j) { // Col iterator
      if(A(i,j) != 0) {
        ++nnz;
      }
    }
  }

  // Initialization
  val.resize(nnz);
  row_ind.resize(nnz);
  col_ptr.resize(n);

  // Store $A$ in CCS format
  int index = 0;
  for(int j=0; j<n; ++j) { // Col iterator

    // Update 'col\_ptr'
    col_ptr(j) = index;

    for(int i=0; i<m; ++i) { // Row iterator
      if(A(i,j) != 0) {
        // Record the value to 'val'
        val(index) = A(i,j);
        // Record the row index to 'row\_ind'
        row_ind(index) = i;
        ++index;
      }
    }
  }
  // END
}
/* SAM_LISTING_END_0 */

#endif