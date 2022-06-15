//
// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
// Contributors: tille, jgacon, dcasati
// This file is part of the NumCSE repository.
//

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cmath>
#include <iostream>

#include "sparseCCS.hpp"

int main() {
  // initialization of Poisson matrix
  constexpr unsigned int n = 6;
  Eigen::MatrixXd A(n, n);
  A << 4, -1, 0, -1, 0, 0, -1, 4, -1, 0, -1, 0, 0, -1, 4, 0, 0, -1, -1, 0, 0, 4,
      -1, 0, 0, -1, 0, -1, 4, -1, 0, 0, -1, 0, -1, 4;

  // Test 'CCS'
  Eigen::VectorXd val_1, row_ind_1, col_ptr_1;
  CCS(A, val_1, row_ind_1, col_ptr_1);

  /**
   * @brief Compute the CCS format of matrix $A$ using Eigen methods
   */
  /* SAM_LISTING_BEGIN_1 */
  double* val_2;
  int* row_ind_2;
  int* col_ptr_2;

  // Cast $A$ to 'SparseMatrix' and compress it to store the matrix in CCS
  // format
  Eigen::SparseMatrix<double> As = A.sparseView();
  As.makeCompressed();

  val_2 = As.valuePtr();           // Pointer to values
  row_ind_2 = As.innerIndexPtr();  // Pointer to indices
  col_ptr_2 =
      As.outerIndexPtr();  // Pointer to first indices of each inner vector
  /* SAM_LISTING_END_1 */

  // Verify that the solutions are the same
  // Compute l2-norm of the differences between the CCS vectors
  double diff_val = 0, diff_row_ind = 0, diff_col_ptr = 0;
  for (unsigned int i = 0; i < val_1.size(); ++i) {
    diff_val += std::pow(val_1(i) - *(val_2 + i), 2);
    diff_row_ind += std::pow(row_ind_1(i) - *(row_ind_2 + i), 2);
  }
  for (unsigned int i = 0; i < col_ptr_1.size(); ++i) {
    diff_col_ptr += std::pow(col_ptr_1(i) - *(col_ptr_2 + i), 2);
  }
  std::cout << "l2-norm of the difference between val = " << std::sqrt(diff_val)
            << std::endl;
  std::cout << "l2-norm of the difference between row_ind = "
            << std::sqrt(diff_row_ind) << std::endl;
  std::cout << "l2-norm of the difference between col_ptr = "
            << std::sqrt(diff_col_ptr) << std::endl;
}