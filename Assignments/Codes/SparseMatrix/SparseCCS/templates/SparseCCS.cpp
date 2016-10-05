#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Compute the CCS format of matrix $A$
 * \param[in] A An $n \times n$ matrix
 * \param[out] val Vector of nonzero values of $A$
 * \param[out] row_ind Row indices of each element of 'val'
 * \param[out] col_ptr Indices of the elements in 'val' which start a column of $A$
 */
/* SAM_LISTING_BEGIN_0 */
void CCS(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)
{
	// Number of rows and columns
	m = A.rows();
	n = A.cols();

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
	col_ptr(0) = 0;

    // TODO: compute the CCS format of matrix $A$

	return A;
}
/* SAM_LISTING_END_0 */

/* @brief Compute the CCS format of matrix $A$ using Eigen methods
 * \param[in] A An $n \times n$ matrix
 * \param[out] val Vector of nonzero values of $A$
 * \param[out] row_ind Row indices of each element of 'val'
 * \param[out] col_ptr Indices of the elements in 'val' which start a column of $A$
 */
/* SAM_LISTING_BEGIN_1 */
void CCS_eigen(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)
{
	A.makeCompressed();
	val = A.valuePtr(); // Pointer to values
	row_ind = A.innerIndextr();  // Pointer to indices
	col_ptr = A.outerIndexPtr(); // Pointer to first indices of each inner vector

	return A;
}
/* SAM_LISTING_END_1 */

int main() {
    // Initialization
    unsigned int n = 6;
    MatrixXd A(n,n);
	// Poisson matrix
    A <<  4, -1,  0, -1,  0,  0,
         -1,  4, -1,  0, -1,  0,
          0, -1,  4,  0,  0, -1,
         -1,  0,  0,  4, -1,  0,
          0, -1,  0, -1,  4, -1,
		  0,  0, -1,  0, -1,  4;
	VectorXd val_1, row_ind_1, col_ptr_1;
	VectorXd val_2, row_ind_2, col_ptr_2;

    // Test 'CCS'
    CCS(A, val_1, row_ind_1, col_ptr_1);

    // Test '$CCS_eigen$'
    CCS_eigen(A, val_2, row_ind_2, col_ptr_2);

    // Verify that the solutions are the same
    std::cout << "l2-norm of the difference between val = " << (val_1 - val_2).norm() << std::endl;
	std::cout << "l2-norm of the difference between row_ind = " << (row_ind_1 - row_ind_2).norm() << std::endl;
	std::cout << "l2-norm of the difference between col_ptr = " << (col_ptr_1 - col_ptr_2).norm() << std::endl;
}
