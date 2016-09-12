#pragma once
#include <vector>
#include <Eigen/Sparse>
#include <algorithm>
using namespace Eigen;
/* SAM_LISTING_BEGIN_0 */
/**
 * @brief Equvalent to MATLAB function A = spdiags(B,d,m,n)
 * spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the 
 * columns of B and placing them along the diagonals specified by d. see MATLAB doc
 * @param B A min(m,n)-by-p matrix, usually (but not necessarily) full, whose columns are the diagonals of A.
 * @param d A vector of length p whose integer components specify the diagonals in A.
 * @param m number of rows of output
 * @param n number of columnss of output
 * @return m x n SparseMatrix by taking the columns of B and placing them along the diagonals specified by d.
 */
template <class numeric_t> 
SparseMatrix<numeric_t> spdiags(const Matrix<numeric_t,-1,-1> &B, 
					const VectorXi &d, const int m, const int n) {					
	typedef Triplet<numeric_t> triplet_t;
	std::vector<triplet_t> triplets;
	triplets.reserve(std::min(m,n)*d.size());
	for (int k = 0; k < d.size(); ++k) {
		int diag = d(k);	// get diagonal
		int i_start = std::max(-diag, 0); // get row of 1st element
		int i_end = std::min(m, m-diag-(m-n)); // get row of last element
		int j = -std::min(0, -diag); // get col of 1st element
		int B_i; // start index i in matrix B
		if(m < n)
			B_i = std::max(-diag,0); // m < n
		else
			B_i = std::max(0,diag); // m >= n
		for(int i = i_start; i < i_end; ++i, ++j, ++B_i){
			triplets.push_back( {i, j,  B(B_i,k)} );
		}
	}
	SparseMatrix<numeric_t> A(m,n);
	A.setFromTriplets(triplets.begin(), triplets.end());
	A.makeCompressed();
	return A;
}
/* SAM_LISTING_END_0 */


/*
// Testcase
/// Matlab example 5A
VectorXi diag_no(3);
diag_no << -2,0,2;
MatrixXd B(5,3);
B << 	1,	6,	11,
		2,	7,	12,
		3,	8,	13,
		4,	9,	14,
		5,	10,	15;
std::cout << B << std::endl;
std::cout << spdiags(B, diag_no, 5,5) << std::endl;
std::cout << spdiags(B, diag_no, 5,4) << std::endl;
std::cout << spdiags(B, diag_no, 4,5) << std::endl;
*/
