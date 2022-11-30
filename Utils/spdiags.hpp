#pragma once
#include <Eigen/Sparse>
#include <algorithm>
#include <vector>

/* SAM_LISTING_BEGIN_0 */
/**
 * @brief Equivalent to MATLAB function A = spdiags(B,d,m,n)
 * spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the 
 * columns of B and placing them along the diagonals specified by d. see MATLAB doc
 * @param B A min(m,n)-by-p matrix, usually (but not necessarily) full, whose columns are the diagonals of A.
 * @param d A vector of length p whose integer components specify the diagonals in A.
 * @param m number of rows of output
 * @param n number of columnss of output
 * @return m x n Eigen::SparseMatrix by taking the columns of B and placing them along the diagonals specified by d.
 */
template <class numeric_t> 
Eigen::SparseMatrix<numeric_t> spdiags(const Eigen::Matrix<numeric_t,-1,-1> &B, 
					const Eigen::VectorXi &d, const Eigen::Index m, const Eigen::Index n) {					
	typedef Eigen::Triplet<numeric_t> triplet_t;
	std::vector<triplet_t> triplets;
	triplets.reserve(std::min(m,n)*d.size());
	for (Eigen::Index k = 0; k < d.size(); ++k) {
		const int diag = d(k);	// get diagonal
		const int i_start = std::max(-diag, 0); // get row of 1st element
		const Eigen::Index i_end = std::min(m, m-diag-(m-n)); // get row of last element
		int j = -std::min(0, -diag); // get col of 1st element
		Eigen::Index B_i = m < n ? std::max(-diag,0) : std::max(0,diag); // start index i in matrix B
		for(int i = i_start; i < i_end; ++i, ++j, ++B_i){
			triplets.emplace_back(i, j, B(B_i,k));
		}
	}
	Eigen::SparseMatrix<numeric_t> A(m,n);
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
