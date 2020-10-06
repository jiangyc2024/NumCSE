#ifndef SPAI_HPP
#define SPAI_HPP

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

using index_t = int;

/* @brief Compute $B = \argmin_{X \in P(A)} |I-AX|_F$
 * @param[in] A An $n \times n$ matrix
 * @param[out] B The $n \times n$ matrix $= \argmin_{X \in P(A)} |I-AX|_F$
 */
/* SAM_LISTING_BEGIN_1 */
SparseMatrix<double> spai(SparseMatrix<double>& A) {
	// Size check
	assert(A.rows() == A.cols() &&
		   "Matrix must be square!");
	unsigned int N = A.rows();
	
	SparseMatrix<double> B = SparseMatrix<double>(N, N);
	
	// TODO: (4-6.b) Compute the sparse approximate inverse of A by exploiting
	// Eigen's internal CCS data structure
	// START
	
	// Needed to make sure ...Ptr functions return
	// arrays specified in CCS format
	A.makeCompressed();
	
	// Obtain pointers to data of A
	double* valPtr = A.valuePtr();
	index_t* innPtr = A.innerIndexPtr();
	index_t* outPtr = A.outerIndexPtr();
	
	// Create vector for triplets of B and reserve enough space
	std::vector<Triplet<double>> triplets;
	triplets.reserve(A.nonZeros());
	
	// Loop over each column of $A$
	for(unsigned int i = 0; i < N; ++i) {
		// Number of non zeros in current column of $A$
		index_t nnz_i = outPtr[i + 1] - outPtr[i];
		// If column is empty: skip column (matrix is not invertible)
		if(nnz_i == 0) continue;
		
		// Temporarily create a (small, dense) matrix on which normal formula
		// will be applied. We project the space $\mathcal{P}(A)$
		// onto $\mathcal{P}(a_i)$
		SparseMatrix<double> C(N, nnz_i);
		std::vector<Triplet<double>> C_triplets;
		C_triplets.reserve(nnz_i * nnz_i);
		
		// Need to build matrix $C$. To this end we remove all columns
		// from $A$ for which $a_i = 0$.
		// To this end: loop over all non-zero entries of the $i$-th column
		// This loop has length $n$
		for(unsigned int k = outPtr[i]; k < outPtr[i + 1]; ++k) {
			// Row of this non-zero entry
			index_t row_k = innPtr[k];
			// Number of non-zero entries for $row_k$-th column
			index_t nnz_k = outPtr[row_k + 1] - outPtr[row_k];
			// Loop over all non-zeros of $row_k$-th-column
			// This loop has length complexity $n$
			for(unsigned int l = 0; l < nnz_k; ++l) {
				unsigned int innIdx = outPtr[row_k] + l;
				C_triplets.emplace_back(Triplet<double>(innPtr[innIdx], k - outPtr[i], valPtr[innIdx]));
			}
		}
		C.setFromTriplets(C_triplets.begin(), C_triplets.end());
		C.makeCompressed();
		
		// Compute $C^\top C$ and solve normal equation.
		// Complexity of product: $O(n^3)$
		// Complexity of solve: $O(n^3)$
		// Size of $b$ is at most $n$.
		SparseMatrix<double> S = C.transpose() * C;
		MatrixXd M = MatrixXd(S);
		VectorXd xt = C.row(i).transpose();
		VectorXd b = M.partialPivLu().solve(xt);
		
		// Loop as length $n$.
		for(unsigned int k = 0; k < b.size(); ++k) {
			triplets.emplace_back(Triplet<double>(innPtr[outPtr[i] + k], i, b(k)));
		}
	}
	
	// Build and return SPAI preconditioner
	B.setFromTriplets(triplets.begin(), triplets.end());
	B.makeCompressed();
	// END
	
	return B;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::vector<std::pair<unsigned int, unsigned int>> testSPAIPrecCG(unsigned int L) {
	std::vector<std::pair<unsigned int, unsigned int>> cg_iterations()
	
	
}
/* SAM_LISTING_END_3 */

#endif
