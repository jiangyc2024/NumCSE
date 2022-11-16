#include <Eigen/Sparse>
#include <assert.h>
#include <vector>


/**
 * @brief Equvalent to MATLAB function A = speye(m,n)
 * spdiags(m,n) creates an m-by-n sparse identity matrix
 * @param m number of rows of output
 * @param n number of columnss of output
 * @return m x n idenity SparseMatrix
 */
template <class T>
Eigen::SparseMatrix<T> speye(Eigen::Index m, Eigen::Index n = -1)
{
	if (n == -1) { 
		n = m;
	}
	assert(m > 0 && n > 0);

	std::vector<Eigen::Triplet<T>> triplets;
	const Eigen::Index ndiag = std::min(m,n);
	triplets.reserve(ndiag);
	for (Eigen::Index i=0; i<ndiag; ++i)
	{
		triplets.push_back(Eigen::Triplet<T>(i,i,1));
	}
	Eigen::SparseMatrix<T> E(m,n);
	E.setFromTriplets(triplets.begin(), triplets.end());
	E.makeCompressed();
	return E;
}
