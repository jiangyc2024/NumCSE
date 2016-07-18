#include <Eigen/Sparse>
#include <vector>
#include <assert.h>

/**
 * @brief Equvalent to MATLAB function A = speye(m,n)
 * spdiags(m,n) creates an m-by-n sparse identity matrix
 * @param m number of rows of output
 * @param n number of columnss of output
 * @return m x n idenity SparseMatrix
 */
template <class T>
Eigen::SparseMatrix<T> speye(int m, int n = -1)
{
	if (n == -1) n = m;
	assert(m > 0 && n > 0);

	std::vector<Eigen::Triplet<T>> triplets;
	int ndiag = std::min(m,n);
	triplets.reserve(ndiag);
	for (int i=0; i<ndiag; ++i)
	{
		triplets.push_back(Eigen::Triplet<T>(i,i,1));
	}
	Eigen::SparseMatrix<T> E(m,n);
	E.setFromTriplets(triplets.begin(), triplets.end());
	E.makeCompressed();
	return E;
}
