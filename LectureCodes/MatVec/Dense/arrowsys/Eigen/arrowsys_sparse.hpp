///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace arrowsys {


using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::VectorXi;

//NOLINTBEGIN(bugprone-easily-swappable-parameters)
/* SAM_LISTING_BEGIN_0 */
template <class solver_t>
VectorXd arrowsys_sparse(const VectorXd &d, 
                         const VectorXd &c,
                         const VectorXd &b, double alpha,
                         const VectorXd &y) {
	const Eigen::Index n = d.size();
	SparseMatrix<double> A(n+1, n+1); // default: column-major
	VectorXi reserveVec = VectorXi::Constant(n+1, 2); // nnz per col
	reserveVec(n) = static_cast<int>(n+1); // last full col
	A.reserve(reserveVec);
	for(int j = 0; j < n; ++j){	// initalize along cols for efficiency
		A.insert(j,j) = d(j);	// diagonal entries
		A.insert(n,j) = b(j);	// bottom row entries
	}
	for(int i = 0; i < n; ++i){
		A.insert(i,n) = c(i);	// last col
	}
	A.insert(n,n) = alpha;		// bottomRight entry
	A.makeCompressed();
	return solver_t(A).solve(y);
}
/* SAM_LISTING_END_0 */
//NOLINTEND(bugprone-easily-swappable-parameters)

} // namespace arrowsys