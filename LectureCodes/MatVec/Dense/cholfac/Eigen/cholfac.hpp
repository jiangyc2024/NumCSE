///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>

#include <Eigen/Dense>

namespace cholfac {


using Eigen::MatrixXd;
using Eigen::StrictlyLower;

inline 
/* SAM_LISTING_BEGIN_0 */
//! simple Cholesky factorization
void cholfac(const MatrixXd &A, MatrixXd &R){
	const Eigen::Index n = A.rows();
	R = A; 
	for(Eigen::Index k = 0; k < n; ++k){
		for(Eigen::Index j = k+1; j < n; ++j) {
			R.row(j).tail(n-j) -= R.row(k).tail(n-j)*R(k,j)/R(k,k);
		}
		R.row(k).tail(n-k) /= std::sqrt(R(k,k));
	}
	R.triangularView<StrictlyLower>().setZero();
}
/* SAM_LISTING_END_0 */


} //namespace cholfac