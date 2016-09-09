#pragma once

#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
//! simple Cholesky factorization
void cholfac(const MatrixXd &A, MatrixXd &R){
	int n = A.rows();
	R = A; 
	for(int k = 0; k < n; ++k){
		for(int j = k+1; j < n; ++j)
			R.row(j).tail(n-j) -= R.row(k).tail(n-j)*R(k,j)/R(k,k);
		R.row(k).tail(n-k) /= std::sqrt(R(k,k));
	}
	R.triangularView<StrictlyLower>().setZero();
}
/* SAM_LISTING_END_0 */
