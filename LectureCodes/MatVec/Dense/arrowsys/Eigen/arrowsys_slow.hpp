#pragma once

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
VectorXd arrowsys_slow(const VectorXd &d, const VectorXd &c, const VectorXd &b, const double alpha, const VectorXd &y){
	int n = d.size();
	MatrixXd A(n + 1,n + 1); A.setZero();
	A.diagonal().head(n) = d;
	A.col(n).head(n) = c;
	A.row(n).head(n) = b;
	A(n, n) = alpha;
	return A.lu().solve(y);
}
/* SAM_LISTING_END_0 */
