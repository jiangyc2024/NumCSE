#pragma once

#include <Eigen/Dense>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
VectorXd arrowsys_fast(const VectorXd &d, const VectorXd &c, const VectorXd &b, const double alpha, const VectorXd &y){
	int n = d.size();
	VectorXd z = c.array() / d.array(); 		// \Blue{$\Vz = \VD^{-1}\Vc$}  
	VectorXd w = y.head(n).array() / d.array();	// \Blue{$\Vw = \VD^{-1}\Vy_1$}
	double xi = (y(n) - b.dot(w)) / (alpha - b.dot(z));
	VectorXd x(n+1);
	x << w - xi*z, xi;
	return x;
}
/* SAM_LISTING_END_0 */
