#include <cmath>
#include <Eigen/Dense>
#include "../../sinetransform2d/Eigen/sinetransform2d.hpp"

/* SAM_LISTING_BEGIN_0 */
void fftbasedsolutionlocal(const Eigen::MatrixXd& B,
 double c, double cx, double cy, Eigen::MatrixXd& X)
{
	int m = B.rows();
	int n = B.cols();
	
	// Eigen's meshgrid
	Eigen::MatrixXd I = Eigen::RowVectorXd::LinSpaced(n,1,n).replicate(m,1);
	Eigen::MatrixXd J = Eigen::VectorXd::LinSpaced(m,1,m).replicate(1,n);
	
	// FFT
	Eigen::MatrixXd X_;
	sinetransform2d(B, X_);
	
	// Translation
	Eigen::MatrixXd T;
	T = c + 2*cx*(M_PI/(n+1)*I).array().cos() +
			2*cy*(M_PI/(m+1)*J).array().cos();
	X_ = X_.cwiseQuotient(T);
	
	sinetransform2d(X_, X);
	X = 4*X/((m+1)*(n+1));
}
/* SAM_LISTING_END_0 */
