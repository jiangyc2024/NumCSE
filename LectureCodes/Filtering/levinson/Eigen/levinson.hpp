#include <Eigen/Dense>

/* SAM_LISTING_BEGIN_0 */
void levinson(const Eigen::VectorXd& u, const Eigen::VectorXd& b,
Eigen::VectorXd& x, Eigen::VectorXd& y)
{
	size_t k = u.size() - 1;
	if(k == 0) {
		x.resize(1); y.resize(1);
		x(0) = b(0); y(0) = u(0);
		return;
	}
	
	Eigen::VectorXd xk, yk;
	levinson(u.head(k), b.head(k), xk, yk);
	
	double sigma = 1 - u.head(k).dot(yk);
	
	double t = (b(k) - u.head(k).reverse().dot(xk))/sigma;
	x = xk - t*yk.head(k).reverse();
	x.conservativeResize(x.size()+1);
	x(x.size()-1) = t;
	
	double s = (u(k) - u.head(k).reverse().dot(yk))/sigma;
	y = yk - s*yk.head(k).reverse();
	y.conservativeResize(y.size()+1);
	y(y.size()-1)= s;
}
/* SAM_LISTING_END_0 */
