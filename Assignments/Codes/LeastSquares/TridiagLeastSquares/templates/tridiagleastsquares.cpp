#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

/* @brief 
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd lsqEst(const VectorXd &z, const VectorXd &c)
{
    // Initialization
	int n = z.size();
    assert( z.size() == c.size() && "z and c must have same size");

	VectorXd x(2);

    // TODO: solve linear regression problem

	return x;
}
/* SAM_LISTING_END_0 */

int main() {
    // Initialization
    unsigned int n = 10;
    VectorXd z(n), c(n);
    for(size_t i=0; i<n; ++i) {
		z(i) = i+1;
		c(i) = n-i;
	}

	VectorXd x = lsqEst(z, c);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = "  << x(1) << std::endl;
}
