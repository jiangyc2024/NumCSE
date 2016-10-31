#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Solve the linear regression problem (fitting a line to data)
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted
 */
/* SAM_LISTING_BEGIN_0 */
VectorXd linReg(const VectorXd &t, const VectorXd &y)
{
    // Initialization
	int n = t.size();
    assert( t.size() == y.size() && "t and y must have same size");

	VectorXd x(2);

	VectorXd b = y;
	
	VectorXd ones(n);
	 ones.setOnes(n);
	 
	MatrixXd A(n,2); // Coefficient matrix of overdetermined linear system
	A.col(0) = ones;
	A.col(1) = t;

	// Normal equations are needed for "Eigen::FullPivLu" solver
	MatrixXd lhs = A.transpose() * A; // Left-hand side
	VectorXd rhs = A.transpose() * b; // Right-hand side
	x = lhs.fullPivLu().solve(rhs); // Least accurate than methods below but fastest
	
	// Alternatives (least squares problem is immediately solved)
  //x = A.colPivHouseholderQr().solve(b); // In between
  //x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b); // Most accurate but slowest

	return x;
}
/* SAM_LISTING_END_0 */

/* @brief Solve the linearized exponential problem
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted to the linearized problem
 */
/* SAM_LISTING_BEGIN_1 */
VectorXd expFit(const VectorXd &t, const VectorXd &y)
{
    // Initialization
	int n = t.size();
    assert( t.size() == y.size() && "t and y must have same size");

	VectorXd x(2);

	// Define the correct input $y$ for 'linReg'
	VectorXd y_(n);
	for(size_t i=0; i<n; ++i) {
		y_(i) = log(y(i));
	}
	
	VectorXd x_ = linReg(t, y_);
	
    // Set the proper fitted parameters
    x(0) = exp(x_(0));
    x(1) = x_(1);

	return x;
}
/* SAM_LISTING_END_1 */

int main() {
    // Initialization
    unsigned int n = 6;
    VectorXd t(n), y(n);
    t << 0, 1, 2, 3, 4, 5;
	y << 1, 3, 7, 20, 55, 148;

	VectorXd x = expFit(t, y);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = "  << x(1) << std::endl;
}
