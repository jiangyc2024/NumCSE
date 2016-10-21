//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Solve the linear regression problem (fitting a line to data)
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted
 */
VectorXd linRegr(const VectorXd &t, const VectorXd &y)
{
    // Initialization
	int n = t.size();
    assert( t.size() == y.size() && "t and y must have same size");

	VectorXd x(2);

    // TODO: solve linear regression problem

	return x;
}

/* @brief Solve the linearized exponential problem
 * for data points $(t_i,y_i)$, $i = 1,\dots,n$
 * @param[in] t An $n$-dimensional vector containing one side of input data
 * @param[in] y An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(x_1,x_2)$, intercept and slope of the line fitted to the linearized problem
 */
VectorXd expFit(const VectorXd &t, const VectorXd &y)
{
    // Initialization
	int n = t.size();
    assert( t.size() == y.size() && "t and y must have same size");

	VectorXd x(2);

    // TODO: solve linearized exponential problem

	return x;
}

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
