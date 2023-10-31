///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>

#include "kronmultv.hpp"

int main () {
	const Eigen::Index m = 3; 
	const Eigen::Index n = 2;
	const Eigen::Index l = 2; 
	const Eigen::Index k = 4;
	Eigen::MatrixXd A(m,n);
	A << 1, 2, 3, 4, 5, 6;
	Eigen::MatrixXd B(l,k);
	B << 5, 6, 7, 8, 9, 10, 11, 12;
	Eigen::VectorXd x(n*k);
	x << 1, 3, 8, 13, 7, 4, 42, 343;
	const Eigen::VectorXd y = kronmultv(A,B,x);
	std::cout << "kron(A,B)*x = " << std::endl << y << std::endl;
	Eigen::VectorXd y_matlab(m*l);
	y_matlab << 6377,9645,12937,19573,19497,29501;
	std::cout << "error norm to matlab result" << std::endl << (y_matlab - y).norm() << std::endl;
  return 0;
}
