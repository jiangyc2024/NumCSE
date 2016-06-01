#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include "kronmultv.hpp"

int main () {
	unsigned int m = 3, n = 2, l = 2, k = 4;
	Eigen::MatrixXd A(m,n);
	A << 1, 2, 3, 4, 5, 6;
	Eigen::MatrixXd B(l,k);
	B << 5, 6, 7, 8, 9, 10, 11, 12;
	Eigen::MatrixXd C;
	Eigen::VectorXd x(n*k);
	x << 1, 3, 8, 13, 7, 4, 42, 343;
	Eigen::VectorXd y;
	kronmultv(A,B,x,y);
	std::cout << "kron(A,B)*x = " << std::endl << y << std::endl;
	Eigen::VectorXd y_matlab(m*l);
	y_matlab << 6377,9645,12937,19573,19497,29501;
	std::cout << "error norm to matlab result" << std::endl << (y_matlab - y).norm() << std::endl;
  return 0;
}
