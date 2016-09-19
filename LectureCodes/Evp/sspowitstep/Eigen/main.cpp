#include <iostream>
#include <Eigen/Dense>
#include "sspowitstep.hpp"


int main()
{
	int n = 6;
	Eigen::VectorXd v, w;
	v = Eigen::VectorXd::Random(n);
	w = Eigen::VectorXd::Random(n);

	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	A = A*A.transpose();

	sspowitstep(A,v,w);

	std::cout << "v:" << std::endl;
	std::cout << v << std::endl;
	std::cout << "w:" << std::endl;
	std::cout << w << std::endl;
}
