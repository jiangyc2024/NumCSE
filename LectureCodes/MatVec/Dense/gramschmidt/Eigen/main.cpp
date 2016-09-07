#include <iostream>
#include <Eigen/Dense>

#include "gramschmidt.hpp"

int main () {
	// Ortho test
	unsigned int n = 9;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	Eigen::MatrixXd Q = gramschmidt( A );

	// Output should be idenity matrix
	std::cout << Q*Q.transpose() << std::endl;
  return 0;
}
