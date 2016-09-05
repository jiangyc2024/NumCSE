#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
using namespace std;
using namespace Eigen;
#include "smw.hpp"

int main () {
	int n = 10;
	MatrixXd A = MatrixXd::Random(n,n);
	VectorXd b = VectorXd::Random(n);
	VectorXd u = VectorXd::Random(n);
	VectorXd v = VectorXd::Random(n);
	Eigen::PartialPivLU<MatrixXd> lu(A);
	MatrixXd LU = lu.matrixLU();
	MatrixXd L = MatrixXd::Identity(n,n);
	L.triangularView<StrictlyLower>() = LU;
	MatrixXd U = LU.triangularView<Upper>();
	MatrixXd P = lu.permutationP();
	// Permutation adjustment is bypassed
	std::cout << smw(L,U,u,v,b) << std::endl << std::endl;
	std::cout << (L*U + u * v.transpose()).lu().solve(b) << std::endl;
	return 0;
}