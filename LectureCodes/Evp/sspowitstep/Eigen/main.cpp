#include <iostream>
#include <Eigen/Dense>
#include "sspowitstep.hpp"
#include "sspowitstep1.hpp"
#include "sspowitstep2.hpp"
#include "sspowitstepg.hpp"

void printVectors(const Eigen::VectorXd &v, const Eigen::VectorXd &w)
{
	std::cout << "v:" << std::endl;
	std::cout << v << std::endl;
	std::cout << "w:" << std::endl;
	std::cout << w << std::endl;
}

int main()
{
	int n = 6;
	Eigen::VectorXd v, w;
	Eigen::VectorXd v0, w0;
	v0 = Eigen::VectorXd::Random(n);
	w0 = Eigen::VectorXd::Random(n);

	v = v0;
	w = w0;
	Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
	A = A*A.transpose();

	v = v0;
	w = w0;
	sspowitstep(A,v,w);
	printVectors(v,w);

	v = v0;
	w = w0;
	sspowitstep1(v,w);
	printVectors(v,w);

	v = v0;
	w = w0;
	sspowitstep2(v,w);
	printVectors(v,w);

	Eigen::MatrixXd B = Eigen::MatrixXd::Random(n,n);
	Eigen::MatrixXd V = sspowitstepg(A, B);
}


