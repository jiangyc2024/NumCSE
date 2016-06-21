#include <Eigen/Dense>
#include <iostream>
#include <libgen.h>
#include "imgsegmat.hpp"
#include "imread.hpp"


int main(int arc, char** argv)
{
	std::string path = std::string(dirname(argv[0])) +	"/eth.bmp";

	auto img = readBMP(path);
	Eigen::MatrixXd P = greyscale(img);

	int N = P.rows() * P.cols();
	Eigen::SparseMatrix<double> D(N,N);
	Eigen::SparseMatrix<double> A(N,N);

	auto Sfun = [](double a, double b) { return 1; };

	std::tie(A,D) = imgsegmat(P, Sfun);
}
