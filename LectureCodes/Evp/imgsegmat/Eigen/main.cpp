#include "imgsegmat.hpp"
#include "imread.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <libgen.h>
#include <span>

int main(int argc, char** argv)
{
	const std::span args( argv, argc );
	const std::string path = std::string(dirname(args[0])) + "/eth.bmp"; //NOLINT(concurrency-mt-unsafe)

	auto img = imread::readBMP(path);
	const Eigen::MatrixXd P = imread::greyscale(img);

	Eigen::SparseMatrix<double> D;
	Eigen::SparseMatrix<double> A;

	auto Sfun = [](double, double) { return 1; };

	std::tie(A,D) = imgsegmat::imgsegmat(P, Sfun);

}
