#include "prbuildA.hpp"
#include "../../utils/graphMarketMatrixLoader.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstring>
#include <iostream>
#include <libgen.h>
#include <span>

int main(int argc, char** argv)
{
	const std::span args(argv, argc);
	Eigen::MatrixXi G;
	const std::string path = std::string(dirname(args[0])) + "/Harvard500.mtx"; //NOLINT(concurrency-mt-unsafe)
	loadGraphMarketMatrix(G, path);

	const Eigen::MatrixXd M = prbuilda::prbuildA(G, 0.15);
}
