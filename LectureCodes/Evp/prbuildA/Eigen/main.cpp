#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "prbuildA.hpp"
#include "../../utils/graphMarketMatrixLoader.hpp"
#include <libgen.h>
#include <string.h>

int main(int argc, char** argv)
{
	Eigen::MatrixXi G;
	std::string path = std::string(dirname(argv[0])) +	"/Harvard500.mtx";
	loadGraphMarketMatrix(G, path);

	Eigen::MatrixXd M = prbuildA(G, 0.15);
}
