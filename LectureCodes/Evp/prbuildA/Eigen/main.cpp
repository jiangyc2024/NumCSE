#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "prbuildA.hpp"
#include "../../utils/graphMarketMatrixLoader.hpp"

int main()
{
	Eigen::MatrixXi G;
	loadGraphMarketMatrix(G, "Harvard500.mtx");

	Eigen::MatrixXd M = prbuildA(G, 0.15);
}
