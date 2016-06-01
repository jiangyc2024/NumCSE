#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ios>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include "expval.hpp"

int main () {
	int n = 20;
	VectorXd x = VectorXd::LinSpaced(21,-20, n);
	MatrixXd res(n+1,4);
	for(int i = 0; i <= n; ++i){
		res(i,0) = x(i);
		res(i,1) = expval(x(i));
		res(i,2) = std::exp(x(i));
		res(i,3) = std::abs(res(i,2)-res(i,1))/res(i,2);
		// printing
		std::cout << std::fixed << std::setprecision(0) << std::setw(5) << res(i,0) << std::setprecision(10)<<
		std::setw(25) << std::scientific << res(i,1) << std::setw(25) << res(i,2) << 
		std::setw(25) << std::setprecision(15)<< std::fixed << res(i,3) << std::endl;
	}
	return 0;
}
