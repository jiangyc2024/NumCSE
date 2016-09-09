#include <iostream>

#include <Eigen/Dense>

#include "gsroundoff.hpp"
int main () {
	unsigned int n = 10;
	Eigen::MatrixXd H(n,n);
	for(int i = 1; i <=n; ++i){
		for(int j = 1; j <=n; ++j){
			H(i-1,j-1) = 1./(i+j-1);
		}
	}
	//std::cout << H << std::endl;
	gsroundoff(H);
	return 0;
}
