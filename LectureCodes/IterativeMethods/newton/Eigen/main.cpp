#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "newton.hpp"

void newton2Ddriver(void)
{
	// Function F defined through lambda function
	auto F = [](const Eigen::Vector2d &x) {
		Eigen::Vector2d z;
		const double x1 = x(0),x2=x(1);
		z << x1*x1-2*x1-x2+1,x1*x1+x2*x2-1;
		return(z);
	};
	// lambda function based on DF
	auto DF = [](const Eigen::Vector2d &x,const Eigen::Vector2d &f) {
		Eigen::Matrix2d J;
		const double x1 = x(0),x2=x(1);
		J << 2*x1-2,-1,2*x1,2*x2;
		Eigen::Vector2d s = J.lu().solve(f);
		return(s);
	};
	Eigen::Vector2d x; x << 2,3; // initial guess
	newton(F,DF,x,1E-6,1E-8);
	std::cout << "||F(x)|| = " << F(x).norm() << std::endl;
}

int main()
{
	newton2Ddriver();
}
