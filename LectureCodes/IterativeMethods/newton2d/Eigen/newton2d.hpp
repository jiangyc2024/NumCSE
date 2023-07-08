#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

inline void newton2d()
{
	// F
	auto F = [](const Eigen::Vector2d &x){ 
   		Eigen::Vector2d res;
		res << pow(x(0),2) - pow(x(1),4), x(0) - pow(x(1),3);
		return res;
	};

	// jacobian of F
	auto DF= [] (const Eigen::Vector2d &x){
		Eigen::Matrix2d J;
		J << 2*x(0)	, -4*pow(x(1),3), 
		  	 1		, -3*pow(x(1),2);
	  	return J;
	};	

	Eigen::Vector2d x;
	Eigen::Vector2d x_ast;
	Eigen::Vector2d s;
	x << 0.7, 0.7; // initial value
	x_ast << 1, 1; // solution
   	const double tol=1E-10;

	std::vector<double> errors;
	errors.push_back((x-x_ast).norm());

	do
	{
		s = DF(x).lu().solve(F(x));
	   	x = x-s; // newton iteration
		errors.push_back((x-x_ast).norm());
	}
	while (s.norm() > tol*x.norm());

	// create eigen vector from std::vector
	const unsigned int n = errors.size();
	Eigen::Map<Eigen::VectorXd> err(errors.data(), n);
	
	std::cout << "Errors:" << std::endl;
	std::cout << err << std::endl;

	// compute the convergence rate of each iteration	
	const Eigen::VectorXd logDiff = err.bottomRows(n-1).array().log() - err.topRows(n-1).array().log();
	const Eigen::VectorXd rates = logDiff.bottomRows(n-2).cwiseQuotient(logDiff.topRows(n-2));
	
	std::cout << "Rates:" << std::endl;
	std::cout << rates << std::endl;
}
