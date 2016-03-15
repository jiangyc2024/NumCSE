#include <Eigen/Dense>
#include <iostream>

void corzero(double x0, Eigen::VectorXd &rates, Eigen::VectorXd &err)
{
	const unsigned int N = 15;

	double x = x0; // initial value
	Eigen::VectorXd y(N);

	for (unsigned int i=0; i<N; ++i)
	{
		x = x + (cos(x)+1)/sin(x);
		y(i) = x;
	}

	err.resize(N); rates.resize(N);	
	err = y - Eigen::VectorXd::Constant(N, x);
	rates = err.bottomRows(N-1).cwiseQuotient(err.topRows(N-1));
}
