#include <Eigen/Dense>
#include <iostream>

int main()
{
	const unsigned int N = 15;

	double x = 0.4; // initial value
	Eigen::VectorXd y(N);

	for (unsigned int i=0; i<N; ++i)
	{
		x = x + (cos(x)+1)/sin(x);
		y(i) = x;
	}
	
	Eigen::VectorXd err = y - Eigen::VectorXd::Constant(N, x);
	Eigen::VectorXd rate = err.bottomRows(N-1).cwiseQuotient(err.topRows(N-1));

	std::cout << rate << std::endl;

	return 0;
}
